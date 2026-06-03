<div id="main" class="col-md-9" role="main">

# Parallel performance: OpenMP threads vs R-side parallelism

<div class="section level2">

## Two kinds of parallelism

`heatwave3` detects marine heatwaves one pixel at a time, and pixels are
completely independent of one another: a pixel’s climatology and events
depend only on its own temperature time series. That independence is
what makes the problem *embarrassingly parallel*, and it can be
exploited in two distinct ways.

1.  **OpenMP threads (in-process).** The C++ core spreads pixels across
    CPU cores using OpenMP. Everything happens inside a single R
    process: the data slab is read once into shared memory, the threads
    compute, and one set of output files is written. Controlled by the
    `n_threads` argument (or `setHW3threads()`). This is the fast path,
    when it is available.

2.  **R-side parallelism (across processes).** Several R processes each
    run `heatwave3` single-threaded on a *spatial tile* of the grid, and
    the results are combined. Controlled with the `parallel` package.
    This works even when OpenMP is not compiled in, and it is the only
    way to scale across more than one machine.

This vignette measures both, on the same dataset, and explains when to
reach for each. It closes with platform-by-platform setup advice (Linux,
macOS, Windows).

> **Why this is relevant on macOS.** Apple’s compiler ships without
> OpenMP, so a default macOS build of `heatwave3` runs single-threaded
> regardless of what you pass to `n_threads`. The
> [Setup](#enabling-openmp) section shows how to enable OpenMP, and the
> [R-side](#method-3-r-side-parallelism-psock) section shows how to get
> parallelism without it.

</div>

<div class="section level2">

## Is OpenMP active?

`getHW3threads()` reports the default thread count. With OpenMP enabled
it returns roughly half your cores. If it returns `1`, the package was
built without OpenMP and `n_threads` will have no effect.

<div id="cb1" class="sourceCode">

``` r
library(heatwave3)
getHW3threads()
#> [1] 9        # OpenMP active (here: 50% of 18 cores)
#> [1] 1        # OpenMP NOT compiled in -> single-threaded
```

</div>

You can also confirm at the binary level (macOS/Linux):

<div id="cb2" class="sourceCode">

``` sh
# Should list an OpenMP runtime (libomp.dylib / libgomp.so / libomp.so)
otool -L $(Rscript -e 'cat(system.file("libs/heatwave3.so", package="heatwave3"))') | grep -i omp   # macOS
ldd       $(Rscript -e 'cat(system.file("libs/heatwave3.so", package="heatwave3"))') | grep -i omp   # Linux
```

</div>

</div>

<div class="section level2">

## The dataset

The OSTIA reanalysis over the South-West Indian Ocean / South African
east coast: a 400 × 200 grid at 0.05° resolution (15–35°E, 38–28°S),
14,276 daily time steps (1982–2021), of which roughly 50,000 pixels are
ocean. The climatological baseline is 1982–2011.

<div id="cb3" class="sourceCode">

``` r
sst_file <- "/Volumes/OceanData/OSTIA_East_Coast_MHW/SWIO_Jan1982-Dec2021.nc"
var_name <- "analysed_sst"
clim_period <- c("1982-01-01", "2011-12-31")
```

</div>

All three methods produce the **same deliverable**, the full event table
as an in-memory `data.frame`, so their wall-clock times are directly
comparable.

</div>

<div class="section level2">

## Method 1: single-threaded

The baseline. `n_threads = 1` forces one core regardless of the package
default.

<div id="cb4" class="sourceCode">

``` r
library(heatwave3)

run_once <- function(nthreads, lon_range = NULL) {
  cf <- tempfile(fileext = ".nc")
  ef <- tempfile(fileext = ".nc")
  on.exit(unlink(c(cf, ef)))
  detect3(
    sst_file,
    cf,
    ef,
    climatologyPeriod = clim_period,
    var_name = var_name,
    lon_range = lon_range,
    category = TRUE,
    return_df = TRUE,
    n_threads = nthreads
  )
}

system.time(ev1 <- run_once(1L))
```

</div>

</div>

<div class="section level2">

## Method 2: OpenMP multithreading

Identical call, more threads. Because pixels are independent, OpenMP
gives each thread a share of them, with no locking and no coordination.

<div id="cb5" class="sourceCode">

``` r
system.time(ev12 <- run_once(12L))

# Results are independent of thread count:
nrow(ev1) == nrow(ev12)
#> [1] TRUE
```

</div>

Sweeping the thread count maps the scaling curve:

<div id="cb6" class="sourceCode">

``` r
for (nt in c(1, 2, 4, 8, 12, 16)) {
  t <- system.time(run_once(nt))[["elapsed"]]
  cat(sprintf("%2d threads: %5.1f s\n", nt, t))
}
```

</div>

</div>

<div class="section level2">

## Method 3: R-side parallelism (PSOCK)

When OpenMP is unavailable, or to scale across machines, split the grid
into longitude tiles and run one single-threaded `heatwave3` per worker
process. Since pixels are independent, the union of the tiles is
**identical** to a whole-domain run.

We use a **PSOCK cluster** (independent R processes) rather than
fork-based `mclapply()`. On macOS, forking after the NetCDF/HDF5 and
Apple Objective-C runtimes have initialised is unsafe (sporadic worker
crashes, and the `+[NSCharacterSet initialize] ... fork()` abort). PSOCK
launches clean processes that each open their own files, like running
several `Rscript` jobs side by side. It is also the only option on
Windows, where `fork()` does not exist.

<div id="cb7" class="sourceCode">

``` r
library(parallel)

# Read the longitude coordinate once (cheap), to cut non-overlapping tiles.
g <- heatwave3:::hw3_read_sst(
  sst_file,
  var_name,
  time_range = c("1982-01-01", "1982-01-02")
)
lon <- sort(unique(g$lon))

make_tiles <- function(n) {
  pad <- 0.45 * min(diff(lon)) # < half the grid spacing
  grp <- cut(seq_along(lon), breaks = n, labels = FALSE)
  unname(lapply(split(lon, grp), function(L) c(min(L) - pad, max(L) + pad)))
}

run_psock <- function(n_workers) {
  cl <- makeCluster(n_workers, setup_strategy = "sequential")
  on.exit(stopCluster(cl))
  clusterEvalQ(cl, suppressPackageStartupMessages(library(heatwave3)))
  clusterExport(cl, c("sst_file", "var_name", "clim_period"))
  tiles <- make_tiles(n_workers)
  parts <- parLapply(cl, tiles, function(lr) {
    cf <- tempfile(fileext = ".nc")
    ef <- tempfile(fileext = ".nc")
    on.exit(unlink(c(cf, ef)))
    detect3(
      sst_file,
      cf,
      ef,
      climatologyPeriod = clim_period,
      var_name = var_name,
      lon_range = lr,
      category = TRUE,
      return_df = TRUE,
      n_threads = 1L
    ) # 1 thread per worker!
  })
  do.call(rbind, parts)
}

system.time(evp <- run_psock(16L))
```

</div>

The `n_threads = 1L` inside each worker is essential: with OpenMP
active, each freshly-loaded `heatwave3` would otherwise default to *all*
its threads, so 16 workers × 9 threads would oversubscribe the machine
badly. See [Don’t oversubscribe](#dont-oversubscribe).

</div>

<div class="section level2">

## Results

| Method            | Workers/threads | Time (s) | Speed-up |
|:------------------|----------------:|---------:|---------:|
| OpenMP (threads)  |               1 |    109.3 |      1.0 |
| OpenMP (threads)  |               2 |     72.5 |      1.5 |
| OpenMP (threads)  |               4 |     50.3 |      2.2 |
| OpenMP (threads)  |               8 |     40.2 |      2.7 |
| OpenMP (threads)  |              12 |     37.0 |      3.0 |
| OpenMP (threads)  |              16 |     35.7 |      3.1 |
| PSOCK (R workers) |               4 |     42.5 |      2.6 |
| PSOCK (R workers) |               8 |     24.9 |      4.4 |
| PSOCK (R workers) |              16 |     17.1 |      6.4 |

Wall-clock time to produce the full event table.

![Speed-up versus number of threads or workers for OpenMP and PSOCK,
against the ideal linear
line.](parallelism_files/figure-html/results-plot-1.png)

On this run (whole domain, ca. 50,000 ocean pixels, warm file cache),
the single-threaded baseline took **109 s**. OpenMP reached **3.1×** (36
s) at 16 threads, and R-side PSOCK reached **6.4×** (17 s) at 16
workers. Every method returned the same 4,863,985 events.

<div class="section level3">

### Reading the curves

Both methods speed the analysis up substantially, and both fall short of
ideal linear scaling, but PSOCK scales further than OpenMP here, for a
reason worth understanding.

The work has three phases, namely **read** the NetCDF slab, **compute**
per pixel, and **assemble** the 4.9-million-row event `data.frame`.
OpenMP parallelises only the middle phase. The read and the assembly
happen single-threaded in the one R process, so by Amdahl’s law that
serial fraction caps the speed-up, and OpenMP plateaus near 3×
regardless of thread count.

PSOCK parallelises **all three** phases. Each worker process reads only
its own longitude slice and builds only its own sub-frame, concurrently,
so the I/O and the assembly scale with the worker count too. With the
file already in the OS cache, those concurrent reads are cheap, and
PSOCK keeps gaining past the point where OpenMP has flattened.

Two caveats keep this from being a blanket “PSOCK wins”:

-   **Endpoint.** This benchmark returns the events as an in-memory
    `data.frame`, which makes the serial assembly a large part of
    OpenMP’s tail. The usual production endpoint (write one event
    NetCDF, `return_df = FALSE`) removes that assembly, shrinking
    OpenMP’s serial fraction and narrowing the gap.
-   **Storage.** PSOCK’s edge depends on those per-worker reads being
    cheap. On a cold cache or slow/network disk, sixteen processes
    reading at once can contend, eroding or reversing the advantage.
    OpenMP reads the slab exactly once.

OpenMP also stays the *simpler* option, with one call, one output file,
no tiling or merge, and no inter-process serialisation. PSOCK earns its
overhead by working with no OpenMP at all, by isolating failures to a
single tile, and by scaling across separate machines.

</div>

</div>

<div class="section level2">

## Discussion

<div class="section level3">

### When to use which

-   **OpenMP** is the default choice on a single machine when the
    package is built with it. It has the lowest overhead, uses shared
    memory, and writes one NetCDF. Set `n_threads`, or leave it at the
    package default.
-   **R-side parallelism (PSOCK)** is for when (a) you cannot enable
    OpenMP, such as on a managed Mac, (b) you want to spread work across
    **multiple machines** in a cluster, or (c) you want process
    isolation so one failing tile cannot bring down the whole run. Each
    worker writes its own tile, and you merge afterwards.
-   The two combine across a cluster. One PSOCK (or MPI) worker per
    node, each using OpenMP across that node’s cores, is the standard
    HPC pattern. On a single workstation, pick one.

</div>

<div class="section level3">

### Correctness

Because pixels are independent, every method returns the **identical**
event table (we assert `nrow()` equality above, and the per-pixel
metrics match to the bit). Tiling changes only *which* process computes
a pixel, never the result. The one exception is `detect_blob3()`, which
finds spatially **connected** heatwave blobs across pixels. That
computation must see the whole grid and **cannot** be tiled this way.

</div>

<div class="section level3">

### Do not oversubscribe

Running *N* worker processes that each launch *M* OpenMP threads creates
*N × M* threads. If that exceeds your physical cores, the threads
contend for them and everything slows down. Practical limits:

-   Pure OpenMP: `n_threads` ≤ physical cores.
-   Pure PSOCK: `n_workers` ≤ physical cores, **`n_threads = 1` in every
    worker**.
-   Hybrid (cluster): `n_workers_per_node × n_threads ≈ cores_per_node`.

`heatwave3` defaults to ca. 50% of cores so that a careless nested call
does not saturate the machine, but explicit is better than implicit.

</div>

</div>

<div class="section level2">

## Enabling OpenMP, by platform

`heatwave3`’s `configure` script probes for a working OpenMP toolchain
at install time (compile → link → load → run a parallel reduction) and
falls back to single-threaded if none is found. The install log prints
which path it took:

    Detecting OpenMP...
      Trying Apple Clang + R's bundled libomp
      -> works
    Using OpenMP:
      CXXFLAGS: -Xclang -fopenmp -I/opt/homebrew/opt/libomp/include
      LIBS:     /Library/Frameworks/R.framework/Resources/lib/libomp.dylib

<div class="section level3">

### Linux: works without extra setup

GCC and Clang support OpenMP natively and R already advertises the
flags. Install NetCDF and a compiler:

<div id="cb9" class="sourceCode">

``` sh
sudo apt install libnetcdf-dev build-essential   # Debian/Ubuntu
sudo dnf install netcdf-devel gcc-c++             # Fedora/RHEL
R CMD INSTALL heatwave3
```

</div>

OpenMP scaling is typically the best of the three platforms (low thread
creation overhead). Fork-based `mclapply()` is also safe here, though
PSOCK remains the portable choice.

</div>

<div class="section level3">

### macOS: provide an omp.h header

The R 4.x framework already ships the OpenMP runtime (`libomp.dylib`),
and `heatwave3` links *that* one so the whole R process shares a single
OpenMP runtime. Apple Clang only lacks the `omp.h` **header**, which
`configure` finds from Homebrew (or `/usr/local`). The one-time setup is
then:

<div id="cb10" class="sourceCode">

``` sh
brew install libomp netcdf      # libomp here supplies omp.h; netcdf is required
R CMD INSTALL heatwave3          # configure links R's own libomp, header from brew
```

</div>

R’s bundled `libomp.dylib` is older and lacks `__kmpc_dispatch_deinit`,
the symbol Apple Clang emits for `schedule(dynamic)`. Linking a
*second*, newer libomp (such as Homebrew’s) appears to work until
another package, or your environment’s library search path, redirects
`libomp.dylib` back to R’s copy at load time, and the package then fails
with `symbol not found ... ___kmpc_dispatch_deinit`. `heatwave3` avoids
this by using `schedule(static, 1)`, which needs only symbols R’s own
libomp already provides, so the process holds one OpenMP runtime and no
conflict arises.

For **R-side** parallelism on macOS, prefer PSOCK. If you must use
`mclapply()`, export `OBJC_DISABLE_INITIALIZE_FORK_SAFETY=YES` before
starting R, and treat concurrent forked NetCDF/HDF5 access as unsafe.

</div>

<div class="section level3">

### Windows: OpenMP works, no fork

Rtools (GCC/MinGW) supports OpenMP without extra setup and
`configure.win` locates NetCDF, so multithreading via `n_threads` works
with no extra steps. There is no `fork()` on Windows, so `mclapply()`
silently runs serially, and R-side parallelism **must** use a PSOCK
cluster (`makeCluster()` / `parLapply()`), exactly as shown above.

</div>

</div>

<div class="section level2">

## Summary

| Approach        | Enable with             | Best when                   | Notes                                         |
|-----------------|-------------------------|-----------------------------|-----------------------------------------------|
| Single-threaded | `n_threads = 1`         | debugging, tiny grids       | reference result                              |
| OpenMP          | `n_threads = N`         | one machine, OpenMP built   | low overhead, shared memory, single output    |
| PSOCK           | `parallel::makeCluster` | no OpenMP, or many machines | tile by longitude, `n_threads = 1` per worker |

All three give the same science. OpenMP is the efficient default on a
single host, and R-side parallelism is the portable fallback and the
route to multi-machine scaling.

</div>

</div>
