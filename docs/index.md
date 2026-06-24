# heatwave3

**heatwave3** detects marine heatwaves (MHWs) and cold-spells directly
on gridded NetCDF data using the Hobday et al. (2016, 2018) definition.
It is a ground-up C++ reimplementation of
[**heatwaveR**](https://robwschlegel.github.io/heatwaveR/), designed for
gridded datasets rather than individual time series.

## Why heatwave3?

Working with large gridded SST datasets (such as OSTIA, OISST, ERA5,
CMIP6) using heatwaveR requires looping over each pixel individually in
R. For a 400 × 200 grid with 40 years of daily data, this takes over an
hour. **heatwave3** addresses this with:

1.  **C++17 implementation** of all core algorithms (climatology, event
    detection, categorisation, block averages)
2.  **OpenMP parallelism** across pixels, since each pixel is
    independent
3.  **Direct NetCDF I/O** via libnetcdf, with no R NetCDF packages
    needed at runtime

### Performance

Benchmarked on the OSTIA South-West Indian Ocean reanalysis (400 × 200
grid, ca. 50,000 ocean pixels, 14,276 daily time steps), Apple M5 Max:

| Method                 | Time       | Speedup |
|------------------------|------------|---------|
| heatwaveR (serial)     | ca. 69 min | 1×      |
| heatwave3 (1 thread)   | 4.3 min    | **16×** |
| heatwave3 (12 threads) | 1.9 min    | **36×** |

For larger grids using daily files (Benguela region, 260 × 360 pixels,
16,049 daily OSTIA files, 12 threads): **3.4 minutes** end-to-end.

## Features

- **Climatology** (`ts2clm3`). Sliding-window percentile with Type-7
  quantile, circular-padded rolling mean smoothing, and optional linear
  detrending (Jacox et al. 2020).
- **Event detection** (`detect_event3`). Threshold exceedance,
  run-length encoding, gap joining, and 19 event metrics. Cold-spell
  support. Optional inline Hobday et al. (2018) severity categories.
- **All-in-one pipeline** (`detect3`). Climatology, detection, and
  optional categories in a single call from one `name` stem; read
  results back with
  [`hw3_export()`](https://robwschlegel.github.io/heatwave3/index.html/reference/hw3_export.md).
- **Spatial blob detection** (`detect_blob3`). 3D connected-component
  labelling with voxel-level footprint output.
- **Event categorisation** (`category3`). Hobday et al. (2018)
  Moderate/Strong/Severe/Extreme categories, for both heatwaves and
  cold-spells.
- **Yearly aggregation** (`block_average3`) and **static threshold
  exceedance** (`exceedance3`).
- **Plotting**. `event_line3`, `geom_flame3`, `geom_lolli3`, and
  `plot_metric3` (C++-backed spatial aggregation).
- **Flexible input**. Single multi-timestep NetCDF, directory of daily
  files, or explicit file vector.
- **Automatic dimension detection**. Identifies lon/lat/time/SST from CF
  attributes, standard names, and units (works with GHRSST, OISST,
  OSTIA, ERA5, CMIP6, NEMO).
- **Progress reporting**. Percentage-complete updates during long
  climatology and detection runs.
- **NetCDF-native output, exported on demand**. The compute functions
  always write gridded NetCDF;
  [`hw3_export()`](https://robwschlegel.github.io/heatwave3/index.html/reference/hw3_export.md)
  reads any product back into R or converts it to CSV, RDS, or Parquet,
  with optional variable and lon/lat/time subsetting.
- **Numerical equivalence**. Climatology and event metrics match
  heatwaveR to rounding precision, validated pixel-by-pixel.

## Installation

### Development version (recommended)

The `dev` branch contains the new C++ implementation:

``` r

# install.packages("remotes")
remotes::install_github("robwschlegel/heatwave3")
```

### System requirements

heatwave3 requires the **netCDF C library** (version 4.0+). The
`configure` script finds it automatically via `nc-config` or
`pkg-config`.

**macOS:**

``` bash
brew install netcdf
```

**Ubuntu / Debian:**

``` bash
sudo apt install libnetcdf-dev
```

**Fedora / RHEL:**

``` bash
sudo dnf install netcdf-devel
```

### OpenMP parallelism

On Linux, OpenMP works without extra setup. On macOS, the `configure`
script probes for a working OpenMP runtime (Apple Clang with R’s own
bundled libomp, Homebrew libomp, or user-supplied flags) using a
compile-link-run test. If none is found, heatwave3 falls back to
single-threaded mode.

To install Homebrew’s libomp (may help on some macOS configurations):

``` bash
brew install libomp
```

### Thread management

heatwave3 defaults to **50% of available cores** (overridable via
`R_HEATWAVE3_NUM_THREADS`). Control threads at the session level or per
function call:

``` r

library(heatwave3)
getHW3threads()       # check current default
setHW3threads(8)      # use 8 threads for all subsequent calls
setHW3threads(0)      # reset to default (50% of cores)

# Or override per call:
detect_event3(..., n_threads = 12)
```

heatwave3 never calls `omp_set_num_threads()`, so it does not change
thread counts for other OpenMP-using packages. It is also safe under
[`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html)
(fork-safe via `pthread_atfork`).

## Quick start

``` r

library(heatwave3)

sst_file <- "path/to/sst.nc"  # or a directory of daily files

# All-in-one: climatology + event detection + categories. A single 'name' stem
# writes benguela_clim.nc and benguela_events.nc.
detect3(
  file_in           = sst_file,
  name              = "benguela",
  climatologyPeriod = c("1991-01-01", "2020-12-31"),
  lon_range         = c(15, 35),
  lat_range         = c(-38, -28),
  category          = TRUE,
  hemisphere        = "south",
  n_threads         = 4
)

# Read the events back as a data.frame (or a quick preview with n =)
events <- hw3_export("benguela_events.nc")
head(events)
table(events$category)
```

### Cold-spell detection

``` r

# Use pctile = 10 for the cold tail
detect3(
  file_in           = sst_file,
  name              = "benguela_cold",
  climatologyPeriod = c("1991-01-01", "2020-12-31"),
  pctile            = 10,
  coldSpells        = TRUE,
  category          = TRUE,
  hemisphere        = "south",
  n_threads         = 12
)

table(hw3_export("benguela_cold_events.nc")$category)
```

### Per-pixel time series plot

``` r

event_line3(sst_file, "benguela_clim.nc", lon = 25.0, lat = -34.0,
            start_date = "2018-01-01", end_date = "2019-12-31")
```

### Spatial map of peak intensity

``` r

# C++-backed per-pixel aggregation, efficient even for millions of events
plot_metric3("benguela_events.nc", metric = "intensity_max", summary = "mean")
```

### Event categories (standalone)

``` r

# Reads pre-computed categories from event file (no clim_file needed)
cats <- category3("benguela_events.nc")
table(cats$category)

# Or compute from scratch for older event files
cats <- category3("benguela_events.nc", "benguela_clim.nc", hemisphere = "south")
```

### Using daily files

``` r

# Pass a directory. All .nc/.nc4 files are read, sorted, and merged
ts2clm3(
  file_in = "/path/to/daily_ostia/",
  name = "ostia",
  climatologyPeriod = c("1991-01-01", "2020-12-31"),
  n_threads = 12
)
```

### Spatial blob detection

``` r

blobs <- detect_blob3(
  sst_file  = sst_file,
  clim_file = "benguela_clim.nc",
  minVoxels = 200,
  topN      = 6,
  return    = c("event", "daily", "voxel")
)

# blobs$event: summary per blob (duration, peak area, cumulative intensity)
# blobs$daily: daily progression (area, centroid, bounding box)
# blobs$voxel: full 3D footprint (for spatial maps)
```

### Detrended climatology

``` r

# Remove linear warming trend before computing climatology
# (Jacox et al. 2020 approach)
ts2clm3(sst_file, name = "benguela_detrended",
        climatologyPeriod = c("1991-01-01", "2020-12-31"),
        detrend = TRUE)
```

## API overview

All functions are suffixed with `3` to avoid namespace conflicts with
heatwaveR:

| Function | Purpose |
|----|----|
| [`detect3()`](https://robwschlegel.github.io/heatwave3/index.html/reference/detect3.md) | **Primary entry point**. Climatology, detection, and optional categories |
| [`ts2clm3()`](https://robwschlegel.github.io/heatwave3/index.html/reference/ts2clm3.md) | Compute climatology (NetCDF → NetCDF) |
| [`detect_event3()`](https://robwschlegel.github.io/heatwave3/index.html/reference/detect_event3.md) | Detect per-pixel events (with optional inline categories) |
| [`detect_blob3()`](https://robwschlegel.github.io/heatwave3/index.html/reference/detect_blob3.md) | 3D spatial blob detection |
| [`category3()`](https://robwschlegel.github.io/heatwave3/index.html/reference/category3.md) | Hobday et al. (2018) event categories (reads or computes) |
| [`block_average3()`](https://robwschlegel.github.io/heatwave3/index.html/reference/block_average3.md) | Yearly aggregation of event metrics |
| [`exceedance3()`](https://robwschlegel.github.io/heatwave3/index.html/reference/exceedance3.md) | Static threshold exceedance |
| [`event_line3()`](https://robwschlegel.github.io/heatwave3/index.html/reference/event_line3.md) | Per-pixel time series plot |
| [`geom_flame3()`](https://robwschlegel.github.io/heatwave3/index.html/reference/geom_flame3.md) | ggplot2 flame polygon geom |
| [`geom_lolli3()`](https://robwschlegel.github.io/heatwave3/index.html/reference/geom_lolli3.md) | ggplot2 lollipop geom |
| [`plot_metric3()`](https://robwschlegel.github.io/heatwave3/index.html/reference/plot_metric3.md) | Spatial map of event metrics (C++-backed aggregation) |
| [`hw3_export()`](https://robwschlegel.github.io/heatwave3/index.html/reference/hw3_export.md) | Read any product into a data.frame, or export to CSV/RDS/Parquet |

## Vignettes

- **Getting started**. Full pipeline walkthrough with spatial blob
  figures.
- **NetCDF output internals**. Output file structure, CF compliance, and
  reading in R/Python/CDO.
- **Performance benchmark**. heatwaveR versus heatwave3 timing
  comparison.
- **Parallel performance**. OpenMP threads versus R-side parallelism,
  with per-platform setup.

## Citation

If you use heatwave3 in published research, please cite both:

- Hobday, A.J., et al. (2016). A hierarchical approach to defining
  marine heatwaves. *Progress in Oceanography*, 141, 227–238.
- Hobday, A.J., et al. (2018). Categorizing and naming marine heatwaves.
  *Oceanography*, 31(2), 162–173.

## Code of Conduct

Please note that the heatwave3 project is released with a [Contributor
Code of
Conduct](https://contributor-covenant.org/version/2/1/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
