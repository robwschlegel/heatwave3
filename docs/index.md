# heatwave3

**`heatwave3`** detects marine heatwaves (MHWs) and cold-spells (MCS)
directly on gridded NetCDF files using the Hobday et al. (2016, 2018)
definition. It is a ground-up C++ reimplementation of
[**`heatwaveR`**](https://robwschlegel.github.io/heatwaveR/), designed
for gridded datasets rather than individual time series.

## Why heatwave3?

Working with large gridded ocean temperature datasets (such as OISST,
OSTIA, GLORYS, CMIP6) using **`heatwaveR`** requires looping over each
pixel individually in R. For a 400 × 200 grid with 40 years of daily
data, this takes over an hour. **`heatwave3`** addresses this with:

1.  **C++17 implementation** of all core algorithms (climatology, event
    detection, categorisation, block averages)
2.  **`std::thread` parallelism** across pixels, since each pixel is
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
  run-length encoding, gap joining, and 19 event metrics. MCS support.
  Optional inline Hobday et al. (2018) categories.
- **All-in-one pipeline** (`detect3`). Climatology, detection, and
  optional categories in a single call from one `name` stem; read
  results back with
  [`hw3_export()`](https://robwschlegel.github.io/heatwave3/index.html/reference/hw3_export.md).
- **Spatial blob detection** (`detect_blob3`). 4D connected-component
  labelling with voxel-level footprint output.
- **Event categorisation** (`category3`). Hobday et al. (2018)
  Moderate/Strong/Severe/Extreme categories, for both MHW and MCS.
- **Yearly aggregation** (`block_average3`) and **static threshold
  exceedance** (`exceedance3`).
- **Plotting**. `event_line3`, `geom_flame3`, `geom_lolli3`, and
  `plot_metric3` (C++-backed spatial aggregation).
- **Flexible input**. Single multi-timestep NetCDF file, directory of
  daily files, or explicit file vector.
- **Automatic dimension detection**. Identifies
  lon/lat/time/depth/temperature from CF attributes, standard names, and
  units (works with GHRSST, OISST, OSTIA, ERA5, CMIP6, NEMO, GLORYS).
- **Progress reporting**. Percentage-complete updates during long
  climatology and detection runs.
- **NetCDF-native output, exported on demand**. The compute functions
  always write gridded NetCDF;
  [`hw3_export()`](https://robwschlegel.github.io/heatwave3/index.html/reference/hw3_export.md)
  reads any product back into R or converts it to CSV, RDS, or Parquet,
  with optional variable and lon/lat/time/depth subsetting.
- **Numerical equivalence**. Climatology and event metrics match
  **`heatwaveR`** to rounding precision, validated pixel-by-pixel.

## Installation

### Development version (recommended)

The `dev` branch contains all of the most cutting-edge things being
developed:

``` r
# install.packages("remotes")
remotes::install_github("robwschlegel/heatwave3@dev")
```

If however you feel like playing it safe, the main branch is the way to
go:

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

### In-process parallelism

**`heatwave3`** uses C++ `std::thread` for in-process parallelism, not
OpenMP. There is no OpenMP toolchain to install or enable – no `libomp`,
no `omp.h`, nothing for `configure` to probe for – since `std::thread`
is part of the C++ standard library and works out of the box on every
platform. Earlier versions used OpenMP; it was replaced because a
`dlopen`-ed `libomp` could crash on macOS when heatwave3 ran alongside
other native-heavy packages (terra, sf, GDAL) or inside IDEs that embed
R (e.g. Positron). See the [Parallel
performance](https://robwschlegel.github.io/heatwave3/articles/parallelism.html)
vignette for the full story and a performance comparison.

### Thread management

**`heatwave3`** defaults to **50% of available cores** (overridable via
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

**`heatwave3`** has no OpenMP dependency, so it can never change thread
counts for other OpenMP-using packages in the same session. It is also
fork-safe under
[`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html): each
parallel region joins its worker threads before returning, so no
heatwave3 threads are ever alive at fork time.

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
**`heatwaveR`**:

| Function                                                                                              | Purpose                                                                  |
|-------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------|
| [`detect3()`](https://robwschlegel.github.io/heatwave3/index.html/reference/detect3.md)               | **Primary entry point**. Climatology, detection, and optional categories |
| [`ts2clm3()`](https://robwschlegel.github.io/heatwave3/index.html/reference/ts2clm3.md)               | Compute climatology (NetCDF → NetCDF)                                    |
| [`detect_event3()`](https://robwschlegel.github.io/heatwave3/index.html/reference/detect_event3.md)   | Detect per-pixel events (with optional inline categories)                |
| [`detect_blob3()`](https://robwschlegel.github.io/heatwave3/index.html/reference/detect_blob3.md)     | 3D spatial blob detection                                                |
| [`category3()`](https://robwschlegel.github.io/heatwave3/index.html/reference/category3.md)           | Hobday et al. (2018) event categories (reads or computes)                |
| [`block_average3()`](https://robwschlegel.github.io/heatwave3/index.html/reference/block_average3.md) | Yearly aggregation of event metrics                                      |
| [`exceedance3()`](https://robwschlegel.github.io/heatwave3/index.html/reference/exceedance3.md)       | Static threshold exceedance                                              |
| [`event_line3()`](https://robwschlegel.github.io/heatwave3/index.html/reference/event_line3.md)       | Per-pixel time series plot                                               |
| [`geom_flame3()`](https://robwschlegel.github.io/heatwave3/index.html/reference/geom_flame3.md)       | ggplot2 flame polygon geom                                               |
| [`geom_lolli3()`](https://robwschlegel.github.io/heatwave3/index.html/reference/geom_lolli3.md)       | ggplot2 lollipop geom                                                    |
| [`plot_metric3()`](https://robwschlegel.github.io/heatwave3/index.html/reference/plot_metric3.md)     | Spatial map of event metrics (C++-backed aggregation)                    |
| [`hw3_export()`](https://robwschlegel.github.io/heatwave3/index.html/reference/hw3_export.md)         | Read any product into a data.frame, or export to CSV/RDS/Parquet         |

## Vignettes

- [**Getting
  started**](https://robwschlegel.github.io/heatwave3/articles/heatwave3.html).
  Full pipeline walkthrough with spatial blob figures.
- [**NetCDF output
  internals**](https://robwschlegel.github.io/heatwave3/articles/netcdf-output.html).
  Output file structure, CF compliance, and reading in R/Python/CDO.
- [**Performance
  benchmark**](https://robwschlegel.github.io/heatwave3/articles/benchmark.html).
  **`heatwaveR`** versus **`heatwave3`** timing comparison.
- [**Parallel
  performance**](https://robwschlegel.github.io/heatwave3/articles/parallelism.html).
  In-process `std::thread` parallelism versus R-side (PSOCK)
  parallelism, with per-platform setup and a comparison against the
  retired OpenMP backend.

## Citation

If you use **`heatwave3`** in published research, please cite:

- Hobday, A.J., et al. (2016). A hierarchical approach to defining
  marine heatwaves. *Progress in Oceanography*, 141, 227–238.
  [https://doi.org/10.1016/j.pocean.2015.12.014](https://doi.org/10.1016/j.pocean.2015.12.014)
- Hobday, A.J., et al. (2018). Categorizing and naming marine heatwaves.
  *Oceanography*, 31(2), 162–173.
  [https://doi.org/10.5670/oceanog.2018.205](https://doi.org/10.5670/oceanog.2018.205)
- Schlegel, R. W., & Smit, A. J. (2018). heatwaveR: A central algorithm
  for the detection of heatwaves and cold-spells. *Journal of Open
  Source Software*, 3(27), 821.
  [https://doi.org/10.21105/joss.00821](https://doi.org/10.21105/joss.00821)

## Code of Conduct

Please note that the **`heatwave3`** project is released with a
[Contributor Code of
Conduct](https://contributor-covenant.org/version/2/1/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
