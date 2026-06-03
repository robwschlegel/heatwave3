<div id="main" class="col-md-9" role="main">

# CLAUDE.md

<div id="claudemd" class="section level1">

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

<div class="section level2">

## What Is This Project?

`heatwave3` is an R package that detects marine heatwaves (MHWs) and
cold-spells directly on gridded NetCDF data using the Hobday et
al. (2016, 2018) definition. All core algorithms (climatology, event
detection, categorisation, block averages) are implemented in C++ with
OpenMP parallelism and direct libnetcdf I/O — no R NetCDF packages
needed at runtime.

This is a ground-up C++ reimplementation of
[heatwaveR](https://robwschlegel.github.io/heatwaveR/), designed for
gridded data rather than per-pixel time series.

</div>

<div class="section level2">

## Build & Development Commands

<div id="cb1" class="sourceCode">

``` bash
# Load package in dev mode
devtools::load_all()

# Regenerate RcppExports and documentation
Rcpp::compileAttributes()
devtools::document()

# Run tests
devtools::test()

# Run a single test file
testthat::test_file("tests/testthat/test-ts2clm3.R")

# Full R CMD check
devtools::check()
```

</div>

</div>

<div class="section level2">

## Architecture

<div class="section level3">

### C++ Core (`src/`)

All compute-intensive work is in C++ with Rcpp bindings:

-   **`netcdf_io.cpp/h`** — Direct NetCDF-C I/O with spatial/temporal
    subsetting, CF time parsing, scale\_factor/add\_offset/\_FillValue
    handling. Supports both NC\_CHAR and NC\_STRING attributes. Reorders
    data to pixel-major layout \[pixel\]\[time\] for cache-friendly
    per-pixel processing. Reads and writes climatology, event, and SST
    NetCDF files.

-   **`climatology.cpp/h`** — Full ts2clm algorithm: gap filling, NA
    interpolation, DOY assignment (non-leap years skip DOY 60), spread
    to DOY×year matrix, sliding-window percentile (Type-7 quantile
    matching R exactly), circular-padded rolling mean smoothing. OpenMP
    parallelism over pixels with progress reporting.

-   **`event_detect.cpp/h`** — Full detect\_event algorithm: threshold
    exceedance, run-length encoding, duration filtering, gap joining, 19
    event metric computations. Cold spell support (threshold inversion).
    Optional inline Hobday et al. (2018) category computation (I
    Moderate through IV Extreme) and season assignment. OpenMP
    parallelism over pixels with progress reporting.

-   **`blob_label.cpp`** — 3D connected-component labeling (union-find,
    6-connectivity, optional dateline wrap). Ported from heatwaveR.

-   **`blob_metrics.cpp`** — Per-(blob, day) spatial metric reduction.
    Ported from heatwaveR.

-   **`heatwave3_init.cpp`** — Rcpp exports bridging C++ to R. Key
    exports: `hw3_compute_clim`, `hw3_detect_events` (with
    category/hemisphere params), `hw3_read_sst`, `hw3_read_event_nc`,
    `hw3_read_clim_nc`, `hw3_category`, `hw3_block_average`,
    `hw3_read_metric_summary`.

-   **`heatwave3_types.h`** — Core structs: `EventResult` (19 metrics +
    category + season), `EventData` (for reading event NetCDF),
    `SubsetSpec`, `ClimData`, `GridData`.

</div>

<div class="section level3">

### R API (`R/`)

Function names are suffixed with `3` to avoid conflicts with heatwaveR:

-   **`detect3()`** — **Primary entry point.** All-in-one: climatology +
    detection + optional inline categories. Supports `return_df = TRUE`
    for interactive use.
-   **`ts2clm3()`** — NetCDF → climatology NetCDF (seas + thresh per DOY
    per pixel)
-   **`detect_event3()`** — SST + climatology → event NetCDF. Supports
    `category = TRUE` for inline categorisation and
    `hemisphere = "south"/"north"` for season naming. `return_df = TRUE`
    returns a data.frame directly.
-   **`category3()`** — Reads pre-computed categories from event file,
    or computes from climatology file for older events. `hemisphere`
    parameter replaces deprecated `S`.
-   **`block_average3()`** — Yearly event metric aggregation (pure C++)
-   **`exceedance3()`** — Static threshold exceedance
-   **`detect_blob3()`** — 3D spatially-connected blob detection
-   **`event_line3()`** — Per-pixel time series plot from NetCDF (C++
    reader, no ncdf4)
-   **`geom_flame3()`** / **`geom_lolli3()`** — ggplot2 custom geoms
-   **`plot_metric3()`** — Spatial map of event metrics (C++-backed
    aggregation via `hw3_read_metric_summary`)

</div>

<div class="section level3">

### Data Flow

    NetCDF (SST) → ts2clm3() → NetCDF (climatology)
                                    ↓
    NetCDF (SST) + clim → detect_event3(category=TRUE) → NetCDF (events + categories)
                                                               ↓
                                            category3() / block_average3() / plotting

Or all-in-one:

    NetCDF (SST) → detect3(category=TRUE, return_df=TRUE) → data.frame

</div>

<div class="section level3">

### Key Design Decisions

1.  **Pixel-major memory layout**: SST data reordered to
    \[pixel\]\[time\] for contiguous per-pixel time series access during
    OpenMP parallel loops.
2.  **Type-7 quantile**: Matches R’s default `quantile()` exactly via
    linear interpolation.
3.  **CF ragged array**: Event NetCDF uses an `event` dimension with
    lon/lat coordinate variables — avoids wasteful padding to
    max\_events.
4.  **No R NetCDF dependency at runtime**: Uses libnetcdf C API directly
    via `configure` script. All four functions that previously used
    ncdf4 (`category3`, `block_average3`, `plot_metric3`, `event_line3`)
    now use the C++ reader.
5.  **Inline categories**: Category computation uses the same
    climatology data already loaded during detection — zero extra I/O.

</div>

</div>

<div class="section level2">

## Build System

-   `configure` — Finds libnetcdf via nc-config, pkg-config, or common
    paths. Adds rpath on macOS.
-   `src/Makevars.in` — Template substituted by configure.
-   `src/Makevars` — Generated; links `-lnetcdf` + OpenMP flags.
-   SystemRequirements: netcdf (&gt;= 4.0), C++17

</div>

<div class="section level2">

## Testing

Tests compare heatwave3 output pixel-by-pixel against heatwaveR on real
OSTIA SST data from
`/Volumes/OceanData/OSTIA_East_Coast_MHW/SWIO_Jan1982-Dec2021.nc`. Tests
verify: - Climatology seas/thresh match to 1e-4 (rounding precision) -
Event count, duration exact match - Event intensity metrics match to
2e-4

</div>

<div class="section level2">

## Performance

Benchmark on 400 pixels (20×20, 14K time steps): - heatwave3: 0.53 sec
(4 threads) — **\~71× faster** than heatwaveR

Benguela analysis (93,600 pixels, 16,049 daily files, 12 threads): -
Climatology: \~1 min - Detection + categories: \~3 min - Total: \~4 min
end-to-end

</div>

</div>

</div>
