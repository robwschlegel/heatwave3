# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

`heatwave3` is an R package that detects marine heatwaves (MHWs) and cold-spells directly on gridded NetCDF data using the Hobday et al. (2016, 2018) definition. All core algorithms (climatology, event detection) are implemented in C++ with OpenMP parallelism and direct libnetcdf I/O — no R NetCDF packages needed at runtime.

This is a ground-up C++ reimplementation of [heatwaveR](https://robwschlegel.github.io/heatwaveR/), designed for gridded data rather than per-pixel time series.

## Build & Development Commands

```bash
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

## Architecture

### C++ Core (`src/`)

All compute-intensive work is in C++ with Rcpp bindings:

- **`netcdf_io.cpp/h`** — Direct NetCDF-C I/O with spatial/temporal subsetting, CF time parsing, scale_factor/add_offset/_FillValue handling. Supports both NC_CHAR and NC_STRING attributes. Reorders data to pixel-major layout [pixel][time] for cache-friendly per-pixel processing.

- **`climatology.cpp/h`** — Full ts2clm algorithm: gap filling, NA interpolation, DOY assignment (non-leap years skip DOY 60), spread to DOY×year matrix, sliding-window percentile (Type-7 quantile matching R exactly), circular-padded rolling mean smoothing. OpenMP parallelism over pixels.

- **`event_detect.cpp/h`** — Full detect_event algorithm: threshold exceedance, run-length encoding, duration filtering, gap joining, 19 event metric computations. Cold spell support (threshold inversion). OpenMP parallelism over pixels.

- **`blob_label.cpp`** — 3D connected-component labeling (union-find, 6-connectivity, optional dateline wrap). Ported from heatwaveR.

- **`blob_metrics.cpp`** — Per-(blob, day) spatial metric reduction. Ported from heatwaveR.

- **`heatwave3_init.cpp`** — Rcpp exports bridging C++ to R. Contains `hw3_compute_clim`, `hw3_detect_events`, `hw3_read_sst`, and helper exports.

### R API (`R/`)

Function names are suffixed with `3` to avoid conflicts with heatwaveR:

- **`ts2clm3()`** — NetCDF → climatology NetCDF (seas + thresh per DOY per pixel)
- **`detect_event3()`** — SST + climatology → event NetCDF (ragged array of events)
- **`detect3()`** — All-in-one convenience wrapper
- **`exceedance3()`** — Static threshold exceedance
- **`category3()`** — MHW categories (Moderate/Strong/Severe/Extreme)
- **`block_average3()`** — Yearly event metric aggregation
- **`detect_blob3()`** — 3D spatially-connected blob detection
- **`event_line3()`** — Per-pixel time series plot from NetCDF
- **`geom_flame3()`** / **`geom_lolli3()`** — ggplot2 custom geoms
- **`plot_metric3()`** — Spatial map of event metrics

### Data Flow

```
NetCDF (SST) → ts2clm3() → NetCDF (climatology)
                                ↓
NetCDF (SST) + clim → detect_event3() → NetCDF (events)
                                              ↓
                              category3() / block_average3() / plotting
```

### Key Design Decisions

1. **Pixel-major memory layout**: SST data reordered to [pixel][time] for contiguous per-pixel time series access during OpenMP parallel loops.
2. **Type-7 quantile**: Matches R's default `quantile()` exactly via linear interpolation.
3. **CF ragged array**: Event NetCDF uses an `event` dimension with lon/lat coordinate variables — avoids wasteful padding to max_events.
4. **No R NetCDF dependency at runtime**: Uses libnetcdf C API directly via `configure` script.

## Build System

- `configure` — Finds libnetcdf via nc-config, pkg-config, or common paths. Adds rpath on macOS.
- `src/Makevars.in` — Template substituted by configure.
- `src/Makevars` — Generated; links `-lnetcdf` + OpenMP flags.
- SystemRequirements: netcdf (>= 4.0), C++17

## Testing

Tests compare heatwave3 output pixel-by-pixel against heatwaveR on real OSTIA SST data from `/Volumes/OceanData/OSTIA_East_Coast_MHW/SWIO_Jan1982-Dec2021.nc`. Tests verify:
- Climatology seas/thresh match to 1e-4 (rounding precision)
- Event count, duration exact match
- Event intensity metrics match to 2e-4

## Performance

Benchmark on 400 pixels (20×20, 14K time steps):
- heatwave3: 0.53 sec (4 threads) — **~71× faster** than heatwaveR
