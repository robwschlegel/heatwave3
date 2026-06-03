# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What Is This Project?

`heatwave3` is an R package that detects marine heatwaves (MHWs) and cold-spells directly on gridded NetCDF data using the Hobday et al. (2016, 2018) definition. All core algorithms (climatology, event detection, categorisation, block averages) are implemented in C++ with OpenMP parallelism and direct libnetcdf I/O — no R NetCDF packages needed at runtime.

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

- **`netcdf_io.cpp/h`** — Direct NetCDF-C I/O with spatial/temporal subsetting, CF time parsing, scale_factor/add_offset/_FillValue handling. Supports both NC_CHAR and NC_STRING attributes. Reorders data to pixel-major layout [pixel][time] for cache-friendly per-pixel processing. Reads and writes climatology, event, and SST NetCDF files.

- **`climatology.cpp/h`** — Full ts2clm algorithm: gap filling, NA interpolation, DOY assignment (non-leap years skip DOY 60), spread to DOY×year matrix, sliding-window percentile (Type-7 quantile matching R exactly), circular-padded rolling mean smoothing. OpenMP parallelism over pixels with progress reporting.

- **`event_detect.cpp/h`** — Full detect_event algorithm: threshold exceedance, run-length encoding, duration filtering, gap joining, 19 event metric computations. Cold spell support (threshold inversion). Optional inline Hobday et al. (2018) category computation (I Moderate through IV Extreme) and season assignment. OpenMP parallelism over pixels with progress reporting.

- **`blob_label.cpp`** — 3D connected-component labeling (union-find, 6-connectivity, optional dateline wrap). Ported from heatwaveR.

- **`blob_metrics.cpp`** — Per-(blob, day) spatial metric reduction. Ported from heatwaveR.

- **`heatwave3_init.cpp`** — Rcpp exports bridging C++ to R. Key exports: `hw3_compute_clim`, `hw3_detect_events` (with category/hemisphere and protoEvent/write_daily params), `hw3_read_sst`, `hw3_read_event_nc`, `hw3_read_clim_nc`, `hw3_read_daily_nc`, `hw3_category`, `hw3_block_average`, `hw3_read_metric_summary`.

- **`heatwave3_types.h`** — Core structs: `EventResult` (19 metrics + category + season), `EventData` (for reading event NetCDF), `SubsetSpec`, `ClimData`, `GridData`.

### R API (`R/`)

Function names are suffixed with `3` to avoid conflicts with heatwaveR:

**Output naming.** `ts2clm3()`, `detect_event3()`, and `detect3()` take a single `name` (path stem) and derive output filenames: `<name>_clim.nc`, `<name>_events.nc`, `<name>_events_daily.nc`, `<name>_protoevents.nc`. They always write their native NetCDF and return the written path(s) invisibly; there is no `return_df`/`save_file`. To read a product into R or export it to CSV/RDS/Parquet, use `hw3_export()`.

- **`detect3()`** — **Primary entry point.** All-in-one: climatology + detection + optional inline categories. Returns the written paths invisibly.
- **`ts2clm3()`** — NetCDF → climatology NetCDF (seas + thresh per DOY per pixel) at `<name>_clim.nc`.
- **`detect_event3()`** — SST + climatology → event NetCDF. `clim_file` defaults to `<name>_clim.nc`. Supports `category = TRUE` (inline categories) and `hemisphere = "south"/"north"`. Per-day output via `daily = c("none","also","only")` (writes `<name>_events_daily.nc` with temp/seas/thresh/intensity/event/event_no + daily category) and `protoEvent = TRUE` (writes `<name>_protoevents.nc` with the heatwaveR proto flags, instead of the event table). All per-day products are gridded `[lon, lat, time]`.
- **`hw3_export()`** — Read/export hub: auto-detects the product (via `hw3_file_meta`), returns a long data.frame (whole or first `n` rows) or writes a flat `.csv`/`.rds`/`.parquet`. Reading is C++-backed.
- **`category3()`** — Reads pre-computed categories from event file, or computes from climatology file for older events. `hemisphere` parameter replaces deprecated `S`.
- **`block_average3()`** — Yearly event metric aggregation (pure C++)
- **`exceedance3()`** — Static threshold exceedance
- **`detect_blob3()`** — 3D spatially-connected blob detection
- **`event_line3()`** — Per-pixel time series plot from NetCDF (C++ reader, no ncdf4)
- **`geom_flame3()`** / **`geom_lolli3()`** — ggplot2 custom geoms
- **`plot_metric3()`** — Spatial map of event metrics (C++-backed aggregation via `hw3_read_metric_summary`)

### Data Flow

```
NetCDF (SST) → ts2clm3(name="X") → X_clim.nc
                                       ↓
NetCDF (SST) + X_clim.nc → detect_event3(name="X", category=TRUE) → X_events.nc
                                                                        ↓
                            category3() / block_average3() / plotting / hw3_export()
```

Or all-in-one:
```
NetCDF (SST) → detect3(name="X", category=TRUE) → X_clim.nc + X_events.nc
hw3_export("X_events.nc") → data.frame  (or .csv/.rds/.parquet)
```

### Key Design Decisions

1. **Pixel-major memory layout**: SST data reordered to [pixel][time] for contiguous per-pixel time series access during OpenMP parallel loops.
2. **Type-7 quantile**: Matches R's default `quantile()` exactly via linear interpolation.
3. **CF ragged array**: Event NetCDF uses an `event` dimension with lon/lat coordinate variables — avoids wasteful padding to max_events.
4. **No R NetCDF dependency at runtime**: Uses libnetcdf C API directly via `configure` script. All four functions that previously used ncdf4 (`category3`, `block_average3`, `plot_metric3`, `event_line3`) now use the C++ reader.
5. **Inline categories**: Category computation uses the same climatology data already loaded during detection — zero extra I/O.

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

Benguela analysis (93,600 pixels, 16,049 daily files, 12 threads):
- Climatology: ~1 min
- Detection + categories: ~3 min
- Total: ~4 min end-to-end
