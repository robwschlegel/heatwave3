# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What Is This Project?

`heatwave3` is an R package that detects marine heatwaves (MHWs) and cold-spells (MCS) directly on gridded NetCDF data using the Hobday et al. (2016, 2018) definition. All core algorithms (climatology, event detection, categorisation, block averages) are implemented in C++ with `std::thread` parallelism and direct libnetcdf I/O — no R NetCDF packages needed at runtime.

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

Documentation is generated with the newest roxygen2 (`Config/roxygen2/version` in `DESCRIPTION`). When a newer roxygen2 changes `DESCRIPTION`, `NAMESPACE` or `man/`, commit those changes (in their own commit) rather than reverting them; do not hand-edit generated files to match an older version.

Python package (in development on the `python-port` branch; same C++ core):

```bash
pip install .            # builds python/src/bindings.cpp + the shared src/ core
python -m pytest         # runs python/tests
python -m build --sdist  # sdist is an explicit allowlist (pyproject.toml)
```

## Architecture

### C++ Core (`src/`)

All compute-intensive work is in C++. The core is language-neutral (no `Rcpp.h`) so the planned Python package can share it: only `heatwave3_init.cpp`, `blob_label.cpp`, `blob_metrics.cpp` and `RcppExports.cpp` use Rcpp, and they only convert types. Core messages go through `hw3::out()`/`hw3::err()` (`hw3_log.h`, defined per binding); core errors are thrown as `std::runtime_error`.

- **`netcdf_io.cpp/h`** — Direct NetCDF-C I/O with spatial/temporal subsetting, CF time parsing, scale_factor/add_offset/_FillValue handling. Supports both NC_CHAR and NC_STRING attributes. Reorders data to pixel-major layout [pixel][time] for cache-friendly per-pixel processing. Reads and writes climatology, event, and SST NetCDF files.

- **`climatology.cpp/h`** — Full ts2clm algorithm: gap filling, NA interpolation, DOY assignment (non-leap years skip DOY 60), spread to DOY×year matrix, sliding-window percentile (Type-7 quantile matching R exactly), circular-padded rolling mean smoothing. `std::thread` parallelism over pixels with progress reporting.

- **`event_detect.cpp/h`** — Full detect_event algorithm: threshold exceedance, run-length encoding, duration filtering, gap joining, 19 event metric computations. Cold spell support (threshold inversion). Optional inline Hobday et al. (2018) category computation (I Moderate through IV Extreme) and season assignment. `std::thread` parallelism over pixels with progress reporting. Per-day output goes to `DailyBuffers` (float for seas/thresh/intensity, any field may be null); the buffers arrive uninitialised and each worker sets its pixel's slice to the defaults first.

- **`hw3_api.cpp/h`**: language-neutral entry points (plain structs in and out) holding all orchestration: read/compute/write climatology and events, `category_daily` (daily categories for a window), `category` fallback, `block_average`, `metric_summary`, date parsing and `SubsetSpec` construction. Missing values are NaN (doubles) and `hw3::NA_INT` (= R's `NA_INTEGER`).

- **`hw3_round.h`**: ports of R's `round()`, `mean()`, `sum()` and `sd()` (long double accumulation, mean refinement pass), used wherever heatwaveR rounds or summarises, so results are bit-identical to heatwaveR.

- **`hw3_h5chunk.cpp/h`**: parallel compression of gridded outputs (clim, daily, protoevents). The writers in `netcdf_io.cpp` define the file with netCDF and close it; `h5_write_chunks()` reopens it with HDF5, converts, shuffles and deflates each chunk on worker threads exactly as HDF5's filters do, and stores them with `H5Dwrite_chunk()`, so files are byte-identical to netCDF's serial write. Sources are `ChunkSrc::DOUBLE` (stored as double or float), `FLOAT`, `SCHAR` or `INT`. Compiled only with `HW3_HAVE_HDF5` (HDF5 >= 1.10.3 + zlib found by `configure`/CMake, headers matching the library libnetcdf loads); otherwise, for any variable it cannot reproduce, or with the environment variable `HEATWAVE3_SERIAL_WRITE=1` (used by the tests to compare both paths), netCDF writes on one thread. Event tables (1D columns) go through the same path (`put_gridded_and_close()`).

- **`blob_core.cpp/h`**: 3D/4D connected-component labelling (union-find, optional dateline wrap) and per-(blob, day) metric reduction. Ported from heatwaveR. `blob_label.cpp` / `blob_metrics.cpp` are the thin R wrappers (`label_components_cpp`, `blob_daily_summary_cpp`).

- **`heatwave3_init.cpp`**: Rcpp type conversion only; each export calls `hw3_api`. Converts core NaN back to `NA_real_` where R expects NA. Key exports: `hw3_compute_clim`, `hw3_detect_events` (with category/hemisphere and protoEvent/write_daily params), `hw3_read_sst`, `hw3_read_event_nc`, `hw3_read_clim_nc`, `hw3_read_daily_nc`, `hw3_category`, `hw3_block_average`, `hw3_read_metric_summary`.

- **`heatwave3_types.h`** — Core structs: `EventResult` (19 metrics + category + season), `EventData` (for reading event NetCDF), `SubsetSpec`, `ClimData`, `GridData`.

### R API (`R/`)

Function names are suffixed with `3` to avoid conflicts with heatwaveR:

**Output naming.** `ts2clm3()`, `detect_event3()`, and `detect3()` take a single `name` (path stem) and derive output filenames: `<name>_clim.nc`, `<name>_events.nc`, `<name>_events_daily.nc`, `<name>_protoevents.nc`. They always write their native NetCDF and return the written path(s) invisibly; there is no `return_df`/`save_file`. To read a product into R or export it to CSV/RDS/Parquet, use `hw3_export()`.

- **`detect3()`** — **Primary entry point.** All-in-one: climatology + detection + optional inline categories. Returns the written paths invisibly.
- **`ts2clm3()`** — NetCDF → climatology NetCDF (seas + thresh per DOY per pixel) at `<name>_clim.nc`.
- **`detect_event3()`** — SST + climatology → event NetCDF. `clim_file` defaults to `<name>_clim.nc`. Supports `category = TRUE` (inline categories) and `hemisphere = "south"/"north"`. Per-day output via `daily = c("none","also","only")` (writes `<name>_events_daily.nc` with temp/seas/thresh/intensity/event/event_no + daily category on event days only; `slim = TRUE` writes intensity + category only and does not allocate the other buffers) and `protoEvent = TRUE` (writes `<name>_protoevents.nc` with the heatwaveR proto flags, instead of the event table). All per-day products are gridded `[lon, lat, time]`.
- **`hw3_export()`** — Read/export hub: auto-detects the product (via `hw3_file_meta`), returns a long data.frame (whole or first `n` rows) or writes a flat `.csv`/`.rds`/`.parquet`. Reading is C++-backed.
- **`category3()`** — Reads pre-computed categories from event file, or computes them (as `heatwaveR::category()`) from `clim_file` + `sst_file` for event files without them. `hemisphere` parameter replaces deprecated `S`.
- **`category_daily3()`** — Per-pixel daily MHW/MCS series for a date window, built in C++ from SST + clim + events (windowed SST hyperslab, no re-detection; event membership read from the events file). Returns the full daily grid with the same columns as the `daily="also"` product (`lon, lat, t, temp, seas, thresh, intensity, event, event_no, category`); `category` is set on event-member exceedance days and `NA` elsewhere, with `ice_thresh` for the cold-spell ice category (event days with `thresh < ice_thresh`). The MHW Tracker's `load_sub_cat_clim()` event-day subset is `subset(x, !is.na(category))`.
- **`block_average3()`** — Yearly event metric aggregation (pure C++)
- **`exceedance3()`** — Static threshold exceedance
- **`detect_blob3()`** — 4D spatially-connected blob detection
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

1. **Pixel-major memory layout**: SST data reordered to [pixel][time] for contiguous per-pixel time series access during the `std::thread` parallel loops.
2. **Type-7 quantile**: Matches R's default `quantile()` exactly via linear interpolation.
3. **CF ragged array**: Event NetCDF uses an `event` dimension with lon/lat coordinate variables — avoids wasteful padding to max_events.
4. **No R NetCDF dependency at runtime**: Uses libnetcdf C API directly via `configure` script. All four functions that previously used ncdf4 (`category3`, `block_average3`, `plot_metric3`, `event_line3`) now use the C++ reader.
5. **Inline categories**: Category computation uses the same climatology data already loaded during detection — zero extra I/O.
6. **Event-table row order is not deterministic with `n_threads > 1` (investigated, deliberately not fixed, 2026-10-11)**: `detect_events_grid()` merges per-worker results worker by worker, and pixels are dispatched dynamically, so the row order of `<name>_events.nc` can differ between runs (event content is identical; single-threaded output is always in pixel order). Gridded products (clim, daily, protoevents) are written from pixel-indexed buffers and are always deterministic. Benchmarked merge options (8 workers, ~100 events/pixel, 1M to 9M events): the current merge, 0.27 to 3.66 s; filling a pixel-indexed slot vector then walking it in pixel order (restores the single-threaded order), no measurable difference; a full stable sort by lon/lat/event_no plus column gather, ~70% slower merge (~0.3 s per million events) and a second in-memory copy of the table. Decision: not implemented, because row order has no operational value; downstream code must not rely on event-table order (sort explicitly if needed).

7. **Results are not bit-identical across CPU types (decided 2026-10-11)**: compilers fuse multiply-adds (FMA) by default on ARM (including Apple Silicon) but not on baseline x86-64, so C++ accumulations (notably detect_blob3's per-day metrics) differ in the last bit (<= ~4e-16 relative). Rounded outputs (climatology, event metrics) are unaffected in practice, but exactly tied blobs can be ranked differently and `topN` can select differently at a tie. `-ffp-contract=off` would make every platform bit-identical (x86 results unchanged); it was considered and not adopted. The R reference files are generated on x86-64; on other CPUs `test_parity_r.py` matches blobs by id with rtol 1e-14 and ignores `rank`.

8. **heatwaveR is the source of truth (decided 2026-10-11, v1.4.0)**: every output (climatology, event metrics, `p_*` and event/daily categories, rates, negative cold-spell intensities and rates) must equal heatwaveR's. Two deliberate exceptions, both heatwaveR bugs to be fixed there: seasons (heatwave3 uses the peak month, Dec-Feb = Summer in the south; `heatwaveR::category()` intends this but its `$mo`/`$day` assignments are ignored, giving calendar quarters) and the series end (heatwaveR joins up to `maxGap` trailing non-exceeding days into the last event; heatwave3 ends it on the last exceeding day). Any other difference is a heatwave3 bug.

## Build System

- `configure` — Finds libnetcdf via nc-config, pkg-config, or common paths. Adds rpath on macOS. No OpenMP probe (parallelism is `std::thread`). Optionally finds the libhdf5 that libnetcdf links (`ldd`/`otool -L`) and matching headers, and test-runs a program that checks header and library versions and `H5Dwrite_chunk`/`compress2`; on success adds `-DHW3_HAVE_HDF5`, the libhdf5 path and `-lz` (skip with `HW3_NO_HDF5=1`). When comparing builds, reinstall with `R CMD INSTALL --preclean`: stale `src/*.o` keep the previous build's `HW3_HAVE_HDF5` setting.
- `src/Makevars.in` — Template substituted by configure.
- `src/Makevars` — Generated; links `-lnetcdf` (+ libhdf5 and `-lz` when found) + `-pthread` (for `std::thread`; no OpenMP runtime). `src/Makevars.win` has no HDF5, so Windows uses the serial writer.
- SystemRequirements: netcdf (>= 4.0), C++17
- Python: `pyproject.toml` (scikit-build-core + pybind11) and `CMakeLists.txt` at the repo root, excluded from the R build via `.Rbuildignore`. CMake compiles only the language-neutral core files plus `python/src/bindings.cpp` (never the Rcpp files) and finds netCDF via `netCDFConfig.cmake`, then `nc-config`, then pkg-config. HDF5 for parallel compression comes from FindHDF5 (`HDF5_ROOT` points at the wheel prefix in `[tool.cibuildwheel]`), checked with `check_cxx_source_runs`; `-DHW3_HDF5=OFF` disables it. The bindings in `heatwave3._core` mirror the R exports name for name (`hw3_compute_clim`, `hw3_detect_events`, ...). The sdist is an allowlist in `[tool.scikit-build]`: add any new core source file there.
- Python API (`python/heatwave3/`): `climatology.py` (`ts2clm3`), `events.py` (`detect_event3`, `detect3`, `exceedance3`), `category.py` (`category3`, `category_daily3`, `block_average3`), `blob.py` (`detect_blob3`: vectorised mask, shared C++ labelling and daily metrics, pandas event aggregation), `plot.py` (matplotlib `event_line3`, `plot_metric3` returning the `Axes`, with a Python-only `ax` argument; `geom_flame3`/`geom_lolli3` draw onto an `Axes` with Python-native signatures; matplotlib is the optional `[plot]` extra), `export.py` (`hw3_export`, console summary), `threads.py`; `_frames.py` holds the vectorised ports of the R data-frame builders. Names, arguments, defaults and output columns match R; `NULL` -> `None`; R `Date` -> `datetime64`; integer columns that can be `NA` (`event_no`, `category`) are pandas `Int32`. Python-only differences: `hw3_export()` writes `.csv`/`.parquet` (no `.rds`); `detect_blob3(return_=)` because `return` is a Python keyword.
- Python tests: `test_parity_r.py` compares every output with R exactly (except `detect_blob3()` event columns aggregated with sum/mean, rtol 1e-14, because R sums in long double) using golden files in `python/tests/data/golden/`; `test_signatures.py` compares argument names, order and defaults with R. `test_plot.py` checks each plot draws exactly R's data (flame polygons from R's own `GeomFlame3$draw_group()`). Regenerate both references with `Rscript python/tests/make_golden.R` after any change that alters outputs or R signatures, and update the Python side to match.
- Python wheels: `[tool.cibuildwheel]` in `pyproject.toml` builds CPython 3.10-3.14 for Linux x86_64/aarch64 (manylinux_2_28) and macOS arm64/x86_64. `python/tools/build_netcdf.sh PREFIX` builds pinned HDF5 1.14.6 + netcdf-c 4.9.3 (SHA-256 checked; DAP, NCZarr, plugins and SZIP off) with recorded rpaths/install names, so auditwheel/delocate bundle them without LD/DYLD_LIBRARY_PATH. To bump versions, change the version and checksum there. CI: `.github/workflows/python-wheels.yml` (sdist + wheels, tests against the installed wheel); `.github/workflows/python-release.yml` publishes via PyPI Trusted Publishing (manual run: TestPyPI; `py-vX.Y.Z` tag on main: PyPI).
- Python docs (`python/docs/`, Sphinx + pydata-sphinx-theme + numpydoc + myst-nb, built by ReadTheDocs from `.readthedocs.yaml` with warnings as errors): `getting-started.md`, `r-to-python.md`, `api.rst` (autodoc of the numpydoc docstrings), and the notebook via a symlink (rendered from saved outputs, never executed). Figures and printed outputs in `_generated/` come from `python/docs/generate_examples.py`; re-run it after changes that alter outputs. Local build: `python -m sphinx -W -b html python/docs <out>` (Sphinx 9 needs Python >= 3.12).
- `python/notebooks/heatwave3_r_and_python.ipynb`: the same analysis in R (via rpy2) and Python, output figures, R-vs-Python difference maps and per-pixel RMSE, and a `detect_blob3()` comparison.

## Testing

`tests/testthat/test-heatwaveR-parity.R` compares heatwave3 with heatwaveR (when installed) on the bundled `inst/extdata/sst_test.nc`, heatwaves and cold-spells, with exact equality: climatology, every event metric, `category3()` (inline and from SST) and daily categories. Series-end events and seasons are excluded (decision 8). Older tests compare against heatwaveR on OSTIA data from `/Volumes/OceanData/OSTIA_East_Coast_MHW/SWIO_Jan1982-Dec2021.nc` and skip when it is absent.

## Performance

Benchmark on 400 pixels (20×20, 14K time steps):
- heatwave3: 0.53 sec (4 threads) — **~71× faster** than heatwaveR

Benguela analysis (93,600 pixels, 16,049 daily files, 12 threads):
- Climatology: ~1 min
- Detection + categories: ~3 min
- Total: ~4 min end-to-end
