# heatwave3 1.3.0 (2026-07-16)

## New features

* **`depth_range`: native depth-resolved (4D) climatology and event
  detection.** `ts2clm3()` and `detect3()` gain
  `depth_range = c(min_depth, max_depth)`, which keeps a contiguous band of
  depth levels as a real dimension instead of squeezing to a single level
  (the existing `depth` index argument is unchanged and still available).
  `detect_event3()` needs no new argument -- it reads the depth range
  straight off the climatology file it already loads, so climatology and
  event detection can never disagree about which levels they are comparing.
  Output NetCDF files write a real `depth` dimension/variable only when more
  than one level is present, so existing 3D output is byte-for-byte
  unchanged. Verified bit-exact against the legacy per-level-squeeze
  approach.
* **Per-day/proto-event output, `hw3_export()`, `category_daily3()`,
  `event_line3()`, and `geom_flame3()` are all depth-aware.** The
  `daily`/`protoEvent` products gain a `depth` dimension when the source
  climatology is depth-resolved; `hw3_export()`'s streaming subset reads
  gain a `depth_range` argument; `category_daily3()` and `event_line3()`
  auto-inherit the depth range from the climatology file the same way
  `detect_event3()` does (`event_line3()` additionally gains an explicit
  `depth` argument, matched to the nearest level actually present in the
  climatology); `geom_flame3()`'s `GeomFlame3` ggproto gains a
  `setup_data()` override that folds a mapped `depth` aesthetic into row
  grouping, so each depth level gets its own flame polygons instead of being
  drawn as one interleaved series.
* **`detect_blob3()` now connects blobs through the water column, not just
  across space and time.** Given a depth-resolved climatology, the
  union-find labeller connects voxels across four axes (longitude,
  latitude, depth, time) instead of three, by depth *index* (adjacent
  level) -- not by a depth range or distance threshold. `event`/`daily`
  gain `depth_min_m`, `depth_max_m`, and `depthOfPeakIntensity_m`, plus the
  volumetric analogues of the existing area metrics: `totalVolume_km3`,
  `peakVolume_km3`, `meanVolume_km3` (volume = footprint area x vertical
  layer thickness, summed over the depth levels a blob occupies that day);
  `voxel` gains a `depth` column, and `rankBy` accepts the new volume
  metrics. The vestigial `connectivity` parameter (which only ever accepted
  `6`) is removed.
* **Two new vignettes.** *Depth-resolved (4D) marine heatwave analysis*
  walks through `depth_range` end to end -- climatology, detection,
  `hw3_export()`, and depth-resolved plotting with
  `geom_flame3()`/`geom_lolli3()`/`event_line3()`. *4D blob connectivity
  through the water column* demonstrates `detect_blob3()`'s depth
  connectivity, including a depth-time "Hovmoeller" view of a single blob's
  vertical structure. Both compile against a small (2x2 pixel, 5
  depth-level) real GLORYS12 fixture bundled with the package
  (`inst/extdata/glorys_depth_test.nc`, plus a 32-file annual split for
  multi-file testing), so they actually run when the documentation site is
  built rather than needing a large external file.
* **New performance benchmark vignette**, *heatwave3 vs xmhw*, comparing
  `heatwave3` against the Python `xmhw` module on an OSTIA Mediterranean
  grid (104 lon x 60 lat, 6,240 ocean pixels, 15,706 daily steps,
  1982-2024). Measured on Ubuntu 24.04 (16 cores, R 4.6.0): `heatwave3`
  completes in 22.5 sec single-threaded (56x faster than the ca. 21 min
  `heatwaveR` serial baseline) and 7.4 sec on 12 threads (170x); `xmhw`
  without `dask` takes ca. 40 min.

## Bug fixes

* **`hw3_export()`'s climatology reader assumed `npixels = nlon*nlat`.** For
  a depth-resolved climatology file this silently misread the data (each
  pixel's depth levels were interleaved into what the reader thought were
  separate horizontal pixels). Fixed alongside the `depth_range` work.
* **Stale `detect3()`/`ts2clm3()` examples were failing `R CMD check
  --as-cran`.** Several `@examples` blocks (`block_average3()`,
  `category3()`, `plot_metric3()`, `event_line3()`, `detect_blob3()`) still
  called `detect3()`/`ts2clm3()` with the pre-1.2.0 positional-argument
  convention (separate `clim_file`/`event_file` paths). Under the current
  `name`-stem signature this let a tempfile path silently land in the
  `var_name` parameter slot, causing a `NetCDF: Variable not found` error
  -- failing the `\donttest{}` examples in GitHub Actions CI without
  failing locally (`devtools::run_examples()` skips `\donttest{}` by
  default).
* **`category3(event_file, clim_file)` silently returned `NA` for every
  category instead of computing them.** `detect_event3()` always writes a
  `category` NetCDF variable to the events file (filled with `0` when
  `category = FALSE` was passed), so `hw3_category()`'s C++ used mere
  presence of that field to decide whether categories were already
  computed -- which was always true -- and so never actually used
  `clim_file`, even when it was supplied precisely because the event file
  had none. Fixed by requiring at least one real (`1`-`4`) category code
  before treating them as pre-computed, matching the check `category3()`
  already used on the R side to raise its "clim_file is required" error.
  Found and fixed while raising test coverage (below).

## Testing

* **Test coverage raised from 77.6% to over 90%** (Codecov), via a
  five-step plan targeting the least-tested code: (1) tests for
  `exceedance3()`, `plot_metric3()`, and `block_average3()`, which
  previously had none at all; (2) tests that force the multi-file/
  directory ingestion path (`hw3_detect_events_multi`,
  `hw3_category_daily_multi`, `hw3_read_sst_multi`) using the
  `glorys_depth_test_annual` fixture, previously untested for
  `detect_event3()`/`category_daily3()`/`detect_blob3()`; (3) direct tests
  for `hw3_version()`, `hw3_pixel_clim()`, `getHW3threads()`/
  `setHW3threads()`, and an explicit `n_threads > 1` run to exercise
  `hw3::parallel_for()`'s worker-pool branch (unreached by default on a
  2-core CI runner, regardless of what `n_threads` a test passed); (4)
  edge-case tests for `export.R` (RDS export, multi-chunk CSV writes,
  empty subset windows, malformed `time_range`, events-product
  `lon_range`/`depth_range` filtering, non-quiet console summaries),
  `event_line3()` (directory input, event-centred window, error paths),
  and `category3()` (file-not-found, deprecated `S` argument). The
  `category3()` bug above was found in the course of this work.

## Documentation

* **Rebuilt pkgdown site** with the two new depth-resolved vignettes and
  the xmhw benchmark vignette indexed and linked from the navbar.
* **README** The Vignettes section now links each title to its pkgdown
  page; the Citation section gained DOIs for the Hobday et al. (2016, 2018)
  references; all external DOI/resource links (README and the pkgdown
  Resources menu) now open in a new tab; the installation section better
  distinguishes the `dev` branch from `main`.

## Internal

* **CI: macOS and Windows `R-CMD-check` jobs failed to find netCDF.**
  Neither runner ships libnetcdf, and `setup-r-dependencies`'s automatic
  system-requirement resolution only covers Linux. Added a
  `brew install netcdf` step for macOS. For Windows, added
  `tools/winlibs.R` plus a rewritten `src/Makevars.win` that downloads a
  prebuilt netcdf-mingw bundle from `rwinlib/netcdf`, mirroring the CRAN
  `ncdf4` package's approach -- ultimately dropped from CI (below) after
  this still failed to link against Rtools' UCRT toolchain (the bundled
  static libraries were built against the legacy MSVCRT runtime).
* **Windows removed from the CI test matrix.** `R-CMD-check` now runs
  macOS and Ubuntu only. Windows remains a best-effort build target via
  `src/Makevars.win` for anyone building locally, but is no longer covered
  by CI.

# heatwave3 1.2.0 (2026-06-04)

## Breaking changes

* **`name`-stem output replaces explicit output paths.** `ts2clm3()`,
  `detect_event3()`, and `detect3()` now take a single `name` argument (a path
  stem) and derive their output filenames from it:
    * `ts2clm3()` writes `paste0(name, "_clim.nc")`.
    * `detect_event3()` writes `paste0(name, "_events.nc")`.
    * the per-day products use `"_events_daily.nc"` and `"_protoevents.nc"`.
  The old `file_out` / `file_out_clim` / `file_out_event` arguments are removed.
  `detect_event3()`'s `clim_file` now defaults to `paste0(name, "_clim.nc")`.
* **`return_df` and `save_file` are removed.** The compute functions always
  write their native NetCDF and (invisibly) return the path(s) written. To read
  a product into R or export it to CSV/RDS/Parquet, use `hw3_export()` (below).

## New features

* **Per-day output, controlled by paths/flags rather than extra return modes.**
    * `detect_event3(daily = "also")` writes a per-day companion NetCDF beside
      the event table; `daily = "only"` writes just the per-day file (no event
      table); `daily = "none"` (default) writes only the event table. The
      per-day product holds `temp`, `seas`, `thresh`, `intensity` (`temp - seas`),
      `event`, `event_no`, and the daily Hobday (2018) `category` (0 = none,
      1 = I Moderate ... 4 = IV Extreme).
    * `protoEvent = TRUE` writes the per-day proto-event series
      (`paste0(name, "_protoevents.nc")`) instead of the event table, mirroring
      `heatwaveR::detect_event(protoEvents = TRUE)`: `temp`, `seas`, `thresh`,
      `threshCriterion`, `durationCriterion`, `event`, `event_no`.
  Both are gridded `[lon, lat, time]` files with a `days since 1970-01-01` time
  axis. Verified against `heatwaveR::detect_event(protoEvents = TRUE)`: identical
  event counts and proto flags pixel-by-pixel.
* **`hw3_export()` is now a read/export hub for every product.** It auto-detects
  the product (climatology, events, daily, proto-events) via the C++ reader and
  either returns a long `data.frame` (whole, or the first `n` rows for a quick
  look) or writes a flat companion file (`.csv`/`.rds`/`.parquet`). It replaces
  the old `return_df`/`save_file` pathways.
* **`category_daily3()`** builds the per-pixel daily marine-heatwave (or
  cold-spell) category series for a date window straight from an SST file, a
  climatology, and an events file, entirely via the C++ readers. It reads only
  the requested SST time-window (a hyperslab) and re-runs no detection: event
  membership is read from the events file, so events that began before the
  window are still recognised (re-detecting on a short window would silently drop
  them). It returns the full daily grid with the same columns as the
  `detect_event3(daily = "also")` product (`lon, lat, t, temp, seas, thresh,
  intensity, event, event_no, category`), with `category` set on event-member
  exceedance days and `NA` elsewhere. This is a NetCDF-native replacement for the
  Marine Heatwave Tracker's tidync-based `load_sub_cat_clim()` (about 240x faster
  per longitude in testing), with an `ice_thresh` argument for the cold-spell
  sea-ice category; the Tracker's event-day subset is `subset(x,
  !is.na(category))`. `sst_file` may be a single multi-time file or a directory /
  vector of daily files (one time step each, e.g. daily-global OSTIA/GHRSST).

## Bug fixes

* **In-process multithreading no longer depends on OpenMP/`libomp`.** The
  parallel loops (climatology and event detection) now use C++ `std::thread`.
  On macOS, a `dlopen`-ed `libomp` could fail to allocate its per-worker
  thread-local storage in a process carrying a large native footprint -- terra /
  raster / sf with their GDAL stack, conda toolchains, or the Positron/`ark`
  kernel that embeds R -- and segfault in `__kmp_suspend` the moment a
  multithreaded run started. `std::thread` uses pthreads and has no such
  dependency, so `n_threads > 1` is now reliable in any session and front-end.
  Results are unchanged (verified pixel-by-pixel). The package no longer detects
  or links an OpenMP runtime; `configure` reports `Parallelism: C++ std::thread`.
* **Noon-stamped daily products are no longer dated a day forward.** CF time
  parsing rounded each timestamp to the nearest day, so daily files stamped at
  12:00:00 (OSTIA, GHRSST, MUR, and similar) were assigned to the *next*
  calendar day. Times are now mapped to the calendar day that contains them
  (floor), matching `heatwaveR`'s `as.Date()` truncation. Midnight-stamped data
  (for example OISST) is unaffected. Climatologies, events, and daily categories
  derived from noon-stamped inputs should be regenerated.
* **Streaming subset extraction in `hw3_export()`.** `vars`, `lon_range`,
  `lat_range`, `time_range`, and `n` select a subset that is read straight from
  the NetCDF as a hyperslab: for the gridded products only the requested
  lon/lat/time window and the chosen variables are pulled from disk, so memory
  and I/O scale with the subset rather than the file. `time_range` takes dates
  for daily/proto-event products (and overlap-filters the events table);
  `lon_range`/`lat_range` filter the events table by per-event coordinate.

## Internal

* **`hw3_read_subset()`** (C++) performs the hyperslab read behind
  `hw3_export()`'s subset path.
* **Post-computation console summary.** `ts2clm3()` and `detect_event3()` print
  the head, tail, and summary statistics of each product written, with a pointer
  to `hw3_export()`. Suppress with `quiet = TRUE`.
* **`hw3_read_daily_nc()`** and **`hw3_file_meta()`** read a per-day NetCDF and a
  product's lightweight metadata (type, dimensions, row count) respectively.

# heatwave3 1.1.4 (2026-06-03)

## Documentation

* **Rebuilt pkgdown site** The public facing documentation website for heatwave3
  was re-compiled. The new vignettes have been added, as well as links to other
  resources and a light/dark theme switch.

## Bug fixes

* **The package now builds on Linux with libstdc++.** The `#include <pthread.h>`
  in `src/hw3_omp.h` sat inside `namespace hw3`, so pthread's symbols
  (`pthread_create`, `pthread_join`, and the rest) were declared as
  `hw3::pthread_*` rather than at global scope. When `hw3_omp.h` was the first
  translation unit to pull in `<pthread.h>` (for instance via `climatology.cpp`
  ahead of Rcpp), libstdc++'s threading shim `<gthr-default.h>` could not
  resolve `::pthread_create` and the build failed on Linux/g++. macOS/libc++
  masked the fault by including `<pthread.h>` globally and early. The guarded
  `#include` now sits at global scope; the fork-handler functions remain inside
  the namespace.

## Internal

* **CI action versions updated.** `actions/upload-artifact@v3` is fully retired
  and hard-failed the `test-coverage` job; `actions/checkout@v3` was on the same
  deprecation path. Both are now pinned to `@v4` across the `R-CMD-check` and
  `test-coverage` workflows. The `r-lib/actions@v2` pins are unchanged, since
  `v2` is the current maintained major.

* **Generated `src/Makevars` is now git-ignored.** The `configure` script writes
  `src/Makevars` with machine-specific absolute paths at build time; only the
  `src/Makevars.in` and `src/Makevars.win` templates are tracked.

# heatwave3 1.1.3 (2026-06-02)

## New features

* **`ts2clm3(return_df = TRUE)`** now returns the computed climatology as a
  long data frame with `lon`, `lat`, `doy`, `seas`, and `thresh`, while still
  writing the primary NetCDF output.

* **Companion exports now use `save_file` only.** `ts2clm3()` and
  `detect_event3()` now infer the requested companion export format from the
  extension of `save_file`. Supported extensions are `.csv`, `.rds`, and
  `.parquet`. If `save_file = NULL`, no companion file is written.

* **`detect_event3()` now implements gridded `threshClim2`.** Users may supply a
  secondary logical criterion as a NetCDF file, vector of files, or directory
  of daily files. Non-zero, non-missing values are treated as `TRUE`, and the
  secondary pass uses `minDuration2` and `maxGap2` with the same event logic as
  `heatwaveR::detect_event()`.

* **CSV and Parquet companion exports are chunked.** `hw3_export()` now streams
  CSV rows and writes Parquet row groups using a configurable `chunk_size`,
  reducing peak memory pressure for large gridded outputs.

## Bug fixes

* **OpenMP now works on macOS.** Multithreading had been silently disabled (the
  package ran single-threaded regardless of `n_threads`) by three faults in the
  `configure` OpenMP probe: the run-test checked the input vector instead of the
  `.C()` return value and so always "failed"; an undefined `SHLIB_CXX17FLAGS`
  spliced an `ERROR:` string into the compile flags; and the probe did not match
  the package's own scheduling. The parallel loops now use `schedule(static, 1)`
  rather than `schedule(dynamic)` (same round-robin load balance, but needing
  only symbols present in every libomp; `schedule(dynamic)` requires
  `__kmpc_dispatch_deinit`, which R's bundled libomp lacks, and that caused
  `dlopen` failures such as `symbol not found ... ___kmpc_dispatch_deinit`). On
  macOS the build now links R's own libomp, so the whole R process shares a
  single OpenMP runtime instead of dragging in a second, conflicting copy. A new
  vignette, *Parallel performance*, documents OpenMP versus R-side (`parallel`)
  scaling and per-platform setup.

* **`detect_event3(return_df = TRUE)` now matches saved event CSV output.**
  The returned data frame now includes the full event-variable set written to
  the event NetCDF and companion CSV, including `pixel_index`, threshold-relative
  intensity metrics, and absolute intensity metrics. Category and season labels
  are returned as character values so that CSV reads and `return_df` output
  compare cleanly.

* **`detect3()` now passes through the current base-function options.** Companion
  output paths can be supplied separately for climatology and event outputs,
  and `detrend` is passed through to `ts2clm3()`.

* **`detect3()` now passes through secondary event-detection options.**
  `threshClim2`, `threshClim2_var_name`, `minDuration2`, and `maxGap2` are
  passed through to `detect_event3()`.

* **NetCDF time calendars and temperature units are handled more defensibly.**
  Non-Gregorian CF calendars now fail with an explicit error instead of being
  silently interpreted as Gregorian, and climatology/event NetCDF outputs now
  preserve the input temperature units instead of hard-coding `degC`.

* **Multi-file input handling is fail-fast by default.** Unreadable files and
  mismatched grids now stop multi-file reads unless `skip_bad_files = TRUE` is
  supplied explicitly.

* **Event categories now use the unrounded peak intensity internally.** Category
  labels and progress values are still written with the requested rounding, but
  category boundaries are no longer affected by rounded `intensity_max` values.

* **`category3()` is now documented and validated as a convenience wrapper.**
  Inline categories from `detect_event3(category = TRUE)` can be read without a
  climatology file; uncategorised event files now require `clim_file` with a
  clear error.

* **`detect_blob3()` now validates fragile arguments.** Unsupported
  connectivity, invalid return components, invalid ranking variables, singleton
  grids, and grid mismatches now fail with explicit errors.

# heatwave3 1.1.2 (2026-06-02)

## Bug fixes

* **`detect_event3()` no longer glues events across calendar gaps in the
  input SST stack.** When the input had missing days (e.g. unavailable
  OSTIA files), the run-length-encoding in `detect_pixel_events()` walked
  array indices rather than calendar dates. An event that had
  above-threshold days on both sides of a missing-data window was joined
  into one event whose reported `date_end − date_start + 1` was larger
  than `duration` (calendar-span > index-count). The fix expands the
  per-pixel series to a dense daily array up front (mirroring what
  `ts2clm3()` has always done in `compute_pixel_clim`), filling missing
  dates with `NA`. `proto_event()` then breaks runs at the gap exactly
  as `heatwaveR::detect_event()` does via `make_whole()`. Verified by
  pixel-by-pixel agreement with heatwaveR (1705 / 1705 events on a
  Cape Town OSTIA region match exactly: same count, same date_start /
  peak / end, same duration, intensities to rounding precision) and by
  a synthetic 30-day-gap reproducer.

  Internal API change: `EventResult` gains three integer fields
  (`jd_start`, `jd_peak`, `jd_end`) holding the absolute Julian Day of
  the start, peak, and end of each event. `index_*` fields are now
  indices into the dense (gap-filled) array rather than the sparse
  input. Downstream code that uses `hw3_read_event_nc()` /
  `category3()` / `block_average3()` is unaffected: the NetCDF
  event-file layout and the R-facing return values are unchanged.

# heatwave3 1.1.1 (2026-05-28)

## Bug fixes

* **`ts2clm3()`** now validates that the parent directory of `file_out`
  exists before running the climatology computation. Previously, an
  invalid path (e.g. a typo such as `dev/test/clim.nc` when only
  `dev/tests/` exists) caused libnetcdf to fail at the *write* step
  with a misleading `Permission denied`, after the full per-pixel
  climatology had already been computed. The new check fires
  immediately with a clear `Output directory does not exist: ...`
  message.

* **`hw3_export()`** now correctly handles climatology files whose grid
  has a singleton longitude or latitude dimension. Previously
  `ncdf4::ncvar_get()` dropped the singleton dimension by default,
  causing `seas[, j, i]` to fail with
  *incorrect number of dimensions*. The variable is now read with
  `collapse_degen = FALSE`.

* **`hw3_export()`** now writes Hobday et al. (2018) **category** and
  **season** values as human-readable labels (`"I Moderate"`,
  `"II Strong"`, …, `"Summer"`, `"Fall"`, …) when exporting an event
  file, matching the labels returned by `category3()`. Previously the
  raw integer codes (1–4) were written, making the exported CSV / RDA /
  Parquet hard to use directly.

# heatwave3 1.1.0 (2026-05-27)

## OpenMP thread management overhaul

Adopts data.table-style OpenMP patterns for safer, more portable
parallelism.

### New features

* **`getHW3threads()` / `setHW3threads()`** — package-level thread
  management. Defaults to 50% of available cores (polite). Overridable
  via the `R_HEATWAVE3_NUM_THREADS` environment variable. The
  per-function `n_threads` parameter still takes precedence when set.

* **Auto-detect OpenMP on macOS** — the `configure` script now performs a
  compile-link-run test (following data.table's approach) that probes
  five OpenMP variants in order: user flags, `-fopenmp`, Apple Clang
  `-Xclang -fopenmp -lomp`, Homebrew libomp (ARM64 and Intel). Each
  test compiles a small program with `schedule(dynamic)` and runs a
  parallel reduction in R, catching runtime symbol conflicts that
  compile-only tests miss. Falls back to single-threaded gracefully.

* **Fork safety** — a `pthread_atfork` handler (registered at package
  load via `.onLoad`) drops to 1 thread before fork and restores
  afterward, preventing deadlocks when heatwave3 is used inside
  `parallel::mclapply()`.

### Internal improvements

* **Never calls `omp_set_num_threads()`** — all parallel regions now use
  `#pragma omp parallel for num_threads(nt)` with a private thread
  count, eliminating the global side effect that changed thread counts
  for other OpenMP-using packages in the R session.

* **`src/hw3_omp.h` compatibility header** — provides stub macros
  (`omp_get_thread_num()`, `omp_get_max_threads()`, etc.) when
  `_OPENMP` is not defined, eliminating `#ifdef _OPENMP` guards from
  source files.

* **Thread throttling** — at least 100 iterations per thread to avoid
  OpenMP overhead on small grids.

---

# heatwave3 1.0.0 (2026-05-27)

## Major: C++ reimplementation

This release is a ground-up rewrite of heatwave3. All core algorithms
(climatology, event detection, categorisation, block averages) are now
implemented in C++17 with OpenMP parallelism and direct libnetcdf I/O.
The previous R/heatwaveR-based implementation is fully replaced.

### New features

* **Inline category computation** — `detect_event3(category = TRUE)` computes
  Hobday et al. (2018) severity categories (I Moderate through IV Extreme)
  during event detection, writing them directly to the event NetCDF. This
  eliminates the separate `category3()` call and the second read of the
  climatology file.

* **Cold-spell categories** — category computation now uses `|seas - thresh|`
  as the category width, correctly handling both heatwaves (thresh > seas)
  and cold-spells (thresh < seas). Previously, `category3()` returned all
  NAs for cold-spell events.

* **`hemisphere` parameter** — replaces the cryptic `S = TRUE` boolean with
  `hemisphere = "south"` (or `"north"`). Used for season naming in category
  output. Backward-compatible: the `S` parameter is still accepted in
  `category3()`.

* **`return_df` parameter** — `detect_event3(return_df = TRUE)` returns the
  event table as a tidy `data.frame` directly, avoiding the write-then-read
  round-trip for interactive use.

* **Progress reporting** — climatology and event detection now emit progress
  updates (percentage complete) from the C++ OpenMP loops, visible in the R
  console during long-running jobs.

* **`hw3_read_metric_summary()`** — new C++ function that reads a single
  metric from the event NetCDF and aggregates per pixel in C++, for efficient
  spatial plotting without reading the full event table into R.

* **`detect3()` as primary entry point** — now accepts `category`,
  `hemisphere`, and `return_df` parameters, making it the single function
  needed for a complete analysis pipeline.

### Eliminated ncdf4 dependency

* All four R functions that previously used the ncdf4 package (`category3`,
  `block_average3`, `plot_metric3`, `event_line3`) now use the package's
  own C++ NetCDF reader. This eliminates a class of segfaults caused by
  symbol collisions when ncdf4 and heatwave3 both linked against libnetcdf
  in the same R session.

### Performance

* `category3()` and `block_average3()` are now pure C++. On a 5.6M-event
  Benguela dataset, `category3()` completes in ~3 seconds (previously
  ~15 minutes in R with ncdf4).
* `plot_metric3()` aggregation is now C++-backed via
  `hw3_read_metric_summary()`.

### Other improvements

* Event NetCDF files now contain `category`, `p_moderate`–`p_extreme`,
  `season`, and `hemisphere` variables when produced with `category = TRUE`.
* `read_event_netcdf()` gracefully handles older event files that lack
  category fields.
* Cleaned up ggplot2 code in `event_line3()` and `plot_metric3()`: proper
  `expression()` axis labels, composable plot objects, no hard-coded colours.

### Breaking changes

* `category3(S = TRUE)` is deprecated in favour of
  `category3(hemisphere = "south")`. The `S` parameter is still accepted
  but will be removed in a future release.
* Event NetCDF files produced with `category = TRUE` contain additional
  variables. This does not break reading by older code, which simply
  ignores unknown variables.

---

# heatwave3 0.0.4 (2024-06-13)

* `detect3()` catches NA pixels before calculating climatologies

# heatwave3 0.0.3 (2024-01-22)

* `detect3()` now converts hourly data to daily before running detection code
* It also checks for pixels with missing data and handles them accordingly

# heatwave3 0.0.2 (2023-12-29)

* Small tweak to get codecov to 100%

# heatwave3 0.0.2 (2023-12-05)

* Pushing codecov up towards 100%
* Changed __`heatwaveR`__ dependence to development version

# heatwave3 0.0.1 (2023-08-31)

* Connecting to codecov

# heatwave3 0.0.1 (2023-08-28)

* Updating examples and adding tests

# heatwave3 0.0.1 (2023-08-27)

* Package now checks and builds correctly
* Added example dataset as CSV for now

# heatwave3 0.0.0.9006 (2023-08-20)

* Added if gates at the start of `detect3()` to catch common errors

# heatwave3 0.0.0.9005 (2023-08-17)

* Added option to output the results as .csv

# heatwave3 0.0.0.9004 (2023-07-30)

* `detect3()` now returns the full event metric output.

# heatwave3 0.0.0.9003 (2023-07-29)

* Initial CRAN submission.
