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
