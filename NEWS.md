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
