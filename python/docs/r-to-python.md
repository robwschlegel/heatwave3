# R to Python

The Python package mirrors the R package: the same functions, the same
argument names and defaults, and the same output columns. Both call the same
C++ core, and the test suite checks that every output matches the R package
exactly. Arguments are camelCase (`climatologyPeriod`, `minDuration`,
`coldSpells`) and outputs snake_case (`intensity_max`, `rate_onset`), as in
the original Python `marineHeatWaves` code.

## Functions

| R | Python | Python returns |
|---|---|---|
| `detect3()` | {func}`~heatwave3.detect3` | dict of the paths written |
| `ts2clm3()` | {func}`~heatwave3.ts2clm3` | path of the climatology file |
| `detect_event3()` | {func}`~heatwave3.detect_event3` | dict of the paths written |
| `exceedance3()` | {func}`~heatwave3.exceedance3` | path of the events file |
| `detect_event_array3()` | {func}`~heatwave3.detect_event_array3` | NumPy array, or tuple of arrays |
| `hw3_export()` | {func}`~heatwave3.hw3_export` | DataFrame, or the path written |
| `category3()` | {func}`~heatwave3.category3` | DataFrame |
| `category_daily3()` | {func}`~heatwave3.category_daily3` | DataFrame |
| `block_average3()` | {func}`~heatwave3.block_average3` | DataFrame |
| `detect_blob3()` | {func}`~heatwave3.detect_blob3` | dict of DataFrames |
| `event_line3()` | {func}`~heatwave3.event_line3` | matplotlib `Axes` |
| `plot_metric3()` | {func}`~heatwave3.plot_metric3` | matplotlib `Axes` |
| `geom_flame3()` | {func}`~heatwave3.geom_flame3` | matplotlib `PolyCollection` |
| `geom_lolli3()` | {func}`~heatwave3.geom_lolli3` | matplotlib artists |
| `getHW3threads()`, `setHW3threads()` | {func}`~heatwave3.getHW3threads`, {func}`~heatwave3.setHW3threads` | int |

## Translating arguments

| R | Python |
|---|---|
| `NULL` | `None` |
| `TRUE`, `FALSE` | `True`, `False` |
| `c("1982-01-01", "2011-12-31")` | `("1982-01-01", "2011-12-31")` or a list |
| `c(15, 35)` | `(15, 35)` |
| a vector of file paths | a list of paths (`str` or `pathlib.Path`) |
| dates as character or `Date` | `str`, `datetime.date`, `numpy.datetime64` or `pandas.Timestamp` |
| `5L` | `5` |

For example, in R:

```r
detect3("sst.nc", name = "run1", climatologyPeriod = c("1982-01-01", "2011-12-31"),
        lon_range = c(15, 35), category = TRUE, n_threads = 4)
```

and in Python:

```python
hw3.detect3("sst.nc", name="run1", climatologyPeriod=("1982-01-01", "2011-12-31"),
            lon_range=(15, 35), category=True, n_threads=4)
```

## Outputs

| R | Python |
|---|---|
| `data.frame` | `pandas.DataFrame` |
| `Date` column | `datetime64[ns]` column |
| integer column | `int32`; nullable `Int32` where values can be missing (`event_no`, `category`) |
| `NA` in a numeric column | `NaN` |
| character `NA` | `None` |
| named list | `dict` |
| paths returned invisibly | paths returned normally |

## Differences

- **`detect_blob3(return_=...)`**: the R argument `return` is a Python
  keyword, so it is `return_` in Python.
- **`detect_event_array3()` axis order**: time is the *last* axis in Python
  and the *first* dimension in R. Each language then holds every pixel's
  series contiguously (NumPy arrays are row-major, R arrays column-major), so
  neither copies its input. Several outputs come back as a tuple in Python
  (as `xarray.apply_ufunc` expects) and as a named list in R; `event` is
  `bool`, `event_no` `int32`, `intensity` `float32` and `category` `int8`.
- **`hw3_export()` file formats**: `.csv` and `.parquet` (the latter needs
  `pyarrow`). The R-only `.rds` format is not available.
- **Plotting**: in R, `geom_flame3()` and `geom_lolli3()` are ggplot2 layers.
  In Python they are functions that draw onto a matplotlib `Axes`; with
  `data=`, their arguments can be column names, and `group=` (for example a
  depth column) plays the role of R's automatic grouping by depth.
  `event_line3()` and `plot_metric3()` take the R arguments plus `ax=`, and
  return the `Axes`.
- **Coastlines** in `plot_metric3()` use `cartopy` when it is installed (R uses
  `rnaturalearth`); otherwise they are skipped, as in R.
- **Threads**: `HEATWAVE3_NUM_THREADS` sets the default thread count;
  `R_HEATWAVE3_NUM_THREADS` is also read.

## Numerical agreement

Both packages reproduce heatwaveR: the climatology, the event metrics, the
event and daily categories are identical to `heatwaveR::ts2clm()`,
`detect_event()` and `category()` (including heatwaveR's rounding and its
sign convention for cold-spells, whose intensities are negative). heatwave3
differs from heatwaveR in two places, where heatwaveR has known issues: the
season is that of the event's peak month (December to February = summer in
the southern hemisphere), and an event is never extended over non-exceeding
days at the very end of the series.

For the same input, the R and Python packages give identical results, with two
exceptions in {func}`~heatwave3.detect_blob3`:

- The per-blob totals and means (`cumI_km2_day`, `meanArea_km2`, the
  centroids and the volume summaries) are summed in R with extended (80-bit)
  precision and in Python with 64-bit precision. They agree to about 1e-15
  relative.
- On ARM processors (including Apple Silicon), compilers fuse multiply-adds,
  so the accumulated blob metrics can differ from x86-64 in the last bit, in
  both languages. Blobs with exactly equal ranking values can then be ordered
  and ranked differently, and `topN` can select a different blob at a tie.
