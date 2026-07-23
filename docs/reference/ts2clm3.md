# Compute climatology for gridded NetCDF data

Calculates seasonal and threshold climatologies for each pixel in a
gridded NetCDF file, following the Hobday et al. (2016) methodology.
This is the gridded equivalent of
[`heatwaveR::ts2clm()`](https://rdrr.io/pkg/heatwaveR/man/ts2clm.html).

## Usage

``` r
ts2clm3(
  file_in,
  name,
  climatologyPeriod,
  lon_range = NULL,
  lat_range = NULL,
  time_range = NULL,
  depth = NULL,
  depth_range = NULL,
  var_name = NULL,
  maxPadLength = FALSE,
  windowHalfWidth = 5L,
  pctile = 90,
  smoothPercentile = TRUE,
  smoothPercentileWidth = 31L,
  var = FALSE,
  detrend = FALSE,
  roundClm = 4L,
  n_threads = 1L,
  skip_bad_files = FALSE,
  quiet = FALSE
)
```

## Arguments

- file_in:

  Path to a single multi-timestep NetCDF file, **or** a character vector
  of daily NetCDF file paths, **or** a directory path containing daily
  NetCDF files (matched by `.nc` or `.nc4` extension).

- name:

  Output base name (a path stem). The climatology is written to
  `paste0(name, "_clim.nc")`, so `name = "results/cape_coast"` produces
  `results/cape_coast_clim.nc`. The directory must exist.

- climatologyPeriod:

  A character vector of length 2 specifying the start and end dates of
  the baseline period, for example `c("1982-01-01", "2011-12-31")`.

- lon_range:

  Optional numeric vector of length 2: `c(min_lon, max_lon)`. If `NULL`,
  all longitudes are used.

- lat_range:

  Optional numeric vector of length 2: `c(min_lat, max_lat)`. If `NULL`,
  all latitudes are used.

- time_range:

  Optional character vector of length 2: `c("start", "end")`. If `NULL`,
  all time steps are read.

- depth:

  Optional integer depth/level index for 4D data. Default `NULL` (no
  depth subsetting). Selects and squeezes a single level – the
  climatology is computed as if the data were 3D. Mutually exclusive
  with `depth_range`.

- depth_range:

  Optional numeric vector of length 2: `c(min_depth, max_depth)` in the
  depth coordinate's own units (typically metres, positive down). Unlike
  `depth`, this keeps the full contiguous band of depth levels falling
  in this range as a real dimension: the climatology is computed
  independently for every `(lon, lat, depth)` triple, and
  `paste0(name, "_clim.nc")` gains a `depth` dimension/coordinate.
  [`detect_event3`](https://robwschlegel.github.io/heatwave3/index.html/reference/detect_event3.md)
  automatically detects and matches this depth range from the
  climatology file – no separate argument is needed there. Mutually
  exclusive with `depth`. Default `NULL` (no depth subsetting).

- var_name:

  Name of the SST variable in the NetCDF file. If `NULL`, the variable
  is auto-detected.

- maxPadLength:

  Maximum number of consecutive missing days to interpolate. Default
  `FALSE` (no interpolation). Set to an integer to enable.

- windowHalfWidth:

  Half-width of the sliding window for climatology calculation. Default
  `5` (11-day window).

- pctile:

  Percentile for the threshold climatology. Default `90`.

- smoothPercentile:

  Logical. Apply rolling mean smoothing to the climatology? Default
  `TRUE`.

- smoothPercentileWidth:

  Width of the rolling mean window for smoothing. Default `31`.

- var:

  Logical. Compute variance climatology? Default `FALSE`.

- detrend:

  Logical. Remove a linear trend from each pixel's time series before
  computing the climatology? Default `FALSE` (fixed-baseline, Hobday et
  al. 2016). Set to `TRUE` to apply the detrended-baseline approach
  (Jacox et al. 2020).

- roundClm:

  Number of decimal places for rounding. Default `4`. Set to `FALSE` to
  disable.

- n_threads:

  Number of threads for parallel computation. Default `1`.

- skip_bad_files:

  Logical. For multi-file inputs, skip unreadable files or files with
  mismatched grids instead of failing. Default `FALSE`.

- quiet:

  Logical. Suppress the post-computation console summary (head, tail,
  and summary statistics of the climatology)? Default `FALSE`.

## Value

Invisibly returns the path to the climatology NetCDF
(`paste0(name, "_clim.nc")`). Use
[`hw3_export`](https://robwschlegel.github.io/heatwave3/index.html/reference/hw3_export.md)
to read it into a `data.frame` or export it to CSV/RDS/Parquet.

## Examples

``` r
if (FALSE) { # \dontrun{
ts2clm3(file_in = "path/to/sst.nc",
        name = file.path(tempdir(), "cape_coast"),
        climatologyPeriod = c("1982-01-01", "2011-12-31"),
        lon_range = c(25, 26), lat_range = c(-34, -33))
} # }
```
