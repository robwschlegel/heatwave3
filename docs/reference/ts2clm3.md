<div id="main" class="col-md-9" role="main">

# Compute climatology for gridded NetCDF data

<div class="ref-description section level2">

Calculates seasonal and threshold climatologies for each pixel in a
gridded NetCDF file, following the Hobday et al. (2016) methodology.
This is the gridded equivalent of `heatwaveR::ts2clm()`.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
ts2clm3(
  file_in,
  file_out,
  climatologyPeriod,
  lon_range = NULL,
  lat_range = NULL,
  time_range = NULL,
  depth = NULL,
  var_name = NULL,
  maxPadLength = FALSE,
  windowHalfWidth = 5L,
  pctile = 90,
  smoothPercentile = TRUE,
  smoothPercentileWidth = 31L,
  var = FALSE,
  detrend = FALSE,
  roundClm = 4L,
  save_file = NULL,
  return_df = FALSE,
  n_threads = 1L,
  skip_bad_files = FALSE
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   file\_in:

    Path to a single multi-timestep NetCDF file, **or** a character
    vector of daily NetCDF file paths, **or** a directory path
    containing daily NetCDF files (matched by `.nc` or `.nc4`
    extension).

-   file\_out:

    Path for the output NetCDF file containing the climatology.

-   climatologyPeriod:

    A character vector of length 2 specifying the start and end dates of
    the baseline period, for example `c("1982-01-01", "2011-12-31")`.

-   lon\_range:

    Optional numeric vector of length 2: `c(min_lon, max_lon)`. If
    `NULL`, all longitudes are used.

-   lat\_range:

    Optional numeric vector of length 2: `c(min_lat, max_lat)`. If
    `NULL`, all latitudes are used.

-   time\_range:

    Optional character vector of length 2: `c("start", "end")`. If
    `NULL`, all time steps are read.

-   depth:

    Optional integer depth/level index for 4D data. Default `NULL` (no
    depth subsetting).

-   var\_name:

    Name of the SST variable in the NetCDF file. If `NULL`, the variable
    is auto-detected.

-   maxPadLength:

    Maximum number of consecutive missing days to interpolate. Default
    `FALSE` (no interpolation). Set to an integer to enable.

-   windowHalfWidth:

    Half-width of the sliding window for climatology calculation.
    Default `5` (11-day window).

-   pctile:

    Percentile for the threshold climatology. Default `90`.

-   smoothPercentile:

    Logical. Apply rolling mean smoothing to the climatology? Default
    `TRUE`.

-   smoothPercentileWidth:

    Width of the rolling mean window for smoothing. Default `31`.

-   var:

    Logical. Compute variance climatology? Default `FALSE`.

-   detrend:

    Logical. Remove a linear trend from each pixel's time series before
    computing the climatology? Default `FALSE` (fixed-baseline, Hobday
    et al. 2016). Set to `TRUE` to apply the detrended-baseline approach
    (Jacox et al. 2020).

-   roundClm:

    Number of decimal places for rounding. Default `4`. Set to `FALSE`
    to disable.

-   save\_file:

    Optional path for an additional output file. The extension
    determines the format and must be one of `.csv`, `.rds`, or
    `.parquet`. If `NULL`, no companion file is written.

-   return\_df:

    Logical. If `TRUE`, return the climatology as a long `data.frame` in
    addition to writing the NetCDF file. Default `FALSE` (returns the
    file path invisibly).

-   n\_threads:

    Number of OpenMP threads for parallel computation. Default `1`.

-   skip\_bad\_files:

    Logical. For multi-file inputs, skip unreadable files or files with
    mismatched grids instead of failing. Default `FALSE`.

</div>

<div class="section level2">

## Value

If `return_df = FALSE` (the default), invisibly returns the path to the
output file. If `return_df = TRUE`, returns a `data.frame` with columns
`lon`, `lat`, `doy`, `seas`, and `thresh`. The NetCDF is still written.

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
if (FALSE) { # \dontrun{
ts2clm3(file_in = "path/to/sst.nc",
        file_out = tempfile(fileext = ".nc"),
        climatologyPeriod = c("1982-01-01", "2011-12-31"),
        lon_range = c(25, 26), lat_range = c(-34, -33))
} # }
```

</div>

</div>

</div>
