<div id="main" class="col-md-9" role="main">

# All-in-one marine heatwave detection for gridded data

<div class="ref-description section level2">

The primary entry point for `heatwave3`. Runs `ts2clm3` followed by
`detect_event3` in a single call, with optional inline category
computation.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
detect3(
  file_in,
  file_out_clim,
  file_out_event,
  climatologyPeriod,
  save_file_clim = NULL,
  save_file_event = NULL,
  var_name = NULL,
  lon_range = NULL,
  lat_range = NULL,
  time_range = NULL,
  depth = NULL,
  maxPadLength = FALSE,
  windowHalfWidth = 5L,
  pctile = 90,
  smoothPercentile = TRUE,
  smoothPercentileWidth = 31L,
  detrend = FALSE,
  minDuration = 5L,
  minDuration2 = minDuration,
  joinAcrossGaps = TRUE,
  maxGap = 2L,
  maxGap2 = maxGap,
  threshClim2 = NULL,
  threshClim2_var_name = NULL,
  coldSpells = FALSE,
  category = FALSE,
  hemisphere = "south",
  roundClm = 4L,
  roundRes = 4L,
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

    Path to input NetCDF file (or directory/vector of daily files).

-   file\_out\_clim:

    Path for the output climatology NetCDF file.

-   file\_out\_event:

    Path for the output event NetCDF file.

-   climatologyPeriod:

    Character vector of length 2 with start and end dates of the
    baseline period, for example `c("1991-01-01", "2020-12-31")`.

-   save\_file\_clim:

    Optional companion output path for the climatology. The extension
    must be one of `.csv`, `.rds`, or `.parquet`.

-   save\_file\_event:

    Optional companion output path for the event table. The extension
    must be one of `.csv`, `.rds`, or `.parquet`.

-   var\_name:

    Name of the SST variable. If `NULL`, auto-detected.

-   lon\_range:

    Optional `c(min, max)` longitude range.

-   lat\_range:

    Optional `c(min, max)` latitude range.

-   time\_range:

    Optional `c("start", "end")` date range.

-   depth:

    Optional depth/level index for 4D data.

-   maxPadLength:

    Max consecutive NAs to interpolate. Default `FALSE`.

-   windowHalfWidth:

    Half-width of climatology window. Default `5`.

-   pctile:

    Percentile for threshold. Default `90` (heatwaves); use `10` for
    cold-spells.

-   smoothPercentile:

    Apply rolling mean smoothing? Default `TRUE`.

-   smoothPercentileWidth:

    Smoothing window width. Default `31`.

-   detrend:

    Logical. Remove a linear trend before climatology calculation?
    Default `FALSE`.

-   minDuration:

    Minimum event duration in days. Default `5`.

-   minDuration2:

    Secondary minimum duration. See `detect_event3`.

-   joinAcrossGaps:

    Join events across short gaps? Default `TRUE`.

-   maxGap:

    Maximum gap length to join. Default `2`.

-   maxGap2:

    Secondary maximum gap length. See `detect_event3`.

-   threshClim2:

    Optional gridded NetCDF logical criterion for secondary event
    detection. See `detect_event3`.

-   threshClim2\_var\_name:

    Name of the `threshClim2` variable. If `NULL`, it is auto-detected.

-   coldSpells:

    Detect cold-spells? Default `FALSE`.

-   category:

    Logical. Compute Hobday et al. (2018) severity categories inline?
    Default `FALSE`.

-   hemisphere:

    Character. `"south"` (default) or `"north"`.

-   roundClm:

    Decimal places for climatology rounding. Default `4`.

-   roundRes:

    Decimal places for event metric rounding. Default `4`.

-   return\_df:

    Logical. Return the event table as a `data.frame`? Default `FALSE`.

-   n\_threads:

    Number of OpenMP threads. Default `1`.

-   skip\_bad\_files:

    Logical. For multi-file inputs, skip unreadable files or files with
    mismatched grids instead of failing. Default `FALSE`.

</div>

<div class="section level2">

## Value

If `return_df = FALSE` (the default), invisibly returns a list with
`clim_file` and `event_file` paths. If `return_df = TRUE`, returns a
`data.frame` of event metrics.

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
# \donttest{
sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
clim_file <- tempfile(fileext = ".nc")
event_file <- tempfile(fileext = ".nc")

detect3(sst_file, clim_file, event_file,
        climatologyPeriod = c("1982-01-01", "2011-12-31"))
#> Reading SST data from /private/var/folders/3w/nmplbnm109b9903rx8z9q0kc0000gn/T/RtmpEGZ5mo/temp_libpath719a3f2010e4/heatwave3/extdata/sst_test.nc...
#> Grid: 2 lon x 3 lat x 14276 time = 6 pixels
#> Computing climatology with 1 thread(s)...
#> 
  1/6 pixels (16%)
  2/6 pixels (33%)
  3/6 pixels (50%)
  4/6 pixels (66%)
  5/6 pixels (83%)
  6/6 pixels (100%)
#> Writing climatology to /var/folders/3w/nmplbnm109b9903rx8z9q0kc0000gn/T//RtmpEGZ5mo/file719a6954d44e.nc...
#> Done.
#> Reading climatology from /var/folders/3w/nmplbnm109b9903rx8z9q0kc0000gn/T//RtmpEGZ5mo/file719a6954d44e.nc...
#> Reading SST data from /private/var/folders/3w/nmplbnm109b9903rx8z9q0kc0000gn/T/RtmpEGZ5mo/temp_libpath719a3f2010e4/heatwave3/extdata/sst_test.nc...
#> Grid: 2 lon x 3 lat x 14276 time = 6 pixels
#> Detecting events with 1 thread(s)...
#> 
  1/6 pixels (16%)
  2/6 pixels (33%)
  3/6 pixels (50%)
  4/6 pixels (66%)
  5/6 pixels (83%)
  6/6 pixels (100%)
#> Found 610 events across 6 pixels
#> Writing events to /var/folders/3w/nmplbnm109b9903rx8z9q0kc0000gn/T//RtmpEGZ5mo/file719a403e82e8.nc...
#> Done.
# }
```

</div>

</div>

</div>
