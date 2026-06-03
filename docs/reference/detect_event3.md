<div id="main" class="col-md-9" role="main">

# Detect marine heatwave events in gridded data

<div class="ref-description section level2">

Detects marine heatwave (or cold-spell) events for each pixel in a
gridded dataset, using climatologies computed by `ts2clm3`. This is the
gridded equivalent of `heatwaveR::detect_event()`.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
detect_event3(
  file_in,
  clim_file,
  file_out,
  var_name = NULL,
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
  roundRes = 4L,
  return_df = FALSE,
  save_file = NULL,
  n_threads = 1L,
  skip_bad_files = FALSE
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   file\_in:

    Path to the input NetCDF file containing SST data. May also be a
    character vector of file paths or a directory path containing daily
    `.nc`/`.nc4` files.

-   clim\_file:

    Path to the climatology NetCDF file produced by `ts2clm3`.

-   file\_out:

    Path for the output NetCDF file containing detected events.

-   var\_name:

    Name of the SST variable. If `NULL`, auto-detected.

-   minDuration:

    Minimum duration (days) for an event. Default `5`.

-   minDuration2:

    Minimum duration for events that also satisfy `threshClim2`. Used
    only when `threshClim2` is supplied.

-   joinAcrossGaps:

    Logical. Join events separated by short gaps? Default `TRUE`.

-   maxGap:

    Maximum gap length (days) to join across. Default `2`.

-   maxGap2:

    Maximum gap length for the secondary `threshClim2` criterion. Used
    only when `threshClim2` is supplied.

-   threshClim2:

    Optional gridded NetCDF logical criterion for the secondary event
    pass, equivalent to `heatwaveR::detect_event()`'s `threshClim2`. The
    file must align with `file_in`; non-zero and non-missing values are
    treated as `TRUE`. May be a single NetCDF file, a vector of files,
    or a directory of daily `.nc`/`.nc4` files.

-   threshClim2\_var\_name:

    Name of the secondary criterion variable. If `NULL`, it is
    auto-detected.

-   coldSpells:

    Logical. Detect cold-spells instead of heatwaves? Default `FALSE`.

-   category:

    Logical. Compute Hobday et al. (2018) severity categories (I
    Moderate through IV Extreme) inline during detection? Categories are
    written to the event NetCDF as additional variables. Default
    `FALSE`.

-   hemisphere:

    Character. Season-naming convention: `"south"` (default, austral:
    DJF = Summer) or `"north"` (boreal: DJF = Winter). Only used when
    `category = TRUE`.

-   roundRes:

    Number of decimal places for rounding event metrics. Default `4`.

-   return\_df:

    Logical. If `TRUE`, return the event table as a `data.frame` in
    addition to writing the NetCDF file. Default `FALSE` (returns the
    file path invisibly).

-   save\_file:

    Optional path for an additional output file. The extension
    determines the format and must be one of `.csv`, `.rds`, or
    `.parquet`. If `NULL`, no companion file is written.

-   n\_threads:

    Number of OpenMP threads. Default `1`.

-   skip\_bad\_files:

    Logical. For multi-file inputs, skip unreadable files or files with
    mismatched grids instead of failing. Default `FALSE`.

</div>

<div class="section level2">

## Value

If `return_df = FALSE` (the default), invisibly returns the path to the
output file. If `return_df = TRUE`, returns a `data.frame` of event
metrics (the NetCDF is still written).

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
# \donttest{
sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
clim_file <- tempfile(fileext = ".nc")
event_file <- tempfile(fileext = ".nc")

ts2clm3(sst_file, clim_file,
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
#> Writing climatology to /var/folders/3w/nmplbnm109b9903rx8z9q0kc0000gn/T//RtmpEGZ5mo/file719a5b0e6202.nc...
#> Done.

detect_event3(sst_file, clim_file, event_file)
#> Reading climatology from /var/folders/3w/nmplbnm109b9903rx8z9q0kc0000gn/T//RtmpEGZ5mo/file719a5b0e6202.nc...
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
#> Writing events to /var/folders/3w/nmplbnm109b9903rx8z9q0kc0000gn/T//RtmpEGZ5mo/file719ad44a002.nc...
#> Done.
# }
```

</div>

</div>

</div>
