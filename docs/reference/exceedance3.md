<div id="main" class="col-md-9" role="main">

# Static threshold exceedance detection for gridded data

<div class="ref-description section level2">

Detects periods where SST exceeds (or falls below) a fixed threshold for
each pixel, without requiring a climatology baseline.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
exceedance3(
  file_in,
  file_out,
  threshold,
  var_name = NULL,
  below = FALSE,
  minDuration = 5L,
  joinAcrossGaps = TRUE,
  maxGap = 2L,
  lon_range = NULL,
  lat_range = NULL,
  time_range = NULL,
  depth = NULL,
  roundRes = 4L,
  n_threads = 1L
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   file\_in:

    Path to input NetCDF file.

-   file\_out:

    Path for output event NetCDF file.

-   threshold:

    The static threshold value.

-   var\_name:

    SST variable name. If `NULL`, auto-detected.

-   below:

    Logical. If `TRUE`, detect exceedances below the threshold. Default
    `FALSE`.

-   minDuration:

    Minimum event duration in days. Default `5`.

-   joinAcrossGaps:

    Join events across short gaps? Default `TRUE`.

-   maxGap:

    Maximum gap to join. Default `2`.

-   lon\_range:

    Optional `c(min, max)` longitude range.

-   lat\_range:

    Optional `c(min, max)` latitude range.

-   time\_range:

    Optional `c("start", "end")` date range.

-   depth:

    Optional depth index.

-   roundRes:

    Decimal places for rounding. Default `4`.

-   n\_threads:

    Number of OpenMP threads. Default `1`.

</div>

<div class="section level2">

## Value

Invisibly returns the output file path.

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
if (FALSE) { # \dontrun{
sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
event_file <- tempfile(fileext = ".nc")

# Detect periods where SST exceeds 295 K
exceedance3(sst_file, event_file, threshold = 295)
} # }
```

</div>

</div>

</div>
