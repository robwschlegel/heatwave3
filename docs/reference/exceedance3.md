# Static threshold exceedance detection for gridded data

Detects periods where SST exceeds (or falls below) a fixed threshold for
each pixel, without requiring a climatology baseline.

## Usage

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
  n_threads = NULL,
  compression = 1L,
  quiet = FALSE
)
```

## Arguments

- file_in:

  Path to input NetCDF file.

- file_out:

  Path for output event NetCDF file.

- threshold:

  The static threshold value.

- var_name:

  SST variable name. If `NULL`, auto-detected.

- below:

  Logical. If `TRUE`, detect exceedances below the threshold. Default
  `FALSE`.

- minDuration:

  Minimum event duration in days. Default `5`.

- joinAcrossGaps:

  Join events across short gaps? Default `TRUE`.

- maxGap:

  Maximum gap to join. Default `2`.

- lon_range:

  Optional `c(min, max)` longitude range.

- lat_range:

  Optional `c(min, max)` latitude range.

- time_range:

  Optional `c("start", "end")` date range.

- depth:

  Optional depth index.

- roundRes:

  Decimal places for rounding. Default `4`.

- n_threads:

  Number of threads for parallel computation. Default `NULL`: the
  session default, half the available cores unless changed with
  [`setHW3threads`](https://robwschlegel.github.io/heatwave3/reference/getHW3threads.md)
  or the `R_HEATWAVE3_NUM_THREADS` environment variable.

- compression:

  Integer from 0 to 9. zlib compression level of the NetCDF output: 1
  (the default) compresses fastest, 9 gives the smallest files (a few
  percent smaller), and 0 writes uncompressed files, which are several
  times larger but much quicker to write. Gridded outputs are compressed
  on `n_threads` threads when heatwave3 was built with HDF5.

- quiet:

  Logical. Suppress the progress messages from the computation? Warnings
  are still shown. Default `FALSE`.

## Value

Invisibly returns the output file path.

## Examples

``` r
# \donttest{
sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
event_file <- tempfile(fileext = ".nc")

# Detect periods where SST exceeds 295 K
exceedance3(sst_file, event_file, threshold = 295)
#> Reading climatology from /var/folders/w8/8kw9l4sx1v57nx6l9s0rrklw0000gp/T//RtmpYxNu1M/file161515f78b33c.nc...
#> Reading SST data from /Users/rschlegel/Library/R/arm64/4.6/library/heatwave3/extdata/sst_test.nc...
#> Grid: 2 lon x 3 lat x 14276 time = 6 pixels
#> Detecting events with 1 thread(s)...
#> 
  1/6 pixels (16%)
  2/6 pixels (33%)
  3/6 pixels (50%)
  4/6 pixels (66%)
  5/6 pixels (83%)
  6/6 pixels (100%)
#> Found 860 events across 6 pixels
#> Writing events to /var/folders/w8/8kw9l4sx1v57nx6l9s0rrklw0000gp/T//RtmpYxNu1M/file1615177c062a.nc...
#> Done.
# }
```
