<div id="main" class="col-md-9" role="main">

# Block average of event metrics by year and pixel

<div class="ref-description section level2">

Computes yearly summary statistics of event metrics for each pixel.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
block_average3(event_file)
```

</div>

</div>

<div class="section level2">

## Arguments

-   event\_file:

    Path to the event NetCDF file from `detect_event3`.

</div>

<div class="section level2">

## Value

A data.frame with yearly aggregated event metrics per pixel, including
count, mean/max duration, mean/max intensity, and total days.

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
#> Writing climatology to /var/folders/3w/nmplbnm109b9903rx8z9q0kc0000gn/T//RtmpEGZ5mo/file719a70833657.nc...
#> Done.
#> Reading climatology from /var/folders/3w/nmplbnm109b9903rx8z9q0kc0000gn/T//RtmpEGZ5mo/file719a70833657.nc...
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
#> Writing events to /var/folders/3w/nmplbnm109b9903rx8z9q0kc0000gn/T//RtmpEGZ5mo/file719a3664c376.nc...
#> Done.

ba <- block_average3(event_file)
head(ba)
#>      lon     lat year count duration_mean duration_max intensity_mean
#> 1 26.525 -34.125 1982     1      19.00000           19       2.209900
#> 2 26.525 -34.125 1983     5       6.20000            7       2.130060
#> 3 26.525 -34.125 1984     2       7.00000            9       2.016500
#> 4 26.525 -34.125 1985     2       5.00000            5       1.632900
#> 5 26.525 -34.125 1986     5       8.40000           14       2.530060
#> 6 26.525 -34.125 1987     3      11.66667           15       2.399233
#>   intensity_max_mean intensity_max_max intensity_cumulative_mean total_days
#> 1           3.153100            3.1531                  41.98740         19
#> 2           2.384700            3.6150                  12.80478         31
#> 3           2.448550            2.7094                  13.67920         14
#> 4           1.723400            1.9355                   8.16455         10
#> 5           3.153800            3.8774                  21.01420         42
#> 6           2.840733            3.1123                  28.25203         35
#>   total_icum
#> 1    41.9874
#> 2    64.0239
#> 3    27.3584
#> 4    16.3291
#> 5   105.0710
#> 6    84.7561
# }
```

</div>

</div>

</div>
