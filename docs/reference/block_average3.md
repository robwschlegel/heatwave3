# Block average of event metrics by year and pixel

Computes yearly summary statistics of event metrics for each pixel.

## Usage

``` r
block_average3(event_file)
```

## Arguments

- event_file:

  Path to the event NetCDF file from
  [`detect_event3`](https://robwschlegel.github.io/heatwave3/index.html/reference/detect_event3.md).

## Value

A data.frame with yearly aggregated event metrics per pixel, including
count, mean/max duration, mean/max intensity, and total days.

## Examples

``` r
# \donttest{
sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
clim_file <- tempfile(fileext = ".nc")
event_file <- tempfile(fileext = ".nc")

detect3(sst_file, clim_file, event_file,
        climatologyPeriod = c("1982-01-01", "2011-12-31"))
#> Reading SST data from /private/var/folders/3w/nmplbnm109b9903rx8z9q0kc0000gn/T/RtmppsylXt/temp_libpathec9b4aa2cdfa/heatwave3/extdata/sst_test.nc...
#> Error: NetCDF error in inq varid /var/folders/3w/nmplbnm109b9903rx8z9q0kc0000gn/T//Rtmp3f4glh/file17dc7b8d9864.nc: NetCDF: Variable not found

ba <- block_average3(event_file)
#> Error: NetCDF error in open /var/folders/3w/nmplbnm109b9903rx8z9q0kc0000gn/T//Rtmp3f4glh/file17dc7b8d9864.nc: No such file or directory
head(ba)
#> Error: object 'ba' not found
# }
```
