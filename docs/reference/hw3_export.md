<div id="main" class="col-md-9" role="main">

# Export heatwave3 output to additional formats

<div class="ref-description section level2">

Reads a heatwave3 NetCDF output file (climatology or events) and writes
it to a companion file. The output format is inferred from the extension
of `file_out`.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
hw3_export(nc_file, file_out, type = c("clim", "event"), chunk_size = 1000000L)
```

</div>

</div>

<div class="section level2">

## Arguments

-   nc\_file:

    Path to the heatwave3 NetCDF output file.

-   file\_out:

    Path for the exported companion file. The extension must be one of
    `.csv`, `.rds`, or `.parquet`.

-   type:

    Character. One of `"clim"` or `"event"`, indicating the type of
    output file.

-   chunk\_size:

    Number of rows per CSV or Parquet write chunk. Default `1000000`.

</div>

<div class="section level2">

## Value

Invisibly returns the path to the exported file.

</div>

<div class="section level2">

## Details

-   `"csv"`: Writes a flat CSV file in row chunks.

-   `"rds"`: Writes the data frame with `saveRDS()`.

-   `"parquet"`: Writes a single Apache Parquet file in row groups.
    Requires the arrow package.

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
# \donttest{
sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
event_file <- tempfile(fileext = ".nc")
clim_file <- tempfile(fileext = ".nc")

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
#> Writing climatology to /var/folders/3w/nmplbnm109b9903rx8z9q0kc0000gn/T//RtmpEGZ5mo/file719a3a37cb96.nc...
#> Done.
#> Reading climatology from /var/folders/3w/nmplbnm109b9903rx8z9q0kc0000gn/T//RtmpEGZ5mo/file719a3a37cb96.nc...
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
#> Writing events to /var/folders/3w/nmplbnm109b9903rx8z9q0kc0000gn/T//RtmpEGZ5mo/file719a34dd80df.nc...
#> Done.

# Export events as CSV
csv_path <- hw3_export(event_file, file_out = tempfile(fileext = ".csv"),
                       type = "event")
head(read.csv(csv_path))
#>      lon     lat pixel_index event_no date_start  date_peak   date_end duration
#> 1 26.525 -34.125           0        1 1982-11-06 1982-11-14 1982-11-24       19
#> 2 26.525 -34.125           0        2 1983-04-19 1983-04-20 1983-04-23        5
#> 3 26.525 -34.125           0        3 1983-05-27 1983-05-30 1983-06-01        6
#> 4 26.525 -34.125           0        4 1983-06-24 1983-06-25 1983-06-30        7
#> 5 26.525 -34.125           0        5 1983-07-10 1983-07-12 1983-07-15        6
#> 6 26.525 -34.125           0        6 1983-07-20 1983-07-21 1983-07-26        7
#>   intensity_mean intensity_max intensity_var intensity_cumulative
#> 1         2.2099        3.1531        0.5616              41.9874
#> 2         3.4849        3.6150        0.1535              17.4243
#> 3         1.9696        2.0416        0.0668              11.8179
#> 4         1.9655        2.4626        0.3339              13.7587
#> 5         1.5890        1.8050        0.1784               9.5342
#> 6         1.6413        1.9993        0.2555              11.4888
#>   intensity_mean_relThresh intensity_max_relThresh intensity_var_relThresh
#> 1                   0.6780                  1.6275                  0.5667
#> 2                   1.6822                  1.8107                  0.1481
#> 3                   0.2178                  0.2885                  0.0702
#> 4                   0.5567                  1.0416                  0.3225
#> 5                   0.2596                  0.4747                  0.1772
#> 6                   0.3435                  0.6965                  0.2499
#>   intensity_cumulative_relThresh intensity_mean_abs intensity_max_abs
#> 1                        12.8813           295.2789            296.22
#> 2                         8.4110           297.9640            298.10
#> 3                         1.3069           295.6150            295.67
#> 4                         3.8969           294.7100            295.27
#> 5                         1.5575           294.0717            294.29
#> 6                         2.4043           294.0400            294.42
#>   intensity_var_abs intensity_cumulative_abs rate_onset rate_decline
#> 1            0.5886                  5610.30     0.1840       0.1471
#> 2            0.1573                  1489.82     0.1007       0.0349
#> 3            0.0677                  1773.69     0.0303       0.0282
#> 4            0.3945                  2062.97     0.2077       0.1885
#> 5            0.1824                  1764.43     0.1363       0.1270
#> 6            0.2772                  2058.28     0.1262       0.1208

# Export climatology as RDS
hw3_export(clim_file, file_out = tempfile(fileext = ".rds"), type = "clim")
# }
```

</div>

</div>

</div>
