# Detect marine heatwave events in gridded data

Detects marine heatwave (or cold-spell) events for each pixel in a
gridded dataset, using climatologies computed by
[`ts2clm3`](https://robwschlegel.github.io/heatwave3/index.html/reference/ts2clm3.md).
This is the gridded equivalent of
[`heatwaveR::detect_event()`](https://rdrr.io/pkg/heatwaveR/man/detect_event.html).

## Usage

``` r
detect_event3(
  file_in,
  name,
  clim_file = NULL,
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
  daily = c("none", "also", "only"),
  protoEvent = FALSE,
  n_threads = 1L,
  skip_bad_files = FALSE,
  quiet = FALSE
)
```

## Arguments

- file_in:

  Path to the input NetCDF file containing SST data. May also be a
  character vector of file paths or a directory path containing daily
  `.nc`/`.nc4` files.

- name:

  Output base name (a path stem). Events are written to
  `paste0(name, "_events.nc")`; the per-day products use
  `"_events_daily.nc"` and `"_protoevents.nc"`. The directory must
  exist.

- clim_file:

  Path to the climatology NetCDF from
  [`ts2clm3`](https://robwschlegel.github.io/heatwave3/index.html/reference/ts2clm3.md).
  When `NULL` (the default), `paste0(name, "_clim.nc")` is used.

- var_name:

  Name of the SST variable. If `NULL`, auto-detected.

- minDuration:

  Minimum duration (days) for an event. Default `5`.

- minDuration2:

  Minimum duration for events that also satisfy `threshClim2`. Used only
  when `threshClim2` is supplied.

- joinAcrossGaps:

  Logical. Join events separated by short gaps? Default `TRUE`.

- maxGap:

  Maximum gap length (days) to join across. Default `2`.

- maxGap2:

  Maximum gap length for the secondary `threshClim2` criterion. Used
  only when `threshClim2` is supplied.

- threshClim2:

  Optional gridded NetCDF logical criterion for the secondary event
  pass, equivalent to
  [`heatwaveR::detect_event()`](https://rdrr.io/pkg/heatwaveR/man/detect_event.html)'s
  `threshClim2`. The file must align with `file_in`; non-zero and
  non-missing values are treated as `TRUE`. May be a single NetCDF file,
  a vector of files, or a directory of daily `.nc`/`.nc4` files.

- threshClim2_var_name:

  Name of the secondary criterion variable. If `NULL`, it is
  auto-detected.

- coldSpells:

  Logical. Detect cold-spells instead of heatwaves? Default `FALSE`.

- category:

  Logical. Compute Hobday et al. (2018) severity categories (I Moderate
  through IV Extreme) inline during detection? Categories are written to
  the event NetCDF as additional variables. Default `FALSE`.

- hemisphere:

  Character. Season-naming convention: `"south"` (default, austral: DJF
  = Summer) or `"north"` (boreal: DJF = Winter). Only used when
  `category = TRUE`.

- roundRes:

  Number of decimal places for rounding event metrics. Default `4`.

- daily:

  Per-day output control, one of `"none"` (default; only the event
  table), `"also"` (event table plus a per-day NetCDF), or `"only"`
  (per-day NetCDF, no event table). The per-day product is written to
  `paste0(name, "_events_daily.nc")` and contains `temp`, `seas`,
  `thresh`, `intensity` (`temp - seas`), `event`, `event_no`, and the
  daily Hobday (2018) `category` (0 = none, 1 = I Moderate ... 4 = IV
  Extreme).

- protoEvent:

  Logical. Write the per-day proto-event series to
  `paste0(name, "_protoevents.nc")` instead of the event table,
  mirroring `heatwaveR::detect_event(protoEvents = TRUE)`. Variables:
  `temp`, `seas`, `thresh`, `threshCriterion`, `durationCriterion`,
  `event`, `event_no`. Mutually exclusive with `daily`. Default `FALSE`.

- n_threads:

  Number of threads for parallel computation. Default `1`.

- skip_bad_files:

  Logical. For multi-file inputs, skip unreadable files or files with
  mismatched grids instead of failing. Default `FALSE`.

- quiet:

  Logical. Suppress the post-computation console summary (head, tail,
  and summary statistics of each product written)? Default `FALSE`.

## Value

Invisibly returns a named character vector of the NetCDF files written
(a subset of `events`, `daily`, `protoevents`). Use
[`hw3_export`](https://robwschlegel.github.io/heatwave3/index.html/reference/hw3_export.md)
to read any of them into a `data.frame` or to export to CSV/RDS/Parquet.

## Examples

``` r
# \donttest{
sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
stem <- file.path(tempdir(), "demo")

ts2clm3(sst_file, name = stem,
        climatologyPeriod = c("1982-01-01", "2011-12-31"))
#> Reading SST data from /tmp/RtmpeUMl7g/temp_libpatheaf5a29c3cf/heatwave3/extdata/sst_test.nc...
#> Grid: 2 lon x 3 lat x 14276 time = 6 pixels
#> Computing climatology with 1 thread(s)...
#>   1/6 pixels (16%)  2/6 pixels (33%)  3/6 pixels (50%)  4/6 pixels (66%)  5/6 pixels (83%)  6/6 pixels (100%)
#> Writing climatology to /tmp/RtmpKE019E/demo_clim.nc...
#> Done.
#> 
#> ------------------------------------------------------------------
#> Climatology written to: /tmp/RtmpKE019E/demo_clim.nc
#> Rows (long format): 2,196   grid: 2 lon x 3 lat
#> 
#> Head:
#>      lon     lat doy     seas   thresh
#> 1 26.525 -34.125   1 294.4208 295.9951
#> 2 26.525 -34.125   2 294.4648 296.0311
#> 3 26.525 -34.125   3 294.5088 296.0720
#> 4 26.525 -34.125   4 294.5524 296.1133
#> 5 26.525 -34.125   5 294.5955 296.1553
#> 
#> Tail:
#>         lon     lat doy     seas   thresh
#> 2192 26.575 -34.025 362 293.5308 295.1848
#> 2193 26.575 -34.025 363 293.5672 295.2141
#> 2194 26.575 -34.025 364 293.6065 295.2460
#> 2195 26.575 -34.025 365 293.6480 295.2776
#> 2196 26.575 -34.025 366 293.6907 295.3100
#> 
#> Summary:
#>   ocean pixels (valid climatology): 6
#>   seas:   291.1 to 295.6
#>   thresh: 292.4 to 297.6
#> 
#> Examine with  hw3_export("/tmp/RtmpKE019E/demo_clim.nc", n = 20)
#> or export with hw3_export("/tmp/RtmpKE019E/demo_clim.nc", file_out = "out.csv")  (.csv/.rds/.parquet)
#> ------------------------------------------------------------------

detect_event3(sst_file, name = stem)                 # demo_events.nc
#> Reading climatology from /tmp/RtmpKE019E/demo_clim.nc...
#> Reading SST data from /tmp/RtmpeUMl7g/temp_libpatheaf5a29c3cf/heatwave3/extdata/sst_test.nc...
#> Grid: 2 lon x 3 lat x 14276 time = 6 pixels
#> Detecting events with 1 thread(s)...
#>   1/6 pixels (16%)  2/6 pixels (33%)  3/6 pixels (50%)  4/6 pixels (66%)  5/6 pixels (83%)  6/6 pixels (100%)
#> Found 610 events across 6 pixels
#> Writing events to /tmp/RtmpKE019E/demo_events.nc...
#> Done.
#> 
#> ------------------------------------------------------------------
#> Events written to: /tmp/RtmpKE019E/demo_events.nc
#> Rows (long format): 610
#> 
#> Head:
#>      lon     lat pixel_index event_no date_start  date_peak   date_end duration
#> 1 26.525 -34.125           0        1 1982-11-06 1982-11-14 1982-11-24       19
#> 2 26.525 -34.125           0        2 1983-04-19 1983-04-20 1983-04-23        5
#> 3 26.525 -34.125           0        3 1983-05-27 1983-05-30 1983-06-01        6
#> 4 26.525 -34.125           0        4 1983-06-24 1983-06-25 1983-06-30        7
#> 5 26.525 -34.125           0        5 1983-07-10 1983-07-12 1983-07-15        6
#>   intensity_mean intensity_max intensity_var intensity_cumulative
#> 1         2.2099        3.1531        0.5616              41.9874
#> 2         3.4849        3.6150        0.1535              17.4243
#> 3         1.9696        2.0416        0.0668              11.8179
#> 4         1.9655        2.4626        0.3339              13.7587
#> 5         1.5890        1.8050        0.1784               9.5342
#>   intensity_mean_relThresh intensity_max_relThresh intensity_var_relThresh
#> 1                   0.6780                  1.6275                  0.5667
#> 2                   1.6822                  1.8107                  0.1481
#> 3                   0.2178                  0.2885                  0.0702
#> 4                   0.5567                  1.0416                  0.3225
#> 5                   0.2596                  0.4747                  0.1772
#>   intensity_cumulative_relThresh intensity_mean_abs intensity_max_abs
#> 1                        12.8813           295.2789            296.22
#> 2                         8.4110           297.9640            298.10
#> 3                         1.3069           295.6150            295.67
#> 4                         3.8969           294.7100            295.27
#> 5                         1.5575           294.0717            294.29
#>   intensity_var_abs intensity_cumulative_abs rate_onset rate_decline
#> 1            0.5886                  5610.30     0.1840       0.1471
#> 2            0.1573                  1489.82     0.1007       0.0349
#> 3            0.0677                  1773.69     0.0303       0.0282
#> 4            0.3945                  2062.97     0.2077       0.1885
#> 5            0.1824                  1764.43     0.1363       0.1270
#> 
#> Tail:
#>        lon     lat pixel_index event_no date_start  date_peak   date_end
#> 606 26.575 -34.025           5       98 2019-07-07 2019-07-11 2019-07-14
#> 607 26.575 -34.025           5       99 2019-08-30 2019-09-02 2019-09-05
#> 608 26.575 -34.025           5      100 2019-10-21 2019-10-26 2019-11-01
#> 609 26.575 -34.025           5      101 2020-07-05 2020-07-06 2020-07-09
#> 610 26.575 -34.025           5      102 2020-08-24 2020-08-25 2020-08-30
#>     duration intensity_mean intensity_max intensity_var intensity_cumulative
#> 606        8         1.9123        2.6167        0.3509              15.2983
#> 607        7         2.6135        3.2218        0.6238              18.2944
#> 608       12         2.7097        4.2806        0.9613              32.5162
#> 609        5         1.9184        2.4610        0.4567               9.5918
#> 610        7         1.7027        1.8887        0.1673              11.9191
#>     intensity_mean_relThresh intensity_max_relThresh intensity_var_relThresh
#> 606                   0.5331                  1.2395                  0.3545
#> 607                   1.2196                  1.8291                  0.6356
#> 608                   1.2481                  2.8252                  0.9574
#> 609                   0.5193                  1.0566                  0.4505
#> 610                   0.3486                  0.5476                  0.1719
#>     intensity_cumulative_relThresh intensity_mean_abs intensity_max_abs
#> 606                         4.2646           293.7725            294.47
#> 607                         8.5369           293.8743            294.48
#> 608                        14.9767           294.3308            295.87
#> 609                         2.5966           293.8140            294.37
#> 610                         2.4402           292.9757            293.17
#>     intensity_var_abs intensity_cumulative_abs rate_onset rate_decline
#> 606            0.3438                  2350.18     0.2090       0.1836
#> 607            0.6214                  2057.12     0.2155       0.5020
#> 608            0.9940                  3531.97     0.5195       0.4084
#> 609            0.4723                  1469.07     0.3247       0.2915
#> 610            0.1695                  2050.83     0.2322       0.0693
#> 
#> Summary:
#>   events: 610   pixels with events: 6
#>   dates:  1982-11-06 to 2020-09-26
#>   duration (days):     5 to    38
#>   intensity_max:   1.314 to 4.911
#> 
#> Examine with  hw3_export("/tmp/RtmpKE019E/demo_events.nc", n = 20)
#> or export with hw3_export("/tmp/RtmpKE019E/demo_events.nc", file_out = "out.csv")  (.csv/.rds/.parquet)
#> ------------------------------------------------------------------
detect_event3(sst_file, name = stem, daily = "also")  # + demo_events_daily.nc
#> Reading climatology from /tmp/RtmpKE019E/demo_clim.nc...
#> Reading SST data from /tmp/RtmpeUMl7g/temp_libpatheaf5a29c3cf/heatwave3/extdata/sst_test.nc...
#> Grid: 2 lon x 3 lat x 14276 time = 6 pixels
#> Detecting events with 1 thread(s)...
#>   1/6 pixels (16%)  2/6 pixels (33%)  3/6 pixels (50%)  4/6 pixels (66%)  5/6 pixels (83%)  6/6 pixels (100%)
#> Found 610 events across 6 pixels
#> Writing per-day series to /tmp/RtmpKE019E/demo_events_daily.nc...
#> Done.
#> Writing events to /tmp/RtmpKE019E/demo_events.nc...
#> Done.
#> 
#> ------------------------------------------------------------------
#> Events written to: /tmp/RtmpKE019E/demo_events.nc
#> Rows (long format): 610
#> 
#> Head:
#>      lon     lat pixel_index event_no date_start  date_peak   date_end duration
#> 1 26.525 -34.125           0        1 1982-11-06 1982-11-14 1982-11-24       19
#> 2 26.525 -34.125           0        2 1983-04-19 1983-04-20 1983-04-23        5
#> 3 26.525 -34.125           0        3 1983-05-27 1983-05-30 1983-06-01        6
#> 4 26.525 -34.125           0        4 1983-06-24 1983-06-25 1983-06-30        7
#> 5 26.525 -34.125           0        5 1983-07-10 1983-07-12 1983-07-15        6
#>   intensity_mean intensity_max intensity_var intensity_cumulative
#> 1         2.2099        3.1531        0.5616              41.9874
#> 2         3.4849        3.6150        0.1535              17.4243
#> 3         1.9696        2.0416        0.0668              11.8179
#> 4         1.9655        2.4626        0.3339              13.7587
#> 5         1.5890        1.8050        0.1784               9.5342
#>   intensity_mean_relThresh intensity_max_relThresh intensity_var_relThresh
#> 1                   0.6780                  1.6275                  0.5667
#> 2                   1.6822                  1.8107                  0.1481
#> 3                   0.2178                  0.2885                  0.0702
#> 4                   0.5567                  1.0416                  0.3225
#> 5                   0.2596                  0.4747                  0.1772
#>   intensity_cumulative_relThresh intensity_mean_abs intensity_max_abs
#> 1                        12.8813           295.2789            296.22
#> 2                         8.4110           297.9640            298.10
#> 3                         1.3069           295.6150            295.67
#> 4                         3.8969           294.7100            295.27
#> 5                         1.5575           294.0717            294.29
#>   intensity_var_abs intensity_cumulative_abs rate_onset rate_decline
#> 1            0.5886                  5610.30     0.1840       0.1471
#> 2            0.1573                  1489.82     0.1007       0.0349
#> 3            0.0677                  1773.69     0.0303       0.0282
#> 4            0.3945                  2062.97     0.2077       0.1885
#> 5            0.1824                  1764.43     0.1363       0.1270
#> 
#> Tail:
#>        lon     lat pixel_index event_no date_start  date_peak   date_end
#> 606 26.575 -34.025           5       98 2019-07-07 2019-07-11 2019-07-14
#> 607 26.575 -34.025           5       99 2019-08-30 2019-09-02 2019-09-05
#> 608 26.575 -34.025           5      100 2019-10-21 2019-10-26 2019-11-01
#> 609 26.575 -34.025           5      101 2020-07-05 2020-07-06 2020-07-09
#> 610 26.575 -34.025           5      102 2020-08-24 2020-08-25 2020-08-30
#>     duration intensity_mean intensity_max intensity_var intensity_cumulative
#> 606        8         1.9123        2.6167        0.3509              15.2983
#> 607        7         2.6135        3.2218        0.6238              18.2944
#> 608       12         2.7097        4.2806        0.9613              32.5162
#> 609        5         1.9184        2.4610        0.4567               9.5918
#> 610        7         1.7027        1.8887        0.1673              11.9191
#>     intensity_mean_relThresh intensity_max_relThresh intensity_var_relThresh
#> 606                   0.5331                  1.2395                  0.3545
#> 607                   1.2196                  1.8291                  0.6356
#> 608                   1.2481                  2.8252                  0.9574
#> 609                   0.5193                  1.0566                  0.4505
#> 610                   0.3486                  0.5476                  0.1719
#>     intensity_cumulative_relThresh intensity_mean_abs intensity_max_abs
#> 606                         4.2646           293.7725            294.47
#> 607                         8.5369           293.8743            294.48
#> 608                        14.9767           294.3308            295.87
#> 609                         2.5966           293.8140            294.37
#> 610                         2.4402           292.9757            293.17
#>     intensity_var_abs intensity_cumulative_abs rate_onset rate_decline
#> 606            0.3438                  2350.18     0.2090       0.1836
#> 607            0.6214                  2057.12     0.2155       0.5020
#> 608            0.9940                  3531.97     0.5195       0.4084
#> 609            0.4723                  1469.07     0.3247       0.2915
#> 610            0.1695                  2050.83     0.2322       0.0693
#> 
#> Summary:
#>   events: 610   pixels with events: 6
#>   dates:  1982-11-06 to 2020-09-26
#>   duration (days):     5 to    38
#>   intensity_max:   1.314 to 4.911
#> 
#> Examine with  hw3_export("/tmp/RtmpKE019E/demo_events.nc", n = 20)
#> or export with hw3_export("/tmp/RtmpKE019E/demo_events.nc", file_out = "out.csv")  (.csv/.rds/.parquet)
#> ------------------------------------------------------------------
#> 
#> ------------------------------------------------------------------
#> Daily series written to: /tmp/RtmpKE019E/demo_events_daily.nc
#> Rows (long format): 85,656   grid: 2 lon x 3 lat
#> 
#> Head:
#>      lon     lat          t   temp     seas   thresh intensity event event_no
#> 1 26.525 -34.125 1982-01-01 295.58 294.4208 295.9951    1.1592 FALSE       NA
#> 2 26.525 -34.125 1982-01-02 295.62 294.4648 296.0311    1.1552 FALSE       NA
#> 3 26.525 -34.125 1982-01-03 295.58 294.5088 296.0720    1.0712 FALSE       NA
#> 4 26.525 -34.125 1982-01-04 295.78 294.5524 296.1133    1.2276 FALSE       NA
#> 5 26.525 -34.125 1982-01-05 295.85 294.5955 296.1553    1.2545 FALSE       NA
#>   category
#> 1       NA
#> 2       NA
#> 3       NA
#> 4       NA
#> 5       NA
#> 
#> Tail:
#>          lon     lat          t   temp     seas   thresh intensity event
#> 85652 26.575 -34.025 2021-01-27 290.65 294.5295 296.3604   -3.8795 FALSE
#> 85653 26.575 -34.025 2021-01-28 291.12 294.5434 296.3859   -3.4234 FALSE
#> 85654 26.575 -34.025 2021-01-29 291.69 294.5551 296.4046   -2.8651 FALSE
#> 85655 26.575 -34.025 2021-01-30 293.52 294.5648 296.4204   -1.0448 FALSE
#> 85656 26.575 -34.025 2021-01-31 295.21 294.5736 296.4334    0.6364 FALSE
#>       event_no category
#> 85652       NA       NA
#> 85653       NA       NA
#> 85654       NA       NA
#> 85655       NA       NA
#> 85656       NA       NA
#> 
#> Summary:
#>   dates: 1982-01-01 to 2021-01-31
#>   event-days: 5195
#>   category: 1=8922  2=439  3=3
#> 
#> Examine with  hw3_export("/tmp/RtmpKE019E/demo_events_daily.nc", n = 20)
#> or export with hw3_export("/tmp/RtmpKE019E/demo_events_daily.nc", file_out = "out.csv")  (.csv/.rds/.parquet)
#> ------------------------------------------------------------------
# }
```
