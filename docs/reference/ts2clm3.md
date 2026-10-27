# Compute climatology for gridded NetCDF data

Calculates seasonal and threshold climatologies for each pixel in a
gridded NetCDF file, following the Hobday et al. (2016) methodology.
This is the gridded equivalent of
[`heatwaveR::ts2clm()`](https://robwschlegel.github.io/heatwaveR/reference/ts2clm.html).

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
  n_threads = NULL,
  compression = 1L,
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
  [`detect_event3`](https://robwschlegel.github.io/heatwave3/reference/detect_event3.md)
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

- skip_bad_files:

  Logical. For multi-file inputs, skip unreadable files or files with
  mismatched grids instead of failing. Default `FALSE`.

- quiet:

  Logical. Suppress all console output: the progress messages from the
  computation and the post-computation summary (head, tail, and summary
  statistics of the climatology)? Warnings, such as skipped files, are
  still shown. Default `FALSE`.

## Value

Invisibly returns the path to the climatology NetCDF
(`paste0(name, "_clim.nc")`). Use
[`hw3_export`](https://robwschlegel.github.io/heatwave3/reference/hw3_export.md)
to read it into a `data.frame` or export it to CSV/RDS/Parquet.

## Examples

``` r
# \donttest{
sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")

ts2clm3(file_in = sst_file,
        name = file.path(tempdir(), "cape_coast"),
        climatologyPeriod = c("1982-01-01", "2011-12-31"))
#> Reading SST data from /Users/rschlegel/Library/R/arm64/4.6/library/heatwave3/extdata/sst_test.nc...
#> Grid: 2 lon x 3 lat x 14276 time = 6 pixels
#> Computing climatology with 1 thread(s)...
#> 
  1/6 pixels (16%)
  2/6 pixels (33%)
  3/6 pixels (50%)
  4/6 pixels (66%)
  5/6 pixels (83%)
  6/6 pixels (100%)
#> Writing climatology to /var/folders/w8/8kw9l4sx1v57nx6l9s0rrklw0000gp/T//Rtmp4xYh4P/cape_coast_clim.nc...
#> Done.
#> 
#> ------------------------------------------------------------------
#> Climatology written to: /var/folders/w8/8kw9l4sx1v57nx6l9s0rrklw0000gp/T//Rtmp4xYh4P/cape_coast_clim.nc
#> Rows (long format): 2,196   grid: 2 lon x 3 lat
#> 
#> Head:
#>      lon     lat doy     seas   thresh
#> 1 26.525 -34.125   1 294.4175 295.9932
#> 2 26.525 -34.125   2 294.4615 296.0291
#> 3 26.525 -34.125   3 294.5055 296.0700
#> 4 26.525 -34.125   4 294.5491 296.1114
#> 5 26.525 -34.125   5 294.5922 296.1533
#> 
#> Tail:
#>         lon     lat doy     seas   thresh
#> 2192 26.575 -34.025 362 293.5284 295.1877
#> 2193 26.575 -34.025 363 293.5648 295.2171
#> 2194 26.575 -34.025 364 293.6041 295.2490
#> 2195 26.575 -34.025 365 293.6456 295.2806
#> 2196 26.575 -34.025 366 293.6883 295.3129
#> 
#> Summary:
#>   ocean pixels (valid climatology): 6
#>   seas:   291.1 to 295.6
#>   thresh: 292.4 to 297.6
#> 
#> Examine with  hw3_export("/var/folders/w8/8kw9l4sx1v57nx6l9s0rrklw0000gp/T//Rtmp4xYh4P/cape_coast_clim.nc", n = 20)
#> or export with hw3_export("/var/folders/w8/8kw9l4sx1v57nx6l9s0rrklw0000gp/T//Rtmp4xYh4P/cape_coast_clim.nc", file_out = "out.csv")  (.csv/.rds/.parquet)
#> ------------------------------------------------------------------
# }
```
