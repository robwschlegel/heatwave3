# Detect spatially connected marine heatwave blobs

Identifies spatially contiguous marine heatwave (or cold-spell) events
("blobs") by performing connected-component labelling on the exceedance
mask across longitude, latitude, (optionally) depth, and time.

## Usage

``` r
detect_blob3(
  sst_file,
  clim_file,
  var_name = NULL,
  wrapDateline = FALSE,
  minVoxels = 1L,
  topN = NULL,
  rankBy = "cumI_km2_day",
  coldSpells = FALSE,
  return = c("event", "daily"),
  skip_bad_files = FALSE,
  connectivity = 6L
)
```

## Arguments

- sst_file:

  Path to the SST NetCDF file (the same file passed to
  [`ts2clm3`](https://robwschlegel.github.io/heatwave3/reference/ts2clm3.md)).

- clim_file:

  Path to the climatology NetCDF produced by
  [`ts2clm3`](https://robwschlegel.github.io/heatwave3/reference/ts2clm3.md).
  The spatial (and, if depth-resolved, vertical) grid of the climatology
  determines the region analysed.

- var_name:

  Name of the SST variable in `sst_file`. If `NULL`, auto-detected from
  CF attributes.

- wrapDateline:

  Logical. Wrap the longitude axis so that blobs can connect across the
  antimeridian? Default `FALSE`. Depth is never wrapped (not periodic).

- minVoxels:

  Minimum number of voxels (space-time, or space-depth-time for
  depth-resolved input) for a blob to be retained. Default `1`. Set
  higher (for example `200`) to filter out small, short-lived events.

- topN:

  Return only the top N blobs ranked by `rankBy`. Default `NULL` (return
  all).

- rankBy:

  Character. Metric to rank blobs by. Default `"cumI_km2_day"`
  (cumulative intensity-weighted area). Other options include
  `"peakArea_km2"`, `"duration"`, `"meanArea_km2"`, and, for
  depth-resolved input only, `"totalVolume_km3"`, `"peakVolume_km3"`,
  `"meanVolume_km3"`.

- coldSpells:

  Logical. If `TRUE`, detect cold-spell blobs (`temp < thresh`).
  Requires a climatology computed with a low percentile (for example
  `pctile = 10` in
  [`ts2clm3`](https://robwschlegel.github.io/heatwave3/reference/ts2clm3.md)).
  Default `FALSE`.

- return:

  Character vector specifying which components to include in the output.
  One or more of:

  `"event"`

  : One row per blob with summary metrics (duration, peak area,
    cumulative intensity, centroid, etc.)

  `"daily"`

  : One row per (blob, date) with daily area, mean/max anomaly,
    centroid, and bounding box

  `"voxel"`

  : One row per (blob, lon, lat, date) (plus depth, if depth-resolved),
    the full footprint. Required for spatial footprint maps and
    persistence analysis. Can be large.

  Default `c("event", "daily")`.

- skip_bad_files:

  Logical. For multi-file SST inputs, skip unreadable files or files
  with mismatched grids instead of failing. Default `FALSE`.

- connectivity:

  Integer. Which neighbouring voxels in the lon/lat/time cube join a
  blob: `6` (default, faces only), `10` (adds horizontal diagonals),
  `18` (adds edges, including space-time diagonals) or `26` (all
  neighbours). Depth levels join only straight up and down, whatever the
  setting. See **Connectivity**.

## Value

A named list containing the requested components (`event`, `daily`,
and/or `voxel` data.frames).

## How it works

This function does **not** read a pre-computed event file. Instead, it
takes the raw SST data and the climatology, and internally constructs a
boolean exceedance mask by comparing each pixel's daily SST against its
threshold climatology (`temp > thresh` for heatwaves, or `temp < thresh`
for cold-spells). Connected voxels in this mask are then labelled as
coherent spatial (or spatio-vertical) events using a union-find
algorithm.

The reason an event file from
[`detect_event3`](https://robwschlegel.github.io/heatwave3/reference/detect_event3.md)
cannot be used directly is that the event file stores per-pixel event
*summaries* (start/end dates, intensity metrics), not the daily
per-pixel boolean mask needed for spatial connectivity analysis.

## Connectivity

Two exceedance voxels belong to the same blob when they are neighbours,
namely when their grid indices differ by at most one step along each
axis. `connectivity` sets which of these neighbours count, named by
their number in the 3 x 3 x 3 lon/lat/time cube:

- `6` (default):

  Faces only: one axis changes. East, west, north and south on the same
  day, and the same pixel on the previous and next day. Pixels that
  touch only at a corner are not connected (4-connectivity in the
  horizontal plane), and a blob continues from one day to the next only
  through pixels that exceed the threshold on both days.

- `10`:

  Faces plus the four horizontal diagonals on the same day
  (8-connectivity in the horizontal plane). Days still join only through
  pixels that exceed on both days.

- `18`:

  Faces and edges: up to two axes change at once. Adds space-time
  diagonals, so a pixel joins its neighbour on the next day even when
  the two never exceed on the same day.

- `26`:

  All neighbours in the cube, including corners.

The setting names the same neighbourhood whether or not the input is
depth-resolved: with depth, the cube applies on each depth level, and a
voxel also joins the same pixel on the same day one level above and
below (see **Depth connectivity**), so 8, 12, 20 or 28 neighbours in
all. Higher settings only merge blobs, never split them: every blob
found with `6` lies inside one blob found with `26`. Adjacency is by
grid index alone: the distance between neighbours in kilometres, metres
or days plays no role, so a space-time diagonal step stands for a
displacement of one grid cell per day, a speed that depends on the grid.
With `wrapDateline = TRUE` the first and last longitude columns are
neighbours under the same rule.

## Depth connectivity

When `clim_file` is depth-resolved (from `ts2clm3(depth_range = ...)`),
`detect_blob3()` detects this automatically – no separate argument is
needed, the same pattern used by
[`detect_event3`](https://robwschlegel.github.io/heatwave3/reference/detect_event3.md)
and
[`category_daily3`](https://robwschlegel.github.io/heatwave3/reference/category_daily3.md).
Voxels are then also connected through the water column: a voxel joins
the same pixel on the same day at the level immediately above/below it,
by depth *index*, whatever `connectivity` is (there are no depth
diagonals), and whatever that level's actual depth in metres happens to
be, with no minimum/maximum depth gap or distance threshold.
`event`/`daily` gain `depth_min_m`, `depth_max_m`,
`depthOfPeakIntensity_m`, and volumetric analogues of the area metrics
(`totalVolume_km3`/`peakVolume_km3`/`meanVolume_km3` at the event level,
`volume_km3` at the daily level, both computed as area x vertical layer
thickness); `voxel` gains a `depth` column. For an ordinary 3D
`clim_file`, none of this changes: connectivity stays 3-axis (lon, lat,
time) and no depth columns are added.

## Ties and CPU type

Blobs with exactly equal values of `rankBy` (common for small blobs) are
ranked in blob order. On CPUs whose compilers fuse multiply-adds by
default (e.g. ARM, including Apple Silicon) the accumulated metrics can
differ from x86-64 in the last bit, so tied blobs can be ordered and
ranked differently, and `topN` can select a different blob at a tie.
Values agree to about 1e-15 relative.

## See also

[`ts2clm3`](https://robwschlegel.github.io/heatwave3/reference/ts2clm3.md)
for computing the climatology,
[`detect_event3`](https://robwschlegel.github.io/heatwave3/reference/detect_event3.md)
for per-pixel event detection (different from spatial blob detection).

## Examples

``` r
# \donttest{
sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
stem <- tempfile()

# Step 1: compute climatology
ts2clm3(sst_file, name = stem,
        climatologyPeriod = c("1982-01-01", "2011-12-31"),
        n_threads = 2)
#> Reading SST data from /Users/rschlegel/Library/R/arm64/4.6/library/heatwave3/extdata/sst_test.nc...
#> Grid: 2 lon x 3 lat x 14276 time = 6 pixels
#> Computing climatology with 2 thread(s)...
#> 
  6/6 pixels (100%)
#> Writing climatology to /var/folders/w8/8kw9l4sx1v57nx6l9s0rrklw0000gp/T//Rtmp4xYh4P/file9c7e2d5a4bfb_clim.nc...
#> Done.
#> 
#> ------------------------------------------------------------------
#> Climatology written to: /var/folders/w8/8kw9l4sx1v57nx6l9s0rrklw0000gp/T//Rtmp4xYh4P/file9c7e2d5a4bfb_clim.nc
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
#> Examine with  hw3_export("/var/folders/w8/8kw9l4sx1v57nx6l9s0rrklw0000gp/T//Rtmp4xYh4P/file9c7e2d5a4bfb_clim.nc", n = 20)
#> or export with hw3_export("/var/folders/w8/8kw9l4sx1v57nx6l9s0rrklw0000gp/T//Rtmp4xYh4P/file9c7e2d5a4bfb_clim.nc", file_out = "out.csv")  (.csv/.rds/.parquet)
#> ------------------------------------------------------------------

# Step 2: detect spatial blobs (heatwaves)
blobs <- detect_blob3(
  sst_file  = sst_file,
  clim_file = paste0(stem, "_clim.nc"),
  topN      = 10,
  return    = c("event", "daily", "voxel")
)

head(blobs$event)   # summary per blob
#>   blob duration date_start   date_end  date_peak cumI_km2_day peakArea_km2
#> 1   42       29 1986-02-25 1986-03-25 1986-02-25     8953.650     153.9682
#> 2  549       28 2018-09-09 2018-10-06 2018-09-10     8799.842     153.9682
#> 3  465       36 2013-09-03 2013-10-08 2013-09-04     7928.461     153.9682
#> 4   10       22 1982-11-06 1982-11-27 1982-11-11     6898.510     153.9682
#> 5  570       23 2019-10-19 2019-11-10 2019-10-21     6753.086     153.9682
#> 6   54       16 1987-01-21 1987-02-05 1987-01-21     6076.573     153.9682
#>   meanArea_km2 centroid_lon centroid_lat rank
#> 1     117.7061     26.54718    -34.06220    1
#> 2     145.7194     26.54984    -34.07530    2
#> 3     116.9140     26.54817    -34.06614    3
#> 4     130.6397     26.54911    -34.07498    4
#> 5     121.6198     26.54794    -34.06994    5
#> 6     144.3490     26.54944    -34.07276    6
head(blobs$daily)   # daily progression
#>   blob t_idx  area_km2 mean_delta max_delta centroid_lon centroid_lat lon_min
#> 1   10   310  51.29244   1.661800    1.7342       26.550    -34.12500  26.525
#> 2   10   311 102.61519   1.996661    2.3603       26.550    -34.09999  26.525
#> 3   10   312 102.61519   1.978866    2.2918       26.550    -34.09999  26.525
#> 4   10   313 128.29170   1.906916    2.2238       26.545    -34.08498  26.525
#> 5   10   314 102.61519   1.755398    1.8569       26.550    -34.09999  26.525
#> 6   10   315 153.96821   1.674754    1.7319       26.550    -34.07498  26.525
#>   lon_max lat_min lat_max n_cells       date
#> 1  26.575 -34.125 -34.125       2 1982-11-06
#> 2  26.575 -34.125 -34.075       4 1982-11-07
#> 3  26.575 -34.125 -34.075       4 1982-11-08
#> 4  26.575 -34.125 -34.025       5 1982-11-09
#> 5  26.575 -34.125 -34.075       4 1982-11-10
#> 6  26.575 -34.125 -34.025       6 1982-11-11
head(blobs$voxel)   # full spatial footprint
#>   blob    lon     lat t_idx       date  delta
#> 1   10 26.525 -34.125   310 1982-11-06 1.5894
#> 2   10 26.575 -34.125   310 1982-11-06 1.7342
#> 3   10 26.525 -34.125   311 1982-11-07 2.3603
#> 4   10 26.575 -34.125   311 1982-11-07 2.0639
#> 5   10 26.525 -34.075   311 1982-11-07 1.9020
#> 6   10 26.575 -34.075   311 1982-11-07 1.6607

# Step 2b: detect cold-spell blobs (use pctile = 10 climatology)
stem_cold <- tempfile()
ts2clm3(sst_file, name = stem_cold,
        climatologyPeriod = c("1982-01-01", "2011-12-31"),
        pctile = 10, n_threads = 2)
#> Reading SST data from /Users/rschlegel/Library/R/arm64/4.6/library/heatwave3/extdata/sst_test.nc...
#> Grid: 2 lon x 3 lat x 14276 time = 6 pixels
#> Computing climatology with 2 thread(s)...
#> 
  6/6 pixels (100%)
#> Writing climatology to /var/folders/w8/8kw9l4sx1v57nx6l9s0rrklw0000gp/T//Rtmp4xYh4P/file9c7e3326300_clim.nc...
#> Done.
#> 
#> ------------------------------------------------------------------
#> Climatology written to: /var/folders/w8/8kw9l4sx1v57nx6l9s0rrklw0000gp/T//Rtmp4xYh4P/file9c7e3326300_clim.nc
#> Rows (long format): 2,196   grid: 2 lon x 3 lat
#> 
#> Head:
#>      lon     lat doy     seas   thresh
#> 1 26.525 -34.125   1 294.4175 292.6069
#> 2 26.525 -34.125   2 294.4615 292.6500
#> 3 26.525 -34.125   3 294.5055 292.6917
#> 4 26.525 -34.125   4 294.5491 292.7345
#> 5 26.525 -34.125   5 294.5922 292.7736
#> 
#> Tail:
#>         lon     lat doy     seas   thresh
#> 2192 26.575 -34.025 362 293.5284 291.6923
#> 2193 26.575 -34.025 363 293.5648 291.7211
#> 2194 26.575 -34.025 364 293.6041 291.7575
#> 2195 26.575 -34.025 365 293.6456 291.8012
#> 2196 26.575 -34.025 366 293.6883 291.8445
#> 
#> Summary:
#>   ocean pixels (valid climatology): 6
#>   seas:   291.1 to 295.6
#>   thresh: 289.7 to 293.5
#> 
#> Examine with  hw3_export("/var/folders/w8/8kw9l4sx1v57nx6l9s0rrklw0000gp/T//Rtmp4xYh4P/file9c7e3326300_clim.nc", n = 20)
#> or export with hw3_export("/var/folders/w8/8kw9l4sx1v57nx6l9s0rrklw0000gp/T//Rtmp4xYh4P/file9c7e3326300_clim.nc", file_out = "out.csv")  (.csv/.rds/.parquet)
#> ------------------------------------------------------------------

mcs <- detect_blob3(
  sst_file   = sst_file,
  clim_file  = paste0(stem_cold, "_clim.nc"),
  coldSpells = TRUE,
  topN       = 6
)

# Step 2c: depth-resolved blobs, connected through the water column too
glorys_file <- system.file("extdata/glorys_depth_test.nc", package = "heatwave3")
stem_4d <- tempfile()
ts2clm3(glorys_file, name = stem_4d,
        climatologyPeriod = c("1993-01-01", "2019-12-31"),
        depth_range = c(0, 10), n_threads = 2)
#> Reading SST data from /Users/rschlegel/Library/R/arm64/4.6/library/heatwave3/extdata/glorys_depth_test.nc...
#> Grid: 2 lon x 2 lat x 5 depth x 11688 time = 20 pixels
#> Computing climatology with 2 thread(s)...
#> 
  20/20 pixels (100%)
#> Writing climatology to /var/folders/w8/8kw9l4sx1v57nx6l9s0rrklw0000gp/T//Rtmp4xYh4P/file9c7e610196a3_clim.nc...
#> Done.
#> 
#> ------------------------------------------------------------------
#> Climatology written to: /var/folders/w8/8kw9l4sx1v57nx6l9s0rrklw0000gp/T//Rtmp4xYh4P/file9c7e610196a3_clim.nc
#> Rows (long format): 7,320   grid: 2 lon x 2 lat
#> 
#> Head:
#>        lon       lat doy    seas  thresh    depth
#> 1 25.16667 -34.91667   1 24.4399 25.7963 0.494025
#> 2 25.16667 -34.91667   2 24.4776 25.8311 0.494025
#> 3 25.16667 -34.91667   3 24.5145 25.8679 0.494025
#> 4 25.16667 -34.91667   4 24.5506 25.9044 0.494025
#> 5 25.16667 -34.91667   5 24.5857 25.9402 0.494025
#> 
#> Tail:
#>        lon       lat doy    seas  thresh    depth
#> 7316 25.25 -34.83333 362 24.0061 25.4965 5.078224
#> 7317 25.25 -34.83333 363 24.0422 25.5364 5.078224
#> 7318 25.25 -34.83333 364 24.0788 25.5742 5.078224
#> 7319 25.25 -34.83333 365 24.1164 25.6112 5.078224
#> 7320 25.25 -34.83333 366 24.1542 25.6515 5.078224
#> 
#> Summary:
#>   ocean pixels (valid climatology): 4
#>   seas:   20.98 to 25.55
#>   thresh:    22 to 26.92
#> 
#> Examine with  hw3_export("/var/folders/w8/8kw9l4sx1v57nx6l9s0rrklw0000gp/T//Rtmp4xYh4P/file9c7e610196a3_clim.nc", n = 20)
#> or export with hw3_export("/var/folders/w8/8kw9l4sx1v57nx6l9s0rrklw0000gp/T//Rtmp4xYh4P/file9c7e610196a3_clim.nc", file_out = "out.csv")  (.csv/.rds/.parquet)
#> ------------------------------------------------------------------

blobs_4d <- detect_blob3(
  sst_file  = glorys_file,
  clim_file = paste0(stem_4d, "_clim.nc"),
  topN      = 10
)
head(blobs_4d$event)  # gains depth_min_m/depth_max_m/totalVolume_km3/...
#>   blob duration date_start   date_end  date_peak cumI_km2_day peakArea_km2
#> 1   67       22 2000-01-25 2000-02-15 2000-01-25     61533.09     1412.072
#> 2   54       28 1999-05-31 1999-06-27 1999-06-02     49120.23     1412.072
#> 3   85       24 2001-04-14 2001-05-07 2001-04-17     49058.69     1412.072
#> 4  294       21 2015-12-08 2015-12-28 2015-12-08     47262.29     1412.072
#> 5  290       21 2015-10-07 2015-10-27 2015-10-08     43830.41     1412.072
#> 6  390       23 2023-05-03 2023-05-25 2023-05-03     42414.75     1412.072
#>   meanArea_km2 centroid_lon centroid_lat depth_min_m depth_max_m
#> 1     1412.072     25.20833    -34.87498    0.494025    5.078224
#> 2     1222.987     25.20790    -34.87283    0.494025    5.078224
#> 3     1235.518     25.20933    -34.87796    0.494025    5.078224
#> 4     1334.739     25.20781    -34.87529    0.494025    5.078224
#> 5     1361.615     25.20885    -34.87652    0.494025    5.078224
#> 6     1273.880     25.20884    -34.87849    0.494025    5.078224
#>   totalVolume_km3 peakVolume_km3 meanVolume_km3 depthOfPeakIntensity_m rank
#> 1        35.64613       1.620279       1.620279               0.494025    1
#> 2        39.29278       1.620279       1.403314               5.078224    2
#> 3        34.02462       1.620279       1.417692               5.078224    3
#> 4        32.14827       1.620279       1.530870               0.494025    4
#> 5        32.81003       1.620279       1.562382               1.541375    5
#> 6        33.61934       1.620279       1.461711               5.078224    6
# }
```
