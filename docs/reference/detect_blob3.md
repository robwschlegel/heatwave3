# Detect spatially connected marine heatwave blobs

Identifies spatially contiguous marine heatwave (or cold-spell) events
("blobs") by performing 3D connected-component labelling on the
exceedance mask across longitude, latitude, and time.

## Usage

``` r
detect_blob3(
  sst_file,
  clim_file,
  var_name = NULL,
  connectivity = 6L,
  wrapDateline = FALSE,
  minVoxels = 1L,
  topN = NULL,
  rankBy = "cumI_km2_day",
  coldSpells = FALSE,
  return = c("event", "daily"),
  skip_bad_files = FALSE
)
```

## Arguments

- sst_file:

  Path to the SST NetCDF file (the same file passed to
  [`ts2clm3`](https://robwschlegel.github.io/heatwave3/index.html/reference/ts2clm3.md)).

- clim_file:

  Path to the climatology NetCDF produced by
  [`ts2clm3`](https://robwschlegel.github.io/heatwave3/index.html/reference/ts2clm3.md).
  The spatial grid of the climatology determines the region analysed.

- var_name:

  Name of the SST variable in `sst_file`. If `NULL`, auto-detected from
  CF attributes.

- connectivity:

  Integer. Voxel connectivity for labelling. Default `6` (face-adjacent
  in 3D: left/right, up/down, forward/backward in time).

- wrapDateline:

  Logical. Wrap the longitude axis so that blobs can connect across the
  antimeridian? Default `FALSE`.

- minVoxels:

  Minimum number of space-time voxels for a blob to be retained. Default
  `1`. Set higher (for example `200`) to filter out small, short-lived
  events.

- topN:

  Return only the top N blobs ranked by `rankBy`. Default `NULL` (return
  all).

- rankBy:

  Character. Metric to rank blobs by. Default `"cumI_km2_day"`
  (cumulative intensity-weighted area). Other options include
  `"peakArea_km2"`, `"duration"`, `"meanArea_km2"`.

- coldSpells:

  Logical. If `TRUE`, detect cold-spell blobs (`temp < thresh`).
  Requires a climatology computed with a low percentile (for example
  `pctile = 10` in
  [`ts2clm3`](https://robwschlegel.github.io/heatwave3/index.html/reference/ts2clm3.md)).
  Default `FALSE`.

- return:

  Character vector specifying which components to include in the output.
  One or more of:

  `"event"`

  :   One row per blob with summary metrics (duration, peak area,
      cumulative intensity, centroid, etc.)

  `"daily"`

  :   One row per (blob, date) with daily area, mean/max anomaly,
      centroid, and bounding box

  `"voxel"`

  :   One row per (blob, lon, lat, date), the full 3D footprint.
      Required for spatial footprint maps and persistence analysis. Can
      be large.

  Default `c("event", "daily")`.

- skip_bad_files:

  Logical. For multi-file SST inputs, skip unreadable files or files
  with mismatched grids instead of failing. Default `FALSE`.

## Value

A named list containing the requested components (`event`, `daily`,
and/or `voxel` data.frames).

## How it works

This function does **not** read a pre-computed event file. Instead, it
takes the raw SST data and the climatology, and internally constructs a
3D boolean exceedance mask by comparing each pixel's daily SST against
its threshold climatology (`temp > thresh` for heatwaves, or
`temp < thresh` for cold-spells). Connected voxels in this 3D mask are
then labelled as coherent spatial events using a union-find algorithm.

This is the same approach used by `heatwaveR::detect_blob3()`. The
reason an event file from
[`detect_event3`](https://robwschlegel.github.io/heatwave3/index.html/reference/detect_event3.md)
cannot be used directly is that the event file stores per-pixel event
*summaries* (start/end dates, intensity metrics), not the daily
per-pixel boolean mask needed for spatial connectivity analysis.

## See also

[`ts2clm3`](https://robwschlegel.github.io/heatwave3/index.html/reference/ts2clm3.md)
for computing the climatology,
[`detect_event3`](https://robwschlegel.github.io/heatwave3/index.html/reference/detect_event3.md)
for per-pixel event detection (different from spatial blob detection).

## Examples

``` r
if (FALSE) { # \dontrun{
sst_file <- "/path/to/sst.nc"
clim_file <- tempfile(fileext = ".nc")

# Step 1: compute climatology
ts2clm3(sst_file, clim_file,
        climatologyPeriod = c("1982-01-01", "2011-12-31"),
        n_threads = 4)

# Step 2: detect spatial blobs (heatwaves)
blobs <- detect_blob3(
  sst_file  = sst_file,
  clim_file = clim_file,
  minVoxels = 200,
  topN      = 10,
  return    = c("event", "daily", "voxel")
)

head(blobs$event)   # summary per blob
head(blobs$daily)   # daily progression
head(blobs$voxel)   # full spatial footprint

# Step 2b: detect cold-spell blobs (use pctile = 10 climatology)
clim_cold <- tempfile(fileext = ".nc")
ts2clm3(sst_file, clim_cold,
        climatologyPeriod = c("1982-01-01", "2011-12-31"),
        pctile = 10, n_threads = 4)

mcs <- detect_blob3(
  sst_file   = sst_file,
  clim_file  = clim_cold,
  coldSpells = TRUE,
  minVoxels  = 200,
  topN       = 6
)
} # }
```
