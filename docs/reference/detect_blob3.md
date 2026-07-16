<div id="main" class="col-md-9" role="main">

# Detect spatially connected marine heatwave blobs

<div class="ref-description section level2">

Identifies spatially contiguous marine heatwave (or cold-spell) events
("blobs") by performing connected-component labelling on the exceedance
mask across longitude, latitude, (optionally) depth, and time.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

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
  skip_bad_files = FALSE
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   sst\_file:

    Path to the SST NetCDF file (the same file passed to `ts2clm3`).

-   clim\_file:

    Path to the climatology NetCDF produced by `ts2clm3`. The spatial
    (and, if depth-resolved, vertical) grid of the climatology
    determines the region analysed.

-   var\_name:

    Name of the SST variable in `sst_file`. If `NULL`, auto-detected
    from CF attributes.

-   wrapDateline:

    Logical. Wrap the longitude axis so that blobs can connect across
    the antimeridian? Default `FALSE`. Depth is never wrapped (not
    periodic).

-   minVoxels:

    Minimum number of voxels (space-time, or space-depth-time for
    depth-resolved input) for a blob to be retained. Default `1`. Set
    higher (for example `200`) to filter out small, short-lived events.

-   topN:

    Return only the top N blobs ranked by `rankBy`. Default `NULL`
    (return all).

-   rankBy:

    Character. Metric to rank blobs by. Default `"cumI_km2_day"`
    (cumulative intensity-weighted area). Other options include
    `"peakArea_km2"`, `"duration"`, `"meanArea_km2"`, and, for
    depth-resolved input only, `"totalVolume_km3"`, `"peakVolume_km3"`,
    `"meanVolume_km3"`.

-   coldSpells:

    Logical. If `TRUE`, detect cold-spell blobs (`temp < thresh`).
    Requires a climatology computed with a low percentile (for example
    `pctile = 10` in `ts2clm3`). Default `FALSE`.

-   return:

    Character vector specifying which components to include in the
    output. One or more of:

    `"event"`

    :   One row per blob with summary metrics (duration, peak area,
        cumulative intensity, centroid, etc.)

    `"daily"`

    :   One row per (blob, date) with daily area, mean/max anomaly,
        centroid, and bounding box

    `"voxel"`

    :   One row per (blob, lon, lat, date) (plus depth, if
        depth-resolved), the full footprint. Required for spatial
        footprint maps and persistence analysis. Can be large.

    Default `c("event", "daily")`.

-   skip\_bad\_files:

    Logical. For multi-file SST inputs, skip unreadable files or files
    with mismatched grids instead of failing. Default `FALSE`.

</div>

<div class="section level2">

## Value

A named list containing the requested components (`event`, `daily`,
and/or `voxel` data.frames).

</div>

<div class="section level2">

## How it works

This function does **not** read a pre-computed event file. Instead, it
takes the raw SST data and the climatology, and internally constructs a
boolean exceedance mask by comparing each pixel's daily SST against its
threshold climatology (`temp > thresh` for heatwaves, or `temp < thresh`
for cold-spells). Connected voxels in this mask are then labelled as
coherent spatial (or spatio-vertical) events using a union-find
algorithm.

This is the same approach used by `heatwaveR::detect_blob3()`. The
reason an event file from `detect_event3` cannot be used directly is
that the event file stores per-pixel event *summaries* (start/end dates,
intensity metrics), not the daily per-pixel boolean mask needed for
spatial connectivity analysis.

</div>

<div class="section level2">

## Depth connectivity

When `clim_file` is depth-resolved (from `ts2clm3(depth_range = ...)`),
`detect_blob3()` detects this automatically – no separate argument is
needed, the same pattern used by `detect_event3` and `category_daily3`.
Voxels are then connected across 4 axes (lon, lat, depth, time) instead
of 3, using face-adjacency by depth *index*: a voxel connects to the
level immediately above/below it, whatever that level's actual depth in
metres happens to be, with no minimum/maximum depth gap or distance
threshold. `event`/`daily` gain `depth_min_m`, `depth_max_m`,
`depthOfPeakIntensity_m`, and volumetric analogues of the area metrics
(`totalVolume_km3`/`peakVolume_km3`/`meanVolume_km3` at the event level,
`volume_km3` at the daily level, both computed as area x vertical layer
thickness); `voxel` gains a `depth` column. For an ordinary 3D
`clim_file`, none of this changes: connectivity stays 3-axis (lon, lat,
time) and no depth columns are added.

</div>

<div class="section level2">

## See also

<div class="dont-index">

`ts2clm3` for computing the climatology, `detect_event3` for per-pixel
event detection (different from spatial blob detection).

</div>

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
if (FALSE) { # \dontrun{
sst_file <- "/path/to/sst.nc"
stem <- tempfile()

# Step 1: compute climatology
ts2clm3(sst_file, name = stem,
        climatologyPeriod = c("1982-01-01", "2011-12-31"),
        n_threads = 4)

# Step 2: detect spatial blobs (heatwaves)
blobs <- detect_blob3(
  sst_file  = sst_file,
  clim_file = paste0(stem, "_clim.nc"),
  minVoxels = 200,
  topN      = 10,
  return    = c("event", "daily", "voxel")
)

head(blobs$event)   # summary per blob
head(blobs$daily)   # daily progression
head(blobs$voxel)   # full spatial footprint

# Step 2b: detect cold-spell blobs (use pctile = 10 climatology)
stem_cold <- tempfile()
ts2clm3(sst_file, name = stem_cold,
        climatologyPeriod = c("1982-01-01", "2011-12-31"),
        pctile = 10, n_threads = 4)

mcs <- detect_blob3(
  sst_file   = sst_file,
  clim_file  = paste0(stem_cold, "_clim.nc"),
  coldSpells = TRUE,
  minVoxels  = 200,
  topN       = 6
)

# Step 2c: depth-resolved blobs, connected through the water column too
stem_4d <- tempfile()
ts2clm3(sst_file, name = stem_4d,
        climatologyPeriod = c("1982-01-01", "2011-12-31"),
        depth_range = c(0, 100), n_threads = 4)

blobs_4d <- detect_blob3(
  sst_file  = sst_file,
  clim_file = paste0(stem_4d, "_clim.nc"),
  minVoxels = 200,
  topN      = 10
)
head(blobs_4d$event)  # gains depth_min_m/depth_max_m/totalVolume_km3/...
} # }
```

</div>

</div>

</div>
