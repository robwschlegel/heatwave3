<div id="main" class="col-md-9" role="main">

# Daily marine heatwave / cold-spell categories for a date window

<div class="ref-description section level2">

Builds the per-pixel, per-day series for a window of dates directly from
an SST file, a heatwave3 climatology, and a heatwave3 events file, using
the C++ readers throughout. It returns the full daily grid for the
window with the same columns as the `detect_event3` `daily = "also"`
product (`lon`, `lat`, `t`, `temp`, `seas`, `thresh`, `intensity`,
`event`, `event_no`, `category`).

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
category_daily3(
  sst_file,
  clim_file,
  event_file,
  time_range,
  var_name = NULL,
  coldSpells = FALSE,
  ice_thresh = -1.7,
  roundRes = 2L,
  skip_bad_files = FALSE
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   sst\_file:

    SST source: a single multi-time NetCDF (for example a per-longitude
    OISST time series), **or** a directory / character vector of daily
    files (one time step each, e.g. daily-global OSTIA/GHRSST). Only
    files covering `time_range` should be passed for a directory/vector.
    The variable is auto-detected unless `var_name` is given.

-   clim\_file:

    Path to the climatology NetCDF from `ts2clm3`.

-   event\_file:

    Path to the events NetCDF from `detect_event3`.

-   time\_range:

    Character vector of length 2, `c(start, end)` dates, for example
    `c("2026-05-10", "2026-05-20")`.

-   var\_name:

    Name of the SST variable. If `NULL`, auto-detected.

-   coldSpells:

    Logical. Categorise marine cold-spells instead of heatwaves? Default
    `FALSE`.

-   ice\_thresh:

    Numeric. For cold-spells only, days whose threshold is below this
    value (in the SST's units) are assigned category 5 ("ice"),
    following the Marine Heatwave Tracker convention. Default `-1.7`.
    Ignored for heatwaves.

-   roundRes:

    Number of decimal places for the returned `intensity`. Default `2`,
    mirroring the Marine Heatwave Tracker's `load_sub_cat_clim()`
    (`round(temp - seas, 2)`). Pass `roundRes = 4` to match
    `detect_event3`'s `daily = "also"` product exactly.

-   skip\_bad\_files:

    Logical. When `sst_file` is a directory or vector, skip unreadable
    files or files with mismatched grids instead of failing. Default
    `FALSE`.

</div>

<div class="section level2">

## Value

A `data.frame` with one row per pixel per day in the window: `lon`,
`lat`, `t` (Date), `temp`, `seas`, `thresh`, `intensity`
(`temp - seas`), `event` (logical; is the day inside a detected event?),
`event_no` (integer, `NA` off-event), and `category` (integer on
event-member exceedance days; 1 = I Moderate ... 4 = IV Extreme, 5 = ice
for cold-spells; `NA` otherwise). Columns match the `detect_event3`
`daily = "also"` product. A `depth` column (metres) is added
automatically when `clim_file`/`event_file` are depth-resolved (from
`ts2clm3(depth_range = ...)`) – one row per pixel per depth level per
day, no separate argument needed here.

</div>

<div class="section level2">

## Details

Crucially, event membership is read from the events file rather than
re-detected on the window, so events that began before `time_range` are
still recognised. Re-running `detect_event3()` on a short recent window
instead would silently drop ongoing events whose in-window portion is
shorter than the minimum duration. Only the requested time window of the
SST file is read (a hyperslab), which makes it cheap to call repeatedly
for the most recent days, as a daily operational pipeline does.

The Marine Heatwave Tracker's `load_sub_cat_clim()` subset (event-member
days with a category) is recovered with
`subset(x, !is.na(category))[c("t", "lon", "lat", "event_no", "intensity", "category")]`.

The SST, climatology, and events files must share the same spatial grid
(as they do when produced from one heatwave3 run for a given region or
longitude). Event membership and `event_no` are taken from `event_file`,
so the events must already cover the requested window.

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
if (FALSE) { # \dontrun{
category_daily3(
  sst_file   = "oisst-avhrr-v02r01.ts.0001.nc",
  clim_file  = "MHW_0001_1982-2011_clim.nc",
  event_file = "MHW_0001_1982-2011_events.nc",
  time_range = c("2026-05-10", "2026-05-20")
)
} # }
```

</div>

</div>

</div>
