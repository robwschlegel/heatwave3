<div id="main" class="col-md-9" role="main">

# Convenience wrapper for marine heatwave event categories

<div class="ref-description section level2">

`category3()` is a convenience wrapper for adding or reading Hobday et
al. (2018) categories (I Moderate, II Strong, III Severe, IV Extreme)
for an event NetCDF file. New analyses should usually set
`category = TRUE` in `detect_event3` or `detect3`, which computes the
same categories during event detection and avoids a second pass over the
files.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
category3(
  event_file,
  clim_file = NULL,
  hemisphere = "south",
  roundVal = 4L,
  S = NULL
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   event\_file:

    Path to the event NetCDF file from `detect_event3`.

-   clim\_file:

    Path to the climatology NetCDF from `ts2clm3`. Only required when
    the event file does not contain pre-computed categories.

-   hemisphere:

    Character. Either `"south"` (default, austral: DJF = Summer) or
    `"north"` (boreal: DJF = Winter).

-   roundVal:

    Decimal places for rounding. Default `4`.

-   S:

    Deprecated. Use `hemisphere` instead.

</div>

<div class="section level2">

## Value

A data.frame with columns: event\_no, lon, lat, peak\_date, category,
intensity\_max, duration, p\_moderate, p\_strong, p\_severe, p\_extreme,
season.

</div>

<div class="section level2">

## Details

If the event file already contains categories, `category3()` reads those
values directly and `clim_file` may be omitted. If the event file does
not contain categories, `clim_file` is required so the category width
can be computed from the seasonal and threshold climatologies.

Works for both marine heatwaves and cold-spells: the category width is
`|seas - thresh|`, which is always positive regardless of whether the
threshold is above (heatwave) or below (cold-spell) the seasonal mean.

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
#> Writing climatology to /var/folders/3w/nmplbnm109b9903rx8z9q0kc0000gn/T//RtmpEGZ5mo/file719a614d755c.nc...
#> Done.
#> Reading climatology from /var/folders/3w/nmplbnm109b9903rx8z9q0kc0000gn/T//RtmpEGZ5mo/file719a614d755c.nc...
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
#> Writing events to /var/folders/3w/nmplbnm109b9903rx8z9q0kc0000gn/T//RtmpEGZ5mo/file719a24541cec.nc...
#> Done.

cats <- category3(event_file, clim_file)
table(cats$category)
#> < table of extent 0 >
# }
```

</div>

</div>

</div>
