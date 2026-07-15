# Convenience wrapper for marine heatwave event categories

`category3()` is a convenience wrapper for adding or reading Hobday et
al. (2018) categories (I Moderate, II Strong, III Severe, IV Extreme)
for an event NetCDF file. New analyses should usually set
`category = TRUE` in
[`detect_event3`](https://robwschlegel.github.io/heatwave3/index.html/reference/detect_event3.md)
or
[`detect3`](https://robwschlegel.github.io/heatwave3/index.html/reference/detect3.md),
which computes the same categories during event detection and avoids a
second pass over the files.

## Usage

``` r
category3(
  event_file,
  clim_file = NULL,
  hemisphere = "south",
  roundVal = 4L,
  S = NULL
)
```

## Arguments

- event_file:

  Path to the event NetCDF file from
  [`detect_event3`](https://robwschlegel.github.io/heatwave3/index.html/reference/detect_event3.md).

- clim_file:

  Path to the climatology NetCDF from
  [`ts2clm3`](https://robwschlegel.github.io/heatwave3/index.html/reference/ts2clm3.md).
  Only required when the event file does not contain pre-computed
  categories.

- hemisphere:

  Character. Either `"south"` (default, austral: DJF = Summer) or
  `"north"` (boreal: DJF = Winter).

- roundVal:

  Decimal places for rounding. Default `4`.

- S:

  Deprecated. Use `hemisphere` instead.

## Value

A data.frame with columns: event_no, lon, lat, peak_date, category,
intensity_max, duration, p_moderate, p_strong, p_severe, p_extreme,
season.

## Details

If the event file already contains categories, `category3()` reads those
values directly and `clim_file` may be omitted. If the event file does
not contain categories, `clim_file` is required so the category width
can be computed from the seasonal and threshold climatologies.

Works for both marine heatwaves and cold-spells: the category width is
`|seas - thresh|`, which is always positive regardless of whether the
threshold is above (heatwave) or below (cold-spell) the seasonal mean.

## Examples

``` r
# \donttest{
sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
clim_file <- tempfile(fileext = ".nc")
event_file <- tempfile(fileext = ".nc")

detect3(sst_file, clim_file, event_file,
        climatologyPeriod = c("1982-01-01", "2011-12-31"))
#> Reading SST data from /private/var/folders/3w/nmplbnm109b9903rx8z9q0kc0000gn/T/RtmppsylXt/temp_libpathec9b4aa2cdfa/heatwave3/extdata/sst_test.nc...
#> Error: NetCDF error in inq varid /var/folders/3w/nmplbnm109b9903rx8z9q0kc0000gn/T//Rtmp3f4glh/file17dc7622485d4.nc: NetCDF: Variable not found

cats <- category3(event_file, clim_file)
#> Error: Event file does not exist: /var/folders/3w/nmplbnm109b9903rx8z9q0kc0000gn/T//Rtmp3f4glh/file17dc7622485d4.nc
table(cats$category)
#> Error: object 'cats' not found
# }
```
