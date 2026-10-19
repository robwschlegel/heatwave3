# Detect marine heatwaves day by day in arrays

The in-memory counterpart of
[`detect_event3`](https://robwschlegel.github.io/heatwave3/reference/detect_event3.md)
with `daily = "only"`: the same detection, but on arrays instead of
files, and with the climatology already matched to each day. It returns
only per-day results, because event metrics such as duration or
cumulative intensity need the whole event, which a short series (a
forecast window, for example) may cut off.

## Usage

``` r
detect_event_array3(
  temp,
  seas,
  thresh,
  minDuration = 5,
  joinAcrossGaps = TRUE,
  maxGap = 2,
  coldSpells = FALSE,
  roundRes = 4,
  outputs = c("event", "intensity", "category"),
  n_threads = NULL
)
```

## Arguments

- temp:

  Numeric vector, matrix or array of daily temperatures with time on the
  **first** dimension (a vector is one series) and any number of further
  dimensions (pixels). The days must be consecutive, with no gaps.
  Missing days are `NA`.

- seas, thresh:

  Seasonal climatology and threshold for each day of `temp`, with the
  same dimensions. Match them to the dates beforehand, for example from
  a
  [`ts2clm3`](https://robwschlegel.github.io/heatwave3/reference/ts2clm3.md)
  climatology by day of year.

- minDuration:

  Minimum duration (days) for an event. Default `5`.

- joinAcrossGaps:

  Logical. Join events separated by short gaps? Default `TRUE`.

- maxGap:

  Maximum gap length (days) to join across. Default `2`.

- coldSpells:

  Logical. Detect cold-spells instead of heatwaves? Default `FALSE`.

- roundRes:

  Decimal places for `intensity`. Default `4`.

- outputs:

  Which results to return, in this order: any of `"event"`,
  `"event_no"`, `"intensity"` and `"category"`.

- n_threads:

  Number of threads for the per-pixel computation. Default `NULL`: the
  session default (see
  [`setHW3threads`](https://robwschlegel.github.io/heatwave3/reference/getHW3threads.md)).
  Use `1` when the calls themselves run in parallel.

## Value

One array per name in `outputs`, each with the dimensions (and dimnames)
of `temp`: the array itself when one output is requested, otherwise a
named list in the order of `outputs`.

- event:

  Logical: the day belongs to an event.

- event_no:

  Integer: the event's number within its pixel, counted from the start
  of this call (0 outside events). It is not an identifier across calls:
  the same event in two overlapping series can have different numbers.

- intensity:

  Numeric: `temp - seas` on every day with data, rounded to `roundRes`;
  `NA` where either is missing.

- category:

  Integer: the daily Hobday et al. (2018) category, 1 (I Moderate) to 4
  (IV Extreme), on event days beyond the threshold; 0 on all other days.

A pixel with no temperature or no `seas` at all gets `FALSE`, 0, `NA`
and 0. Results are identical to the per-day file of
[`detect_event3`](https://robwschlegel.github.io/heatwave3/reference/detect_event3.md)
for the same days and climatology.

## Details

Whether a day belongs to an event depends on the days around it, so
events are cut off at both ends of the series: a run of exceedances that
started before the first day is shorter here, and one that continues
past the last day may not reach `minDuration`. Prepend earlier days when
that matters. `intensity` does not depend on the other days, and
`category` only through event membership.

Time comes first because R stores arrays column by column, so each
pixel's series is then contiguous and is read without copying. Gridded
data from packages such as terra or stars can be passed as plain arrays:
for example `t(terra::values(r))` for a SpatRaster `r` with one layer
per day, or `aperm(x, c(3, 1, 2))` for a `[x, y, time]` array.

## See also

[`detect_event3`](https://robwschlegel.github.io/heatwave3/reference/detect_event3.md),
[`ts2clm3`](https://robwschlegel.github.io/heatwave3/reference/ts2clm3.md)

## Examples

``` r
# Two pixels, 30 days, a constant climatology
temp <- matrix(20, nrow = 30, ncol = 2)
temp[6:15, 1] <- 25
seas <- matrix(20, nrow = 30, ncol = 2)
thresh <- matrix(21, nrow = 30, ncol = 2)
res <- detect_event_array3(temp, seas, thresh)
which(res$event[, 1])
#>  [1]  6  7  8  9 10 11 12 13 14 15
table(res$category)
#> 
#>  0  4 
#> 50 10 
```
