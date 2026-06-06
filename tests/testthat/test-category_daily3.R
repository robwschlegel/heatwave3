# category_daily3(): windowed daily categories from SST + clim + events,
# checked for internal consistency against the per-day product (which carries
# temp/seas/thresh) on the bundled fixture.

test_that("category_daily3 reproduces the daily-category logic for a window", {
  sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
  skip_if(sst_file == "", "sst_test.nc not available")

  stem <- tempfile()
  cp <- c("1982-01-01", "2011-12-31")
  ts2clm3(sst_file, name = stem, climatologyPeriod = cp, quiet = TRUE)
  detect_event3(sst_file, name = stem, daily = "also", quiet = TRUE)
  clim  <- paste0(stem, "_clim.nc")
  event <- paste0(stem, "_events.nc")
  daily <- paste0(stem, "_events_daily.nc")

  # A window near the end of the series (guaranteed to have data)
  dd <- hw3_export(daily, vars = "event", n = 1)  # cheap, just to confirm readable
  full_t <- range(hw3_export(daily, vars = "temp")$t)
  win <- c(format(full_t[2] - 120), format(full_t[2]))

  res <- category_daily3(sst_file, clim, event, time_range = win)

  # Full daily grid with the same columns as the detect_event3 daily product.
  expect_s3_class(res, "data.frame")
  expect_named(res, c("lon", "lat", "t", "temp", "seas", "thresh",
                      "intensity", "event", "event_no", "category"))
  expect_s3_class(res$t, "Date")
  expect_type(res$event, "logical")
  expect_true(all(res$t >= as.Date(win[1]) & res$t <= as.Date(win[2])))

  # MHW: category is 1-4 on event-member exceedance days, NA elsewhere.
  expect_true(all(is.na(res$category) | res$category %in% 1:4))
  expect_true(all(res$event[!is.na(res$category)]))       # category implies in-event
  expect_true(all(is.na(res$event_no) == !res$event))     # event_no set iff in-event

  # Join to the per-day product. detect_event3 ran on the full series here, so
  # the daily product carries the same full-record events: event membership must
  # match exactly. temp/seas/thresh agree (same SST/clim reads, float storage),
  # and intensity differs only by rounding (2 dp here vs 4 dp stored).
  dpro <- hw3_export(daily, time_range = win)
  m <- merge(res, dpro, by = c("t", "lon", "lat"))
  expect_equal(nrow(m), nrow(res))                        # same full grid
  expect_equal(m$event.x, m$event.y)                      # same event membership
  expect_equal(m$event_no.x, m$event_no.y)                # same event numbering
  expect_true(all(abs(m$temp.x - m$temp.y) < 1e-2, na.rm = TRUE))
  expect_true(all(abs(m$seas.x - m$seas.y) < 1e-2, na.rm = TRUE))
  expect_true(all(abs(m$thresh.x - m$thresh.y) < 1e-2, na.rm = TRUE))
  expect_true(all(abs(m$intensity.x - m$intensity.y) <= 0.01, na.rm = TRUE))

  # On event-member exceedance days the category equals the daily product's.
  ev_days <- m[m$event.x & !is.na(m$category.x), ]
  expect_equal(ev_days$category.x, ev_days$category.y)
})

test_that("category_daily3 validates inputs", {
  sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
  skip_if(sst_file == "", "sst_test.nc not available")
  stem <- tempfile()
  ts2clm3(sst_file, name = stem, climatologyPeriod = c("1982-01-01", "2011-12-31"),
          quiet = TRUE)
  detect_event3(sst_file, name = stem, quiet = TRUE)

  expect_error(
    category_daily3(sst_file, paste0(stem, "_clim.nc"), paste0(stem, "_events.nc"),
                    time_range = "2005-01-01"),
    "two dates")
  expect_error(
    category_daily3("nope.nc", paste0(stem, "_clim.nc"), paste0(stem, "_events.nc"),
                    time_range = c("2005-01-01", "2005-02-01")),
    "does not exist")
  expect_error(
    category_daily3(sst_file, paste0(stem, "_clim.nc"), paste0(stem, "_events.nc"),
                    time_range = c("2005-01-01", "2005-02-01"), lon_range = 26),
    "lon_range must be")
  expect_error(
    category_daily3(sst_file, paste0(stem, "_clim.nc"), paste0(stem, "_events.nc"),
                    time_range = c("2005-01-01", "2005-02-01"), lat_range = c(1, 2, 3)),
    "lat_range must be")
})

test_that("category_daily3 lon_range/lat_range subset equals the filtered full result", {
  sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
  skip_if(sst_file == "", "sst_test.nc not available")
  stem <- tempfile()
  ts2clm3(sst_file, name = stem, climatologyPeriod = c("1982-01-01", "2011-12-31"),
          quiet = TRUE)
  detect_event3(sst_file, name = stem, quiet = TRUE)
  clim <- paste0(stem, "_clim.nc"); evt <- paste0(stem, "_events.nc")
  win <- c("2015-01-01", "2015-12-31")
  cd <- hw3_read_clim_nc(clim)

  full <- category_daily3(sst_file, clim, evt, time_range = win)
  lon_r <- c(cd$lon[1] - 1e-4, cd$lon[1] + 1e-4)              # first longitude
  lat_r <- c(min(cd$lat[1:2]) - 1e-4, max(cd$lat[1:2]) + 1e-4) # first two latitudes
  sub <- category_daily3(sst_file, clim, evt, time_range = win,
                         lon_range = lon_r, lat_range = lat_r)

  # subset stays within the requested window
  expect_true(all(sub$lon >= lon_r[1] & sub$lon <= lon_r[2]))
  expect_true(all(sub$lat >= lat_r[1] & sub$lat <= lat_r[2]))

  # and is identical to the full result filtered to the same window
  ref <- full[full$lon >= lon_r[1] & full$lon <= lon_r[2] &
              full$lat >= lat_r[1] & full$lat <= lat_r[2], ]
  expect_equal(nrow(sub), nrow(ref))
  m <- merge(sub, ref, by = c("lon", "lat", "t"), suffixes = c(".s", ".r"))
  expect_equal(nrow(m), nrow(sub))                  # exact row-set match
  for (col in c("temp", "seas", "thresh", "intensity", "event_no", "category")) {
    expect_equal(m[[paste0(col, ".s")]], m[[paste0(col, ".r")]], info = col)
  }
})
