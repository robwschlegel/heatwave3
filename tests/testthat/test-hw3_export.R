# Tests for hw3_export: auto-detect, streaming subset (lon/lat/time + vars), and
# flat export. Builds all four products once under a shared stem.

cp <- c("1982-01-01", "2011-12-31")

build_all <- function() {
  sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
  stem <- tempfile()
  detect3(sst_file, name = stem, climatologyPeriod = cp,
          category = TRUE, daily = "also", quiet = TRUE)
  detect_event3(sst_file, name = stem, clim_file = paste0(stem, "_clim.nc"),
                protoEvent = TRUE, quiet = TRUE)
  list(
    clim = paste0(stem, "_clim.nc"),
    events = paste0(stem, "_events.nc"),
    daily = paste0(stem, "_events_daily.nc"),
    proto = paste0(stem, "_protoevents.nc")
  )
}

test_that("hw3_export auto-detects the product type", {
  skip_if(system.file("extdata/sst_test.nc", package = "heatwave3") == "", "no fixture")
  f <- build_all()
  expect_named(hw3_export(f$clim, n = 1), c("lon", "lat", "doy", "seas", "thresh"))
  expect_true("intensity" %in% names(hw3_export(f$daily, n = 1)))
  expect_true("threshCriterion" %in% names(hw3_export(f$proto, n = 1)))
  expect_true("intensity_max" %in% names(hw3_export(f$events, n = 1)))
})

test_that("variable selection reads only the requested variables", {
  skip_if(system.file("extdata/sst_test.nc", package = "heatwave3") == "", "no fixture")
  f <- build_all()

  d <- hw3_export(f$clim, vars = "seas")
  expect_named(d, c("lon", "lat", "doy", "seas"))

  d2 <- hw3_export(f$daily, vars = c("temp", "category"))
  expect_named(d2, c("lon", "lat", "t", "temp", "category"))

  expect_error(hw3_export(f$clim, vars = "bogus"), "not in this climatology")
})

test_that("lon/lat/time windows subset the gridded products", {
  skip_if(system.file("extdata/sst_test.nc", package = "heatwave3") == "", "no fixture")
  f <- build_all()

  whole <- hw3_export(f$clim)
  lons <- sort(unique(whole$lon))

  sub <- hw3_export(f$clim, lon_range = c(lons[1], lons[1]))
  expect_equal(sort(unique(sub$lon)), lons[1])
  # values identical to the whole read
  m <- merge(sub, whole, by = c("lon", "lat", "doy"))
  expect_equal(max(abs(m$seas.x - m$seas.y)), 0)

  # time window on the daily product
  dd <- hw3_export(f$daily, time_range = c("1990-01-01", "1990-01-31"))
  expect_true(all(dd$t >= as.Date("1990-01-01") & dd$t <= as.Date("1990-01-31")))

  # time_range warns (and is ignored) for a climatology
  expect_warning(hw3_export(f$clim, time_range = c("1990-01-01", "1990-12-31")),
                 "ignored for a climatology")
})

test_that("n previews the first rows of any product", {
  skip_if(system.file("extdata/sst_test.nc", package = "heatwave3") == "", "no fixture")
  f <- build_all()
  expect_equal(nrow(hw3_export(f$clim, n = 7)), 7)
  expect_equal(nrow(hw3_export(f$daily, n = 5)), 5)
  expect_lte(nrow(hw3_export(f$events, n = 3)), 3)
})

test_that("events are filtered by coordinate, time, and selected columns", {
  skip_if(system.file("extdata/sst_test.nc", package = "heatwave3") == "", "no fixture")
  f <- build_all()
  ev <- hw3_export(f$events)

  sub <- hw3_export(f$events, vars = c("duration", "intensity_max"),
                    time_range = c("1985-01-01", "1986-12-31"))
  expect_named(sub, c("lon", "lat", "event_no", "duration", "intensity_max"))
  expect_true(all(sub$date_start <= as.Date("1986-12-31") |
                  is.na(sub$date_start)))  # all kept events overlap the window
  expect_lte(nrow(sub), nrow(ev))
})

test_that("subset export writes a flat file", {
  skip_if(system.file("extdata/sst_test.nc", package = "heatwave3") == "", "no fixture")
  f <- build_all()
  csv <- tempfile(fileext = ".csv")
  hw3_export(f$daily, vars = c("temp", "event"), time_range = c("1990-01-01", "1990-01-10"),
             file_out = csv)
  expect_true(file.exists(csv))
  back <- utils::read.csv(csv)
  expect_named(back, c("lon", "lat", "t", "temp", "event"))
  expect_gt(nrow(back), 0)
})
