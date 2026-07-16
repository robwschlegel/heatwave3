# Additional hw3_export() / .hw3_console_summary() branches not exercised by
# test-hw3_export.R or test-ts2clm3.R: missing-file error, the type= override,
# RDS export (whole and subset), an out-of-range subset that returns zero
# rows, events-product lon/lat/depth_range/vars filtering, invalid time_range
# strings, multi-chunk CSV writes, and the non-quiet console summary for the
# events/daily/protoevents products (only ever exercised with quiet = TRUE
# elsewhere).

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

test_that("hw3_export errors for a missing file", {
  expect_error(hw3_export("nonexistent_product.nc"), "does not exist")
})

test_that("type= overrides auto-detection and agrees with it", {
  skip_if(system.file("extdata/sst_test.nc", package = "heatwave3") == "", "no fixture")
  f <- build_all()

  expect_equal(hw3_export(f$clim, n = 1), hw3_export(f$clim, n = 1, type = "clim"))
  expect_equal(hw3_export(f$events, n = 1), hw3_export(f$events, n = 1, type = "event"))
  expect_equal(hw3_export(f$daily, n = 1), hw3_export(f$daily, n = 1, type = "daily"))
  expect_equal(hw3_export(f$proto, n = 1), hw3_export(f$proto, n = 1, type = "protoevents"))
})

test_that("hw3_export writes RDS, whole and subset", {
  skip_if(system.file("extdata/sst_test.nc", package = "heatwave3") == "", "no fixture")
  f <- build_all()

  rds_whole <- tempfile(fileext = ".rds")
  hw3_export(f$clim, file_out = rds_whole)
  expect_equal(readRDS(rds_whole), hw3_export(f$clim))

  rds_sub <- tempfile(fileext = ".rds")
  hw3_export(f$clim, vars = "seas", file_out = rds_sub)
  expect_named(readRDS(rds_sub), c("lon", "lat", "doy", "seas"))
})

test_that("hw3_export validates chunk_size and the save_file extension", {
  skip_if(system.file("extdata/sst_test.nc", package = "heatwave3") == "", "no fixture")
  f <- build_all()

  expect_error(
    hw3_export(f$clim, file_out = tempfile(fileext = ".csv"), chunk_size = 0L),
    "chunk_size must be a positive integer"
  )
  expect_error(
    hw3_export(f$clim, file_out = tempfile(fileext = ".txt")),
    "\\.csv, \\.rds, or \\.parquet"
  )
})

test_that("hw3_export writes multiple CSV chunks for clim, daily, and events products", {
  skip_if(system.file("extdata/sst_test.nc", package = "heatwave3") == "", "no fixture")
  f <- build_all()

  for (product in c("clim", "daily", "events")) {
    whole <- hw3_export(f[[product]])
    csv <- tempfile(fileext = ".csv")
    hw3_export(f[[product]], file_out = csv, chunk_size = 100L)
    written <- utils::read.csv(csv)
    expect_equal(nrow(written), nrow(whole))
  }
})

test_that("hw3_export subset returns an empty data.frame for an out-of-range window", {
  skip_if(system.file("extdata/sst_test.nc", package = "heatwave3") == "", "no fixture")
  f <- build_all()

  empty <- hw3_export(f$clim, lon_range = c(100, 101))
  expect_s3_class(empty, "data.frame")
  expect_equal(nrow(empty), 0L)
  expect_named(empty, c("lon", "lat", "doy"))
})

test_that("hw3_export events subset supports lon/lat and depth_range (with a warning when absent)", {
  skip_if(system.file("extdata/sst_test.nc", package = "heatwave3") == "", "no fixture")
  f <- build_all()
  ev <- hw3_export(f$events)
  one_lon <- ev$lon[1]; one_lat <- ev$lat[1]

  sub <- hw3_export(f$events, lon_range = c(one_lon, one_lon), lat_range = c(one_lat, one_lat))
  expect_true(all(sub$lon == one_lon & sub$lat == one_lat))
  expect_lt(nrow(sub), nrow(ev))

  expect_warning(
    hw3_export(f$events, depth_range = c(0, 10)),
    "depth_range is ignored"
  )
})

test_that("hw3_export rejects an unparseable time_range for daily and events products", {
  skip_if(system.file("extdata/sst_test.nc", package = "heatwave3") == "", "no fixture")
  f <- build_all()

  expect_error(hw3_export(f$daily, time_range = c(NA_character_, "2010-01-01")),
               "must be dates")
  expect_error(hw3_export(f$events, time_range = c(NA_character_, "2010-01-01")),
               "must be dates")
})

test_that("hw3_export console summary (quiet = FALSE) prints for events, daily, and protoevents", {
  skip_if(system.file("extdata/sst_test.nc", package = "heatwave3") == "", "no fixture")
  sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
  stem <- tempfile()
  ts2clm3(sst_file, name = stem, climatologyPeriod = cp, quiet = TRUE)

  out_ev <- utils::capture.output(
    detect_event3(sst_file, name = stem, category = TRUE, daily = "also")
  )
  joined_ev <- paste(out_ev, collapse = "\n")
  expect_match(joined_ev, "Events written to")
  expect_match(joined_ev, "pixels with events:")
  expect_match(joined_ev, "duration \\(days\\):")
  expect_match(joined_ev, "category: ")
  expect_match(joined_ev, "Daily series written to")
  expect_match(joined_ev, "event-days:")

  out_proto <- utils::capture.output(
    detect_event3(sst_file, name = tempfile(), clim_file = paste0(stem, "_clim.nc"),
                 protoEvent = TRUE)
  )
  joined_proto <- paste(out_proto, collapse = "\n")
  expect_match(joined_proto, "Proto-events written to")
  expect_match(joined_proto, "threshCriterion days:")
})
