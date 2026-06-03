test_that("detect3 writes clim + events from one name stem and returns their paths", {
  sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
  skip_if(sst_file == "", "sst_test.nc not available")

  stem <- tempfile()
  out <- detect3(
    sst_file,
    name = stem,
    climatologyPeriod = c("1982-01-01", "2011-12-31"),
    detrend = TRUE,
    category = TRUE,
    daily = "also",
    quiet = TRUE
  )

  expect_type(out, "character")
  expect_named(out, c("clim", "events", "daily"))
  expect_equal(unname(out["clim"]), paste0(stem, "_clim.nc"))
  expect_equal(unname(out["events"]), paste0(stem, "_events.nc"))
  expect_equal(unname(out["daily"]), paste0(stem, "_events_daily.nc"))
  expect_true(all(file.exists(out)))

  # category is carried into the events product
  ev <- hw3_export(out[["events"]])
  expect_true("category" %in% names(ev))
})

test_that("detect3 protoEvent writes only clim + protoevents", {
  sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
  skip_if(sst_file == "", "sst_test.nc not available")

  stem <- tempfile()
  out <- detect3(sst_file, name = stem,
                 climatologyPeriod = c("1982-01-01", "2011-12-31"),
                 protoEvent = TRUE, quiet = TRUE)

  expect_named(out, c("clim", "protoevents"))
  expect_true(file.exists(paste0(stem, "_protoevents.nc")))
  expect_false(file.exists(paste0(stem, "_events.nc")))
})
