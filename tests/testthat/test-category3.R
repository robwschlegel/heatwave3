test_that("category3 reads inline categories and validates missing climatology", {
  sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
  skip_if(sst_file == "", "sst_test.nc not available")

  cp <- c("1982-01-01", "2011-12-31")

  # events with inline categories
  stem <- tempfile()
  ts2clm3(sst_file, name = stem, climatologyPeriod = cp, quiet = TRUE)
  detect_event3(sst_file, name = stem, category = TRUE, quiet = TRUE)

  cats <- category3(paste0(stem, "_events.nc"))
  expect_s3_class(cats, "data.frame")
  expect_true("category" %in% names(cats))
  expect_gt(nrow(cats), 0)

  # plain events (no inline categories) require a climatology file
  stem2 <- tempfile()
  ts2clm3(sst_file, name = stem2, climatologyPeriod = cp, quiet = TRUE)
  detect_event3(sst_file, name = stem2, quiet = TRUE)
  ev2 <- paste0(stem2, "_events.nc")

  expect_error(category3(ev2), "clim_file is required")
  expect_s3_class(category3(ev2, paste0(stem2, "_clim.nc")), "data.frame")
})

test_that("category3 validates its file arguments", {
  sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
  skip_if(sst_file == "", "sst_test.nc not available")

  expect_error(category3("nonexistent_events.nc"), "Event file does not exist")

  stem <- tempfile()
  ts2clm3(sst_file, name = stem, climatologyPeriod = c("1982-01-01", "2011-12-31"),
          quiet = TRUE)
  detect_event3(sst_file, name = stem, quiet = TRUE)
  expect_error(
    category3(paste0(stem, "_events.nc"), "nonexistent_clim.nc"),
    "Climatology file does not exist"
  )
})

test_that("category3's deprecated S argument still maps to hemisphere", {
  sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
  skip_if(sst_file == "", "sst_test.nc not available")

  stem <- tempfile()
  ts2clm3(sst_file, name = stem, climatologyPeriod = c("1982-01-01", "2011-12-31"),
          quiet = TRUE)
  detect_event3(sst_file, name = stem, category = TRUE, quiet = TRUE)
  ev <- paste0(stem, "_events.nc")

  expect_warning(cats_s <- category3(ev, S = TRUE), "S is deprecated")
  expect_warning(cats_n <- category3(ev, S = FALSE), "S is deprecated")
  expect_equal(cats_s, category3(ev, hemisphere = "south"))
  expect_equal(cats_n, category3(ev, hemisphere = "north"))
})
