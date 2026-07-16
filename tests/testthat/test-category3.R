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
  expect_false(anyNA(cats$category))

  # plain events (no inline categories) require a climatology file
  stem2 <- tempfile()
  ts2clm3(sst_file, name = stem2, climatologyPeriod = cp, quiet = TRUE)
  detect_event3(sst_file, name = stem2, quiet = TRUE)
  ev2 <- paste0(stem2, "_events.nc")

  expect_error(category3(ev2), "clim_file is required")

  # Regression test: category3(event_file, clim_file) must actually compute
  # categories from the climatology rather than silently returning NA. The
  # "category" NetCDF variable is always written by detect_event3() (filled
  # with 0 when category = FALSE), so has_precomputed can't rely on mere
  # presence -- it previously did, and always took the "already computed"
  # branch, returning NA for every event no matter what.
  cats2 <- category3(ev2, paste0(stem2, "_clim.nc"))
  expect_s3_class(cats2, "data.frame")
  expect_false(anyNA(cats2$category))
  expect_true(all(cats2$category %in% c("I Moderate", "II Strong", "III Severe", "IV Extreme")))

  # Same SST/climatology period -> the after-the-fact computation must agree
  # with categories computed inline during detection.
  ord <- order(cats$lon, cats$lat, cats$peak_date)
  ord2 <- order(cats2$lon, cats2$lat, cats2$peak_date)
  expect_equal(cats$category[ord], cats2$category[ord2])
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
