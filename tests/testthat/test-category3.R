test_that("category3 reads inline categories and validates missing climatology", {
  sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
  skip_if(sst_file == "", "sst_test.nc not available")

  out_dir <- tempfile()
  dir.create(out_dir)
  clim_file <- file.path(out_dir, "clim.nc")
  event_cat_file <- file.path(out_dir, "events_cat.nc")
  event_plain_file <- file.path(out_dir, "events_plain.nc")

  ts2clm3(
    file_in = sst_file,
    file_out = clim_file,
    climatologyPeriod = c("1982-01-01", "2011-12-31")
  )

  detect_event3(
    file_in = sst_file,
    clim_file = clim_file,
    file_out = event_cat_file,
    category = TRUE
  )
  cats <- category3(event_cat_file)
  expect_s3_class(cats, "data.frame")
  expect_true("category" %in% names(cats))
  expect_gt(nrow(cats), 0)

  detect_event3(
    file_in = sst_file,
    clim_file = clim_file,
    file_out = event_plain_file
  )
  expect_error(category3(event_plain_file), "clim_file is required")
  expect_s3_class(category3(event_plain_file, clim_file), "data.frame")
})
