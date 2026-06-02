test_that("detect3 passes through companion outputs and returns events", {
  sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
  skip_if(sst_file == "", "sst_test.nc not available")

  out_dir <- tempfile()
  dir.create(out_dir)
  clim_file <- file.path(out_dir, "clim.nc")
  event_file <- file.path(out_dir, "events.nc")
  clim_csv <- file.path(out_dir, "clim.csv")
  event_csv <- file.path(out_dir, "events.csv")

  event_df <- detect3(
    sst_file,
    clim_file,
    event_file,
    c("1982-01-01", "2011-12-31"),
    save_file_clim = clim_csv,
    save_file_event = event_csv,
    detrend = TRUE,
    category = TRUE,
    return_df = TRUE
  )

  expect_s3_class(event_df, "data.frame")
  expect_true(file.exists(clim_csv))
  expect_true(file.exists(event_csv))
  expect_true("category" %in% names(event_df))
})
