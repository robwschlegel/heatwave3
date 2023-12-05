# Tests for the primary function: detect3()

test_that("Error messages return as expected", {
  sys_file <- system.file("extdata/oisst_short.nc", package = "heatwave3")
  expect_error(detect3())
  expect_error(detect3(file_in = sys_file))
  expect_error(detect3(file_in = system.file("extdata/oisst_short.nc", package = "heatwave3"),
                       clim_period = "1900-01-01"))
  expect_error(detect3(file_in = system.file("extdata/oisst_short.nc", package = "heatwave3"),
                       clim_period = c("1982-01-01", "2011-12-31")))
  expect_error(detect3(file_in = system.file("extdata/oisst_short.nc", package = "heatwave3"),
                       clim_period = c("1982-01-01", "2011-12-31"), return_type = "banana"))
  expect_error(detect3(file_in = system.file("extdata/oisst_short.nc", package = "heatwave3"),
                       clim_period = c("1982-01-01", "2011-12-31"), save_to_file = "banana"))
  expect_error(detect3(file_in = system.file("extdata/oisst_short.nc", package = "heatwave3"),
                       clim_period = c("1982-01-01", "2011-12-31"),
                       return_type = "rast", save_to_file = "csv"))
})

test_that("Arguments can be passed to ts2clm and detect_event via the ...", {
  sys_file <- system.file("extdata/oisst_short.nc", package = "heatwave3")

  # Just to test that raster output works
  res_rast <- detect3(file_in = sys_file, clim_period = c("1982-01-01", "2011-12-31"), return_type = "rast")
  expect_s4_class(res_rast, "SpatRasterDataset")

  # Compare different argument outputs
  res_base <- detect3(file_in = sys_file, clim_period = c("1982-01-01", "2011-12-31"), return_type = "df")
  res_99 <- detect3(file_in = sys_file, clim_period = c("1982-01-01", "2011-12-31"), return_type = "df", pctile = 99)
  res_MCS <- detect3(file_in = sys_file, clim_period = c("1982-01-01", "2011-12-31"), return_type = "df",
                     pctile = 10, coldSpells = TRUE)
  expect_s3_class(res_base, "data.frame")
  expect_equal(ncol(res_base), 22)
  expect_equal(nrow(res_base), 340)
  expect_equal(nrow(res_99), 35)
})

test_that("Saving types works correctly", {
  sys_file <- system.file("extdata/oisst_short.nc", package = "heatwave3")
  detect3(file_in = sys_file, clim_period = c("1982-01-01", "2011-12-31"),
          save_to_file = "nc", file_out = paste0(tempdir(),"/test.nc"))
  detect3(file_in = sys_file, clim_period = c("1982-01-01", "2011-12-31"),
          save_to_file = "csv", file_out = paste0(tempdir(),"/test.csv"))
  expect_true(file.exists(paste0(tempdir(),"/test.nc")))
  expect_true(file.exists(paste0(tempdir(),"/test.csv")))
})
