# Tests for the primary function: detect3()

test_that("The function runs as expected", {
  res_df <- detect3(file_in = system.file("extdata/oisst_short.nc", package = "heatwave3"),
                    return_type = "df")
  res_rast <- detect3(file_in = system.file("extdata/oisst_short.nc", package = "heatwave3"),
                      return_type = "rast")
  expect_s3_class(res_df, "data.frame")
  expect_equal(ncol(res_df), 22)
  expect_equal(nrow(res_df), 340)
  expect_s4_class(res_rast, "SpatRasterDataset")
})

test_that("Error messages return as expected", {
  expect_error(detect3(file_in = system.file("extdata/oisst_short.nc", package = "heatwave3")))
  expect_error(detect3(file_in = system.file("extdata/oisst_short.nc", package = "heatwave3"),
                       return_type = "banana"))
  expect_error(detect3(file_in = system.file("extdata/oisst_short.nc", package = "heatwave3"),
                       save_to_file = "banana"))
  expect_error(detect3(file_in = system.file("extdata/oisst_short.nc", package = "heatwave3"),
                       return_type = "rast", save_to_file = "csv"))
})

test_that("Saving types works correctly", {
  detect3(file_in = system.file("extdata/oisst_short.nc", package = "heatwave3"),
                 save_to_file = "nc", file_out = paste0(tempdir(),"/test.nc"))
  detect3(file_in = system.file("extdata/oisst_short.nc", package = "heatwave3"),
          save_to_file = "csv", file_out = paste0(tempdir(),"/test.csv"))
  expect_true(file.exists(paste0(tempdir(),"/test.nc")))
  expect_true(file.exists(paste0(tempdir(),"/test.csv")))
})
