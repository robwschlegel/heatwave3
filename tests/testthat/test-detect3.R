# Tests for the primary function: detect3()

test_that("The function runs as expected", {
  res <- detect3(file_in = system.file("extdata/oisst_short.nc", package = "heatwave3"),
                 return_type = "df")
  expect_s3_class(res, "data.frame")
  expect_equal(ncol(res), 22)
  expect_equal(nrow(res), 340)
})

test_that("Error messages return as expected", {
  expect_error(detect3(file_in = system.file("extdata/oisst_short.nc", package = "heatwave3")),
               "No output selected for the function.")
  expect_error(detect3(file_in = system.file("extdata/oisst_short.nc", package = "heatwave3"),
                       return_type = "banana"), "Invalid return_type.")
  expect_error(detect3(file_in = system.file("extdata/oisst_short.nc", package = "heatwave3"),
                       return_type = "df", save_to_file = "banana"), "Invalid saving option.")
  expect_error(detect3(file_in = system.file("extdata/oisst_short.nc", package = "heatwave3"),
                       return_type = "df", file_out = "mango"), "Invalid saving option.")
  expect_error(detect3(file_in = system.file("extdata/oisst_short.nc", package = "heatwave3"),
                       return_type = "df", file_out = "mango", save_to_file = "banana"), "Invalid saving option.")
})
