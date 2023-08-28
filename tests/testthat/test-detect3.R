# Tests for the primary function: detect3()

test_that("The function runs as expected", {
  # test_nc <- system.file("extdata/oisst_short.nc", package = "heatwave3")
  res <- detect3(file_in = system.file("extdata/oisst_short.nc", package = "heatwave3"),
                 return_type = "df")
  expect_s3_class(res, "data.frame")
  expect_equal(ncol(res), 22)
  expect_equal(nrow(res), 340)
})
