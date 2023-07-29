test_that("multiplication works", {
  file_name <- "data/oisst_short.nc"
  detect_res <- detect3(file_in = file_name, return_rast = TRUE)
  expect_s4_class(detect_res, "SpatRasterDataset")
})
