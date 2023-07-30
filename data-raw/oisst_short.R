## code to prepare `oisst_short` dataset goes here

oisst_short <- terra::rast("data-raw/oisst_short.nc")

usethis::use_data(oisst_short, overwrite = TRUE)
