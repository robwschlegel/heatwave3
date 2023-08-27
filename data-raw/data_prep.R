# This script prepares the four OISST pixels used in the package

# Libraries ---------------------------------------------------------------

# Rather call packages explicitly


# Load data ---------------------------------------------------------------

sst_q <- readr::read_csv("data-raw/sst_q.csv")


# Add to package ----------------------------------------------------------

usethis::use_data(sst_q, overwrite = TRUE)

