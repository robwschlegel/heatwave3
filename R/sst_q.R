#' NOAA Optimally Interpolated (OI) v2.1 daily 1/4 degree SST
#'
#' A dataset containing the sea surface temperature (°C)
#' from 1982-01-01 to 2022-12-31.
#'
#' longitude: 112.625 & 112.875
#' latitude: -29.375 &  -29.125
#'
#' @format A CSV file with 14975 time steps for 4 pixels:
#' \describe{
#'   \item{temp}{sea surface temperature (SST) (°C)}
#'   \item{time}{days since 1970-01-01}
#'   \item{longitude}{degrees easting}
#'   \item{latitude}{degrees northing}
#'   ...
#' }
#' @source \url{https://www.ncei.noaa.gov/products/optimum-interpolation-sst}
"sst_q"
