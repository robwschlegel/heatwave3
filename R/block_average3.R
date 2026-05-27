#' Block average of event metrics by year and pixel
#'
#' Computes yearly summary statistics of event metrics for each pixel.
#'
#' @param event_file Path to the event NetCDF file from \code{\link{detect_event3}}.
#'
#' @return A data.frame with yearly aggregated event metrics per pixel,
#'   including count, mean/max duration, mean/max intensity, and total days.
#'
#' @export
#'
#' @examples
#' \donttest{
#' sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
#' clim_file <- tempfile(fileext = ".nc")
#' event_file <- tempfile(fileext = ".nc")
#'
#' detect3(sst_file, clim_file, event_file,
#'         climatologyPeriod = c("1982-01-01", "2011-12-31"))
#'
#' ba <- block_average3(event_file)
#' head(ba)
#' }
block_average3 <- function(event_file) {
  hw3_block_average(event_file)
}
