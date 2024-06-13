#' detect3clim
#'
#' Determines the climatologies within a raster file.
#' NB: Will not work with pixels with all NA values. E.g. land pixels.
#'
#' @keywords internal
#'
#' @param x Data used for detection.
#' @param time_dim The time dimension.
#' @param clim_period The climatology baseline period provided as two date values.
#' E.g. \code{c("1982-01-01", "2011-12-31")}
#'
#' @return Given the necessary file structure this will create climatologies within a raster.
#' @export
#'
detect3clim <- function(x, time_dim, clim_period, ...){

  # Calculate climatology
  df_seas <- heatwaveR::ts2clm3(data.frame(t = as.Date(time_dim), temp = x),
                                climatologyPeriod = clim_period, ...)

  # Return as matrix
  as.matrix(df_seas[,c("temp", "seas", "thresh")])
}
