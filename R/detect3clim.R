#' detect3clim
#'
#' @param x Data used for detection.
#' @param time_dim The time dimension.
#'
#' @return Given the necessary file structure this will create climatologies within a raster.
#' @export
#'
#' @examples
detect3clim <- function(x, time_dim){

  # Create data.frame
  df_sub <- data.frame(t = as.Date(time_dim), temp = x)

  # Calculate climatology
  df_seas <- heatwaveR::ts2clm(df_sub, climatologyPeriod = c("1982-01-01", "2011-12-31"))

  # Return as matrix
  as.matrix(df_seas[,c("temp", "seas", "thresh")])
}
