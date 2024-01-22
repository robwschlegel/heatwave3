#' detect3clim
#'
#' Determines the climatologies within a raster file.
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

  # Create data.frame
  df_sub <- data.frame(t = as.Date(time_dim), temp = x)

  # Check for missing pixels
  df_check <- df_sub[!is.na(df_sub$temp),]
  if (nrow(df_check) < nrow(df_sub)) {
    df_seas <- data.frame(t = time, temp = NA, seas = NA, thresh = NA)
  } else {
    # Calculate climatology
    df_seas <- heatwaveR::ts2clm3(df_sub, climatologyPeriod = clim_period, ...)
  }

  # Return as matrix
  as.matrix(df_seas[,c("temp", "seas", "thresh")])
}
