#' detect3clim
#'
#' @param x Data used for detection.
#' @param time_dim The time dimension.
#'
#' @return
#' @export
#'
#' @examples
detect3clim <- function(x, time_dim){

  df_sub <- data.frame(t = as.Date(time_dim), temp = x)

  df_seas <- heatwaveR::ts2clm(df_sub, climatologyPeriod = c("1982-01-01", "2011-12-31"))

  as.matrix(df_seas[,c("temp", "seas", "thresh")])

}
