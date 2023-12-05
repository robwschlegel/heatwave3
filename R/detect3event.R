#' detect3event
#'
#' Detects the marine heatwaves within a raster file based on the climatologies.
#'
#' @keywords internal
#'
#' @param x Data used for detection.
#' @param time_dim The time dimension.
#' @param min_dur The minimum duration for acceptance of detected events.
#' The default is \code{5} days.
#' @param max_gap The maximum length of gap allowed for the joining of MHWs. The
#' default is \code{2} days.
#'
#' @return Given the necessary file structure this will detect events within a raster.
#' @export
#'
detect3event <- function(x, time_dim, min_dur, max_gap, ...){

  # Get the maximum possible number of layers we need
  max_layers <- round(length(time_dim)/(min_dur+max_gap+1), digits = 0)

  # Get columns
  temp_vec <- x[grepl(pattern = "temp", x = names(x))]
  seas_vec <- x[grepl(pattern = "seas", x = names(x))]
  thresh_vec <- x[grepl(pattern = "thresh", x = names(x))]

  # Create data.frame
  df_sub <- data.frame(t = as.Date(time_dim), temp = temp_vec, seas = seas_vec, thresh = thresh_vec)

  # Detect events
  df_event <- heatwaveR::detect_event3(df_sub, minDuration = min_dur, maxGap = max_gap, ...)$event

  # keep only the date indices and not date characters
  as.matrix(df_event[1:max_layers, c(1:5, 9:22)])
}
