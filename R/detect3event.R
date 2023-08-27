#' detect3event
#'
#' @param x Data used for detection.
#' @param time_dim The time dimension.
#' @param dur The minimum duration (days) for an event to be detected.
#' @param gap The maximum gap (days) across which two events will be considered one.
#'
#' @return Given the necessary file structure this will detect events within a raster.
#' @export
#'
#' @examples
detect3event <- function(x, time_dim, dur = 5, gap = 2){

  # Get the maximum possible number of layers we need
  max_layers <- round(length(time_dim)/(dur+gap+1), digits = 0)

  # Get columns
  temp_vec <- x[grepl(pattern = "temp", x = names(x))]
  seas_vec <- x[grepl(pattern = "seas", x = names(x))]
  thresh_vec <- x[grepl(pattern = "thresh", x = names(x))]

  # Create data.frame
  df_sub <- data.frame(t = as.Date(time_dim), temp = temp_vec, seas = seas_vec, thresh = thresh_vec)

  # Detect events
  df_event <- heatwaveR::detect_event(df_sub, minDuration = dur, maxGap = gap)$event

  # keep only the date indices and not date characters
  as.matrix(df_event[1:max_layers, c(1:5, 9:22)])
}
