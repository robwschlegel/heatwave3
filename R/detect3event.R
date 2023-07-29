#' detect3event
#'
#' @param x
#' @param y
#'
#' @return
#' @export
#'
#' @examples
detect3event <- function(x, time_dim){

  # Get columns
  # temp_vec <- x[1:n]
  temp_vec <- x[grepl(pattern = "temp", x = names(x))]
  seas_vec <- x[grepl(pattern = "seas", x = names(x))]
  thresh_vec <- x[grepl(pattern = "thresh", x = names(x))]
  # seas_vec <- x[n+1:n*2]
  # thresh_vec <- x[n*2+1:n*3]

  # Create data.frame
  df_sub <- data.frame(t = as.Date(time_dim), temp = temp_vec, seas = seas_vec, thresh = thresh_vec)

  df_event <- heatwaveR::detect_event(df_sub)$event
  as.matrix(df_event[1:200, 9:10])
}
