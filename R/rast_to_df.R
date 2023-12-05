#' rast_to_df
#'
#' Function for converting a raster to a dataframe.
#'
#' @keywords internal
#'
#' @param x Rast stack resulted from detect3event and re-organized until 'nc_no_NA' (see detect3 function)
#' @param time_dim The time dimension. It should be a single value to use as the origin for the date indexes
#'
#' @return This function will convert a raster object to a dataframe.
#' @export
#'
rast_to_df <- function(x, time_dim){

  nc_csv <- terra::as.data.frame(x, xy = T, cells = T)
  nc_csv <- nc_csv[,!colnames(nc_csv) %in% colnames(nc_csv)[grepl("event_no.", colnames(nc_csv))]]
  nc_csv <- tidyr::pivot_longer(data = nc_csv, cols = !c("cell", "x", "y"),
                                names_to = c(".value", "event_no"), names_sep = "[.]")
  nc_csv <- nc_csv[stats::complete.cases(nc_csv),]

  nc_csv$index_start <- as.Date(nc_csv$index_start, origin = time_dim)
  nc_csv$index_peak <- as.Date(nc_csv$index_peak, origin = time_dim)
  nc_csv$index_end <- as.Date(nc_csv$index_end, origin = time_dim)

  return(nc_csv)
}
