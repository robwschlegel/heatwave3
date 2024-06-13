#' detect3
#'
#' Function to detect events within a NetCDF file using the Hobday et al. 2016 definition.
#'
#' @param file_in A NetCDF file
#' @param clim_period The climatology baseline period provided as two date values.
#' E.g. \code{c("1982-01-01", "2011-12-31")}
#' @param min_dur The minimum duration for acceptance of detected events.
#' The default is \code{5} days.
#' @param max_gap The maximum length of gap allowed for the joining of MHWs. The
#' default is \code{2} days.
#' @param file_out File output location with name and extension of the file
#' @param return_type Default NULL will prevent the data being saved in memory.
#'                    Other options are "rast", to return a SpatRasterDataset, and "df",
#'                    to return a data.frame with the events organized by raster cell
#' @param save_to_file Default NULL will prevent the data being saved in memory.
#'                    Other options are "nc", to save a NetCDF, and "csv", to save as a csv file.
#' @param ... One may pass any arguments to this functions that would be used via
#' \code{heatwaveR::ts2clm()} or \code{heatwaveR::detect_event()}.
#'
#' @return Depending on the arguments set, this function will return the heatwaves detected in
#' the NetCDF file provide. It may also output the results as a CSV file.
#'
#' @export
#'
#' @examples
#' \donttest{
#' mhw_cube <- detect3(file_in = system.file("extdata/oisst_short.nc", package = "heatwave3"),
#'                     return_type = "df", clim_period = c("1982-01-01", "2011-12-31"))
#' }
#'
detect3 <- function(file_in, clim_period, min_dur = 5, max_gap = 2, file_out = NULL, return_type = NULL, save_to_file = NULL, ...){

  # Test that required arguments are given
  if (missing(file_in))
    stop("Please provide an input file.", call. = FALSE)
  if (missing(clim_period))
    stop("Please provide a climatology period.", call. = FALSE)
  if (length(clim_period) != 2)
    stop("Please provide BOTH start and end dates for the climatology period.")

  # Test if output types are all empty
  if (is.null(file_out) & is.null(return_type))
    stop("No output selected for the function.\nPlease enter file_out and/or return_type.", call. = FALSE)

  # Test if the output type is correct
  if (!is.null(return_type)){
    if (!(return_type %in% c("rast", "df"))) {
      stop("Invalid return_type.\nPlease enter a valid return_type (\'rast\' or \'df\')", call. = FALSE)
    }
  }

  # Test if the save format and output is correct
  if (!is.null(save_to_file)){
    if (!(save_to_file %in% c("nc", "csv"))) {
      stop("Invalid saving option.\nPlease enter a valid save_to_file (\'nc\' or \'csv\')", call. = FALSE)
    }
    if (is.null(file_out)){
      stop("Missing file destination and name.", call. = FALSE)
    }
  }

  # Load NetCDF as terra::rast
  nc_rast <- terra::rast(file_in)

  # Check for daily data
  time_dim <- terra::time(nc_rast)
  daily_dim <- unique(as.Date(terra::time(nc_rast)))
  if(length(time_dim) > length(daily_dim)){
    nc_rast_daily <- terra::tapp(nc_rast, "days", mean)
  } else {
    nc_rast_daily <- nc_rast
  }

  # Check that lon/lat range exists
  # range(nc_rast_daily@cpp$range_max, na.rm = TRUE)

  # Remove missing pixels
  nc_daily_no_NA <- nc_rast_daily[[!is.na(terra::global(nc_rast_daily, sum, na.rm = TRUE))]]

  # Create temp+seas+clim rasters
  nc_seas <- terra::app(x = nc_daily_no_NA, fun = detect3clim,
                        time_dim = terra::time(nc_daily_no_NA),
                        clim_period = clim_period, ...)

  # Add correct names
  names(nc_seas) <- c(rep(paste0("temp.", 1:terra::nlyr(nc_daily_no_NA))),
                      rep(paste0("seas.", 1:terra::nlyr(nc_daily_no_NA))),
                      rep(paste0("thresh.", 1:terra::nlyr(nc_daily_no_NA)))
  )

  # Calculate event metrics
  nc_event <- terra::app(x = nc_seas, fun = detect3event,
                         time_dim = terra::time(nc_rast_daily),
                         min_dur = min_dur, max_gap = max_gap, ...)

  ## Get the number of layers of each MHW metric
  max_layers <- terra::nlyr(nc_event)/19

  # Assign names to each MHW metric
  names(nc_event) <- c(rep(paste0("event_no.", 1:max_layers)),
                       rep(paste0("index_start.", 1:max_layers)),
                       rep(paste0("index_peak.", 1:max_layers)),
                       rep(paste0("index_end.", 1:max_layers)),
                       rep(paste0("duration.", 1:max_layers)),
                       rep(paste0("intensity_mean.", 1:max_layers)),
                       rep(paste0("intensity_max.", 1:max_layers)),
                       rep(paste0("intensity_var.", 1:max_layers)),
                       rep(paste0("intensity_cumulative.", 1:max_layers)),
                       rep(paste0("intensity_mean_relThresh.", 1:max_layers)),
                       rep(paste0("intensity_max_relThresh.", 1:max_layers)),
                       rep(paste0("intensity_var_relThresh.", 1:max_layers)),
                       rep(paste0("intensity_cumulative_relThresh.", 1:max_layers)),
                       rep(paste0("intensity_mean_abs.", 1:max_layers)),
                       rep(paste0("intensity_max_abs.", 1:max_layers)),
                       rep(paste0("intensity_var_abs.", 1:max_layers)),
                       rep(paste0("intensity_cumulative_abs.", 1:max_layers)),
                       rep(paste0("rate_onset.", 1:max_layers)),
                       rep(paste0("rate_decline.", 1:max_layers))
  )

  # Remove layers with no data
  # nc_no_NA <- nc_event[[1]]
  nc_no_NA <- nc_event[[!is.na(terra::global(nc_event, sum, na.rm = TRUE))]]

  # Create sds object
  nc_sds <- terra::sds(
    nc_no_NA[[grepl("event_no.", names(nc_no_NA))]],
    nc_no_NA[[grepl("index_start.", names(nc_no_NA))]],
    nc_no_NA[[grepl("index_peak.", names(nc_no_NA))]],
    nc_no_NA[[grepl("index_end.", names(nc_no_NA))]],
    nc_no_NA[[grepl("duration.", names(nc_no_NA))]],
    nc_no_NA[[grepl("intensity_mean.", names(nc_no_NA), fixed = T)]],
    nc_no_NA[[grepl("intensity_max.", names(nc_no_NA), fixed = T)]],
    nc_no_NA[[grepl("intensity_var.", names(nc_no_NA), fixed = T)]],
    nc_no_NA[[grepl("intensity_cumulative.", names(nc_no_NA), fixed = T)]],
    nc_no_NA[[grepl("intensity_mean_relThresh.", names(nc_no_NA))]],
    nc_no_NA[[grepl("intensity_max_relThresh.", names(nc_no_NA))]],
    nc_no_NA[[grepl("intensity_var_relThresh.", names(nc_no_NA))]],
    nc_no_NA[[grepl("intensity_cumulative_relThresh.", names(nc_no_NA))]],
    nc_no_NA[[grepl("intensity_mean_abs.", names(nc_no_NA))]],
    nc_no_NA[[grepl("intensity_max_abs.", names(nc_no_NA))]],
    nc_no_NA[[grepl("intensity_var_abs.", names(nc_no_NA))]],
    nc_no_NA[[grepl("intensity_cumulative_abs.", names(nc_no_NA))]],
    nc_no_NA[[grepl("rate_onset.", names(nc_no_NA))]],
    nc_no_NA[[grepl("rate_decline.", names(nc_no_NA))]]
  )

  names(nc_sds) <- c("event_no", "index_start", "index_peak" ,"index_end",
                     "duration", "intensity_mean", "intensity_max" , "intensity_var",
                     "intensity_cumulative", "intensity_mean_relThresh", "intensity_max_relThresh",
                     "intensity_var_relThresh","intensity_cumulative_relThresh", "intensity_mean_abs",
                     "intensity_max_abs" , "intensity_var_abs", "intensity_cumulative_abs" ,"rate_onset", "rate_decline")

  # Save as NetCDF
  if(!is.null(save_to_file)) {
    if (save_to_file == "nc") {
      terra::longnames(nc_sds) <- c("A sequential number indicating the ID and order of the events",
                                    paste0("Start date of event [date], as days since ", min(terra::time(nc_rast_daily))),
                                    paste0("Date of event peak [date], as days since ", min(terra::time(nc_rast_daily))),
                                    paste0("End date of event [date], as days since ", min(terra::time(nc_rast_daily))),
                                    "Duration of event [days]",
                                    "Mean intensity [deg. C]",
                                    "Maximum (peak) intensity [deg. C]",
                                    "Intensity variability (standard deviation) [deg. C]",
                                    "Cumulative intensity [deg. C x days]",
                                    "Mean intensity [deg. C] relative to the threshold (e.g., 90th percentile)",
                                    "Maximum (peak) intensity [deg. C] relative to the threshold (e.g., 90th percentile)",
                                    "Intensity variability (standard deviation) [deg. C] relative to the threshold (e.g., 90th percentile)",
                                    "Cumulative intensity [deg. C x days] relative to the threshold (e.g., 90th percentile)",
                                    "Mean absolute intensity [deg. C]",
                                    "Maximum (peak) absolute intensity [deg. C]",
                                    "Absolute intensity variability (standard deviation) [deg. C]",
                                    "Absolute cumulative intensity [deg. C x days]",
                                    "Onset rate of event [deg. C / day]",
                                    "Decline rate of event [deg. C / day]")
      terra::writeCDF(nc_sds, file_out, overwrite = TRUE)
    }
  }

  # Save as csv
  if(!is.null(save_to_file)) {
    if (save_to_file == "csv") {
      nc_csv <- rast_to_df(x = nc_no_NA, time_dim = min(terra::time(nc_rast_daily)))
      utils::write.csv(nc_csv, file = file_out, row.names = F)
    }
  }

  # Output results to R environment
  if(!is.null(return_type)){

    if(return_type == "rast"){
      #comment(nc_sds) <- paste0("Dates on indexes based on days since ",  min(terra::time(nc_rast_daily)))
      return(nc_sds)
    }

    # Transform the SpatRasterDataset into a data.frame
    if(return_type == "df"){
      nc_csv <- rast_to_df(x = nc_no_NA, time_dim = min(terra::time(nc_rast_daily)))
      return(nc_csv)
    }

  } else {
    print("Finished. Good job team :)")
  }
}
