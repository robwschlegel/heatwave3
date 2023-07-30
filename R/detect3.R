#' detect3
#'
#' Function to detect events within a NetCDF file using the Hobday et al. 2016 definition.
#'
#' @param file_in A NetCDF file
#' @param file_out File output location
#' @param return_rast Default FALSE will prevent the data being saved in memory.
#'
#' @return
#' @export
#'
#' @examples
detect3 <- function(file_in, file_out = NULL, return_rast = FALSE){

  # file_in <- "data/oisst_short.nc"
  # rm(file_in, file_out, return_rast, y, nc_seas); gc()

  # Load NetCDF as terra::rast
  nc_rast <- terra::rast(file_in)

  # Create temp+seas+clim rasters
  nc_seas <- terra::app(x = nc_rast, fun = detect3clim, time_dim = terra::time(nc_rast))

  # Add correct names
  names(nc_seas) <- c(rep(paste0("temp.", 1:length(names(nc_rast)))),
                      rep(paste0("seas.", 1:length(names(nc_rast)))),
                      rep(paste0("thresh.", 1:length(names(nc_rast)))))

  # Calculate event metrics
  nc_event <- terra::app(x = nc_seas, fun = detect3event,
                         time_dim = terra::time(nc_rast))#, ndays=terra::nlyr(nc_rast))

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
  nc_no_NA <- nc_event[[!is.na(terra::global(nc_event, sum, na.rm=TRUE))]]

  # Create sds object
  nc_sds <- terra::sds(
    nc_no_NA[[grepl("event_no.", names(nc_no_NA))]],
    nc_no_NA[[grepl("index_start.", names(nc_no_NA))]],
    nc_no_NA[[grepl("index_peak.", names(nc_no_NA))]],
    nc_no_NA[[grepl("index_end.", names(nc_no_NA))]],
    nc_no_NA[[grepl("duration.", names(nc_no_NA))]],
    nc_no_NA[[grepl("intensity_mean.", names(nc_no_NA), fixed=T)]],
    nc_no_NA[[grepl("intensity_max.", names(nc_no_NA), fixed=T)]],
    nc_no_NA[[grepl("intensity_var.", names(nc_no_NA), fixed=T)]],
    nc_no_NA[[grepl("intensity_cumulative.", names(nc_no_NA), fixed=T)]],
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

  # Save as desired
  if(!is.null(file_out)) terra::writeCDF(nc_sds, file_out, overwrite = TRUE)

  # Output results
  if(return_rast){
    return(nc_sds)
  } else {
    print("Finished. Good job team :)")
  }
}
