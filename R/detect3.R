#' detect3
#'
#' Function to detect events within a NetCDF file using the Hobday et al. 2016 definition.
#'
#' @param file_in A NetCDF file
#' @param file_out File output location
#' @param return_rast Default FALSE will prevent the data being saved in memory.
#'                    Other options are "rast", to return a SpatRasterDataset, and "df", to return a data.frame with the events organized by raster cell
#'
#' @return
#' @export
#'
#' @examples

detect3 <- function(file_in, file_out = NULL, return_type = NULL){


  # Test if output types are all empty
  if (is.null(file_out) & is.null(return_type)) {
    stop("No output selected for the function.\nPlease enter file_out and/or return_type.", call. = FALSE)
  }

  # Test if the output format if correct
  if (!(return_type %in% c("rast", "df"))) {
    stop(shQuote("Invalid return_type.\nPlease enter a valid return_type (\'rast\' or \'df\')"), call. = FALSE)
  }

  # file_in <- "data/oisst_short.nc"
  # rm(file_in, file_out, return_rast, y, nc_seas); gc()

  # Load NetCDF as terra::rast
  nc_rast <- terra::rast(file_in)

  # Create temp+seas+clim rasters
  nc_seas <- terra::app(x = nc_rast, fun = detect3clim, time_dim = terra::time(nc_rast))

  # Add correct names
  names(nc_seas) <- c(rep(paste0("temp.", 1:terra::nlyr(nc_rast))),
                      rep(paste0("seas.", 1:terra::nlyr(nc_rast))),
                      rep(paste0("thresh.", 1:terra::nlyr(nc_rast)))
  )

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

  # Output results to R environment
  if(!is.null(return_type)){

    if(return_type == "rast"){
      #comment(nc_sds) <- paste0("Dates on indexes based on days since ",  min(terra::time(nc_rast)))
      return(nc_sds)
    }

    # Transform the SpatRasterDataset into a data.frame
    if(return_type == "df"){
      nc_csv <- terra::as.data.frame(nc_no_NA, xy = T, cells = T)
      nc_csv <- nc_csv[,!colnames(nc_csv) %in% colnames(nc_csv)[grepl("event_no.", colnames(nc_csv))]]
      nc_csv <- tidyr::pivot_longer(data = nc_csv, cols = !c("cell", "x", "y"), names_to = c(".value", "event_no"), names_sep = "[.]")
      nc_csv <- nc_csv[complete.cases(nc_csv),]

      nc_csv$index_start <- as.Date(nc_csv$index_start, origin = min(terra::time(nc_rast)))
      nc_csv$index_peak <- as.Date(nc_csv$index_peak, origin = min(terra::time(nc_rast)))
      nc_csv$index_end <- as.Date(nc_csv$index_end, origin = min(terra::time(nc_rast)))

      #write.csv2(nc_csv, file = file_out, row.names = F)
      return(nc_csv)
    }

  } else {
    print("Finished. Good job team :)")
  }
}
