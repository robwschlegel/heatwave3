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
  nc_event <- terra::app(nc_seas, detect3event, time_dim = terra::time(nc_rast))

  # Remove NA layers
  # nc_no_NA <- nc_event[[!terra::global(is.na(nc_event), sum) == 4]]

  # Create sds object
  nc_sds <- terra::sds(nc_event[[1:200]], nc_event[[201:400]])
  names(nc_sds) <- c("int_mean", "int_max")


  # Save as desired
  if(!is.null(file_out)) terra::writeCDF(nc_sds, file_out, overwrite = TRUE)

  # Output results
  if(return_rast){
    return(nc_sds)
  } else {
    print("Finished. Good job team :)")
  }
}
