#' Block average of event metrics by year and pixel
#'
#' Computes yearly summary statistics of event metrics for each pixel.
#'
#' @param event_file Path to the event NetCDF file from \code{\link{detect_event3}}.
#'
#' @return A data.frame with yearly aggregated event metrics per pixel,
#'   including count, mean/max duration, mean/max intensity, and total days.
#'
#' @export
#'
#' @examples
#' \donttest{
#' sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
#' clim_file <- tempfile(fileext = ".nc")
#' event_file <- tempfile(fileext = ".nc")
#'
#' detect3(sst_file, clim_file, event_file,
#'         climatologyPeriod = c("1982-01-01", "2011-12-31"))
#'
#' ba <- block_average3(event_file)
#' head(ba)
#' }
block_average3 <- function(event_file) {

  if (!requireNamespace("ncdf4", quietly = TRUE))
    stop("ncdf4 package is required for block_average3()", call. = FALSE)

  nc <- ncdf4::nc_open(event_file)
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  ev_lon <- ncdf4::ncvar_get(nc, "lon")
  ev_lat <- ncdf4::ncvar_get(nc, "lat")
  dur <- ncdf4::ncvar_get(nc, "duration")
  i_mean <- ncdf4::ncvar_get(nc, "intensity_mean")
  i_max <- ncdf4::ncvar_get(nc, "intensity_max")
  i_cum <- ncdf4::ncvar_get(nc, "intensity_cumulative")

  ds <- ncdf4::ncvar_get(nc, "date_start")
  ds_units <- ncdf4::ncatt_get(nc, "date_start", "units")$value
  ref_date <- as.Date(sub("days since ", "", ds_units))
  start_dates <- ref_date + ds

  # Create per-event data
  pixel_key <- paste(ev_lon, ev_lat)
  year <- as.integer(format(start_dates, "%Y"))

  df <- data.frame(
    lon = ev_lon, lat = ev_lat, pixel = pixel_key,
    year = year, duration = dur,
    intensity_mean = i_mean, intensity_max = i_max,
    intensity_cumulative = i_cum,
    stringsAsFactors = FALSE
  )

  # Aggregate by pixel and year
  result <- do.call(rbind, lapply(split(df, interaction(df$pixel, df$year, drop = TRUE)), function(sub) {
    data.frame(
      lon = sub$lon[1],
      lat = sub$lat[1],
      year = sub$year[1],
      count = nrow(sub),
      duration_mean = mean(sub$duration),
      duration_max = max(sub$duration),
      intensity_mean = mean(sub$intensity_mean),
      intensity_max_mean = mean(sub$intensity_max),
      intensity_max_max = max(sub$intensity_max),
      intensity_cumulative_mean = mean(sub$intensity_cumulative),
      total_days = sum(sub$duration),
      total_icum = sum(sub$intensity_cumulative),
      stringsAsFactors = FALSE
    )
  }))

  rownames(result) <- NULL
  result[order(result$lon, result$lat, result$year), ]
}
