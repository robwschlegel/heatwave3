#' Plot a per-pixel MHW time series from NetCDF output
#'
#' Extracts a single pixel's time series from the SST and climatology NetCDF
#' files and produces an event_line-style plot with flame polygons.
#'
#' @param sst_file Path to the SST NetCDF file.
#' @param clim_file Path to the climatology NetCDF from \code{\link{ts2clm3}}.
#' @param lon Longitude of the pixel to plot.
#' @param lat Latitude of the pixel to plot.
#' @param var_name SST variable name. If \code{NULL}, auto-detected.
#' @param start_date Optional start date for the plot window.
#' @param end_date Optional end date for the plot window.
#' @param spread Number of days around the peak event to show. Default \code{150}.
#'   Only used when \code{start_date} and \code{end_date} are not set.
#' @param metric Event metric to use for selecting the peak event.
#'   Default \code{"intensity_cumulative"}.
#' @param event_file Optional path to event NetCDF for event highlighting.
#' @param coldSpells Logical. Plot cold spells? Default \code{FALSE}.
#'
#' @return A ggplot object.
#'
#' @export
#'
#' @examples
#' \donttest{
#' sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
#' clim_file <- tempfile(fileext = ".nc")
#'
#' ts2clm3(sst_file, clim_file,
#'         climatologyPeriod = c("1982-01-01", "2011-12-31"))
#'
#' event_line3(sst_file, clim_file,
#'             lon = 26.525, lat = -34.125,
#'             start_date = "2010-01-01", end_date = "2012-12-31")
#' }
event_line3 <- function(sst_file, clim_file, lon, lat,
                        var_name = NULL,
                        start_date = NULL, end_date = NULL,
                        spread = 150,
                        metric = "intensity_cumulative",
                        event_file = NULL,
                        coldSpells = FALSE) {

  # Read pixel data
  eps <- 0.01
  vn <- if (is.null(var_name)) "" else var_name

  gd <- hw3_read_sst(sst_file, vn,
                     lon_range = c(lon - eps, lon + eps),
                     lat_range = c(lat - eps, lat + eps))

  if (gd$nlon * gd$nlat == 0) stop("No data found at the specified lon/lat.", call. = FALSE)

  temp <- gd$sst[1:gd$ntime]
  dates <- as.Date("1970-01-01") + (gd$time_days - 2440588L)

  # Read climatology for this pixel
  cd <- hw3_read_clim_nc(clim_file)
  # Find matching pixel
  lon_idx <- which.min(abs(cd$lon - lon))
  lat_idx <- which.min(abs(cd$lat - lat))
  px <- (lon_idx - 1L) * cd$nlat + (lat_idx - 1L)

  # Map DOY to get daily seas/thresh
  doy_vals <- vapply(gd$time_days, hw3_jd_to_doy, integer(1))
  seas <- cd$seas[px * 366 + doy_vals]
  thresh <- cd$thresh[px * 366 + doy_vals]

  df <- data.frame(t = dates, temp = temp, seas = seas, thresh = thresh)
  df <- df[complete.cases(df), ]

  # Determine plot window
  if (!is.null(start_date) && !is.null(end_date)) {
    df <- df[df$t >= as.Date(start_date) & df$t <= as.Date(end_date), ]
  } else if (!is.null(event_file)) {
    # Use the most intense event
    nc <- ncdf4::nc_open(event_file)
    ev_lon <- ncdf4::ncvar_get(nc, "lon")
    ev_lat <- ncdf4::ncvar_get(nc, "lat")
    ev_ic <- ncdf4::ncvar_get(nc, metric)
    dp <- ncdf4::ncvar_get(nc, "date_peak")
    dp_units <- ncdf4::ncatt_get(nc, "date_peak", "units")$value
    ref_date <- as.Date(sub("days since ", "", dp_units))
    peak_dates <- ref_date + dp
    ncdf4::nc_close(nc)

    px_mask <- abs(ev_lon - lon) < eps & abs(ev_lat - lat) < eps
    if (any(px_mask)) {
      best <- which(px_mask)[which.max(ev_ic[px_mask])]
      center <- peak_dates[best]
      df <- df[df$t >= (center - spread) & df$t <= (center + spread), ]
    }
  }

  if (nrow(df) == 0) stop("No data in the specified date range.", call. = FALSE)

  # Build plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$t)) +
    geom_flame3(ggplot2::aes(y = .data$temp, y2 = .data$thresh),
                fill = if (coldSpells) "steelblue3" else "salmon") +
    ggplot2::geom_line(ggplot2::aes(y = .data$temp), linewidth = 0.3) +
    ggplot2::geom_line(ggplot2::aes(y = .data$seas), colour = "grey40",
                       linetype = "dashed", linewidth = 0.3) +
    ggplot2::geom_line(ggplot2::aes(y = .data$thresh), colour = "darkgreen",
                       linewidth = 0.3) +
    ggplot2::labs(x = "Date", y = "Temperature",
                  title = sprintf("lon = %.3f, lat = %.3f", lon, lat)) +
    ggplot2::theme_minimal()

  p
}
