#' Plot a per-pixel MHW time series from NetCDF output
#'
#' Extracts a single pixel's time series from the SST and climatology NetCDF
#' files and produces an event_line-style plot with flame polygons. The SST
#' and climatology data are read via the C++ backend; no external NetCDF
#' packages are required.
#'
#' @param sst_file Path to the SST NetCDF file (or directory of daily files).
#' @param clim_file Path to the climatology NetCDF from \code{\link{ts2clm3}}.
#' @param lon Longitude of the pixel to plot.
#' @param lat Latitude of the pixel to plot.
#' @param depth Optional depth (metres) to plot, for a depth-resolved
#'   \code{clim_file} (from \code{ts2clm3(depth_range = ...)}). Matched to the
#'   nearest depth level actually present in \code{clim_file}; the matched
#'   value is used to read \code{sst_file} at the same level and appears in
#'   the plot title. Required when \code{clim_file} is depth-resolved;
#'   must be \code{NULL} (the default) for an ordinary 3D \code{clim_file}.
#' @param var_name SST variable name. If \code{NULL}, auto-detected.
#' @param start_date,end_date Optional date range for the plot window
#'   (character, for example \code{"2010-01-01"}). If both are \code{NULL}, the
#'   plot is centred on the most intense event (see \code{spread}).
#' @param spread Number of days before and after the peak event to show.
#'   Default \code{150}. Only used when \code{start_date}/\code{end_date}
#'   are not set.
#' @param metric Event metric to use for selecting the peak event when
#'   \code{event_file} is supplied. Default \code{"intensity_cumulative"}.
#' @param event_file Optional path to event NetCDF for centring the window
#'   on the most intense event.
#' @param coldSpells Logical. Render cold-spell (blue) or heatwave (red)
#'   flames? Default \code{FALSE}.
#'
#' @return A ggplot object that can be further customised with the
#'   standard \code{+} operator.
#'
#' @export
#'
#' @examples
#' \donttest{
#' sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
#' stem <- file.path(tempdir(), "demo")
#'
#' ts2clm3(sst_file, name = stem,
#'         climatologyPeriod = c("1982-01-01", "2011-12-31"))
#'
#' event_line3(sst_file, paste0(stem, "_clim.nc"),
#'             lon = 26.525, lat = -34.125,
#'             start_date = "2010-01-01", end_date = "2012-12-31")
#' }
event_line3 <- function(sst_file, clim_file, lon, lat,
                        depth = NULL,
                        var_name = NULL,
                        start_date = NULL, end_date = NULL,
                        spread = 150,
                        metric = "intensity_cumulative",
                        event_file = NULL,
                        coldSpells = FALSE) {

  eps <- 0.01
  vn <- if (is.null(var_name)) "" else var_name

  cd <- hw3_read_clim_nc(clim_file)
  lon_idx <- which.min(abs(cd$lon - lon))
  lat_idx <- which.min(abs(cd$lat - lat))
  ndepth <- max(1L, cd$ndepth)

  if (ndepth > 1L && is.null(depth)) {
    stop("clim_file is depth-resolved (ndepth = ", ndepth, "); specify 'depth' ",
         "(metres) to pick a level.", call. = FALSE)
  }
  if (ndepth <= 1L && !is.null(depth)) {
    stop("clim_file has no depth dimension; omit 'depth'.", call. = FALSE)
  }

  # Resolve to the exact depth value stored in clim_file, so the SST read's
  # depth_range = c(resolved, resolved) is guaranteed to match (avoids the
  # user's rough target value, e.g. 50, missing the file's exact 47.37).
  idepth <- 0L
  resolved_depth <- NULL
  depth_range <- NULL
  if (!is.null(depth)) {
    idepth <- which.min(abs(cd$depth - depth)) - 1L
    resolved_depth <- cd$depth[idepth + 1L]
    depth_range <- c(resolved_depth, resolved_depth)
  }
  px <- ((lon_idx - 1L) * cd$nlat + (lat_idx - 1L)) * ndepth + idepth

  # Handle directory or multi-file input
  if (length(sst_file) == 1 && dir.exists(sst_file)) {
    files <- sort(list.files(sst_file, pattern = "\\.(nc|nc4)$",
                             full.names = TRUE))
    if (length(files) == 0)
      stop("No .nc/.nc4 files in ", sst_file, call. = FALSE)
    gd <- hw3_read_sst_multi(files, vn,
                              lon_range = c(lon - eps, lon + eps),
                              lat_range = c(lat - eps, lat + eps),
                              depth_range = depth_range)
  } else {
    gd <- hw3_read_sst(sst_file, vn,
                       lon_range = c(lon - eps, lon + eps),
                       lat_range = c(lat - eps, lat + eps),
                       depth_range = depth_range)
  }
  if (gd$nlon * gd$nlat == 0)
    stop("No data found at lon = ", lon, ", lat = ", lat, call. = FALSE)

  temp <- gd$sst[seq_len(gd$ntime)]
  dates <- as.Date("1970-01-01") + (gd$time_days - 2440588L)

  doy_vals <- vapply(gd$time_days, hw3_jd_to_doy, integer(1))
  seas <- cd$seas[px * 366 + doy_vals]
  thresh <- cd$thresh[px * 366 + doy_vals]

  # Convert kelvin to Celsius if values are clearly in kelvin
  if (stats::median(temp, na.rm = TRUE) > 200) {
    temp <- temp - 273.15
    seas <- seas - 273.15
    thresh <- thresh - 273.15
  }

  df <- data.frame(t = dates, temp = temp, seas = seas, thresh = thresh)
  df <- df[stats::complete.cases(df), ]

  # Determine plot window
  if (!is.null(start_date) && !is.null(end_date)) {
    df <- df[df$t >= as.Date(start_date) & df$t <= as.Date(end_date), ]
  } else if (!is.null(event_file)) {
    ev <- hw3_read_event_nc(event_file)
    ref_date <- as.Date("1970-01-01") + (ev$ref_date_jd - 2440588L)
    peak_dates <- ref_date + ev$date_peak
    ev_metric <- ev[[metric]]

    px_mask <- abs(ev$lon - lon) < eps & abs(ev$lat - lat) < eps
    if (!is.null(resolved_depth) && !is.null(ev$depth)) {
      px_mask <- px_mask & abs(ev$depth - resolved_depth) < 1e-3
    }
    if (any(px_mask)) {
      best <- which(px_mask)[which.max(ev_metric[px_mask])]
      center <- peak_dates[best]
      df <- df[df$t >= (center - spread) & df$t <= (center + spread), ]
    }
  }

  if (nrow(df) == 0)
    stop("No data in the specified date range.", call. = FALSE)

  flame_fill <- if (coldSpells) "steelblue3" else "salmon"

  ggplot2::ggplot(df, ggplot2::aes(x = .data$t)) +
    geom_flame3(ggplot2::aes(y = .data$temp, y2 = .data$thresh),
                fill = flame_fill, reverse = coldSpells) +
    ggplot2::geom_line(ggplot2::aes(y = .data$temp), linewidth = 0.3) +
    ggplot2::geom_line(ggplot2::aes(y = .data$seas),
                       colour = "grey40", linetype = "dashed",
                       linewidth = 0.3) +
    ggplot2::geom_line(ggplot2::aes(y = .data$thresh),
                       colour = "darkgreen", linewidth = 0.3) +
    ggplot2::labs(
      x = NULL, y = expression("Temperature (" * degree * "C)"),
      title = if (is.null(resolved_depth)) {
        sprintf("lon = %.3f, lat = %.3f", lon, lat)
      } else {
        sprintf("lon = %.3f, lat = %.3f, depth = %.1f m", lon, lat, resolved_depth)
      }
    ) +
    ggplot2::theme_minimal()
}
