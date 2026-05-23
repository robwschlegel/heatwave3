#' Spatial map of event metrics
#'
#' Creates a map showing the spatial distribution of a chosen event metric
#' from the event NetCDF output.
#'
#' @param event_file Path to the event NetCDF file from \code{\link{detect_event3}}.
#' @param metric Character. The event metric to map. Options include
#'   \code{"intensity_max"}, \code{"intensity_mean"}, \code{"duration"},
#'   \code{"intensity_cumulative"}, etc. Default \code{"intensity_max"}.
#' @param summary Character. How to aggregate across events per pixel.
#'   One of \code{"mean"}, \code{"max"}, \code{"count"}, or \code{"sum"}.
#'   Default \code{"mean"}.
#' @param coastline Logical. Add a coastline layer? Requires the
#'   \code{rnaturalearth} package. Default \code{TRUE}.
#' @param ... Additional arguments passed to \code{ggplot2::scale_fill_viridis_c}.
#'
#' @return A ggplot object.
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
#' plot_metric3(event_file, metric = "intensity_max", summary = "mean")
#' }
plot_metric3 <- function(event_file,
                         metric = "intensity_max",
                         summary = "mean",
                         coastline = TRUE,
                         ...) {

  if (!requireNamespace("ncdf4", quietly = TRUE))
    stop("ncdf4 package is required for plot_metric3()", call. = FALSE)

  nc <- ncdf4::nc_open(event_file)
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  ev_lon <- ncdf4::ncvar_get(nc, "lon")
  ev_lat <- ncdf4::ncvar_get(nc, "lat")
  vals <- ncdf4::ncvar_get(nc, metric)

  df <- data.frame(lon = ev_lon, lat = ev_lat, value = vals)

  # Aggregate per pixel
  agg_fun <- switch(summary,
    "mean" = mean,
    "max"  = max,
    "min"  = min,
    "sum"  = sum,
    "count" = length,
    stop("Unknown summary function: ", summary, call. = FALSE)
  )

  agg <- stats::aggregate(value ~ lon + lat, data = df, FUN = agg_fun)
  names(agg)[3] <- "value"

  label <- paste0(metric, " (", summary, ")")

  p <- ggplot2::ggplot(agg, ggplot2::aes(x = .data$lon, y = .data$lat, fill = .data$value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(name = label, ...) +
    ggplot2::coord_quickmap() +
    ggplot2::labs(x = "Longitude", y = "Latitude") +
    ggplot2::theme_minimal()

  if (coastline) {
    if (requireNamespace("rnaturalearth", quietly = TRUE) &&
        requireNamespace("sf", quietly = TRUE)) {
      coast <- rnaturalearth::ne_coastline(scale = "medium", returnclass = "sf")
      p <- p + ggplot2::geom_sf(data = coast, inherit.aes = FALSE,
                                colour = "grey30", linewidth = 0.3)
    }
  }

  p
}
