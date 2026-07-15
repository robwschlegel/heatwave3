#' Spatial map of event metrics
#'
#' Creates a map showing the spatial distribution of a chosen event metric
#' from the event NetCDF output. The per-pixel aggregation is performed in
#' C++ for efficiency; the R layer handles only the ggplot2 rendering.
#'
#' @param event_file Path to the event NetCDF file from \code{\link{detect_event3}}.
#' @param metric Character. The event metric to map. Options include
#'   \code{"intensity_max"} (default), \code{"intensity_mean"},
#'   \code{"intensity_cumulative"}, \code{"duration"}, \code{"rate_onset"},
#'   \code{"rate_decline"}, and all \code{relThresh}/\code{abs} variants.
#' @param summary Character. How to aggregate across events per pixel.
#'   One of \code{"mean"} (default), \code{"max"}, \code{"min"},
#'   \code{"sum"}, or \code{"count"}.
#' @param coastline Logical. Add a coastline layer? Requires the
#'   \code{rnaturalearth} and \code{sf} packages. Default \code{TRUE}.
#' @param ... Additional arguments passed to
#'   \code{ggplot2::scale_fill_viridis_c}.
#'
#' @return A ggplot object. The underlying data is accessible via
#'   \code{ggplot2::layer_data()} or by calling
#'   \code{hw3_read_metric_summary()} directly.
#'
#' @export
#'
#' @examples
#' \donttest{
#' sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
#' stem <- file.path(tempdir(), "demo")
#'
#' detect3(sst_file, name = stem,
#'         climatologyPeriod = c("1982-01-01", "2011-12-31"))
#'
#' plot_metric3(paste0(stem, "_events.nc"), metric = "intensity_max", summary = "mean")
#' }
plot_metric3 <- function(event_file,
                         metric = "intensity_max",
                         summary = "mean",
                         coastline = TRUE,
                         ...) {

  agg <- hw3_read_metric_summary(event_file, metric = metric, summary = summary)

  label <- paste0(metric, " (", summary, ")")

  p <- ggplot2::ggplot(agg, ggplot2::aes(
    x = .data$lon, y = .data$lat, fill = .data$value
  )) +
    ggplot2::geom_raster() +
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
