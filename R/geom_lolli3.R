#' Lollipop geom for event metric visualisation
#'
#' A ggplot2 geom that draws lollipop (stem + point) plots for event metrics.
#'
#' @param mapping Aesthetic mappings. Requires \code{x} and \code{y}.
#' @param data The data to display.
#' @param stat Statistical transformation. Default \code{"identity"}.
#' @param position Position adjustment. Default \code{"identity"}.
#' @param ... Additional arguments.
#' @param na.rm Remove NAs? Default \code{FALSE}.
#' @param show.legend Show legend? Default \code{NA}.
#' @param inherit.aes Inherit aesthetics? Default \code{TRUE}.
#'
#' @return A ggplot2 layer.
#'
#' @export
#'
#' @examples
#' library(ggplot2)
#' df <- data.frame(x = as.Date("2020-01-01") + seq(0, 300, by = 30),
#'                  y = c(1.2, 2.1, 0.8, 3.5, 1.9, 2.7, 0.5, 1.1, 2.3, 1.6, 0.9))
#' ggplot(df, aes(x = x, y = y)) + geom_lolli3()
#'
geom_lolli3 <- function(mapping = NULL, data = NULL, stat = "identity",
                        position = "identity", ...,
                        na.rm = FALSE, show.legend = NA,
                        inherit.aes = TRUE) {
  ggplot2::layer(
    data = data, mapping = mapping, stat = stat,
    geom = GeomLolli3, position = position,
    show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}

#' @rdname geom_lolli3
#' @format NULL
#' @usage NULL
#' @export
GeomLolli3 <- ggplot2::ggproto("GeomLolli3", ggplot2::Geom,

  required_aes = c("x", "y"),

  default_aes = ggplot2::aes(
    colour = "grey35", shape = 19, size = 1,
    fill = NA, alpha = NA, stroke = 0.5
  ),

  draw_group = function(data, panel_params, coord, na.rm = FALSE) {
    # Stems
    stem_data <- data.frame(
      x = data$x, xend = data$x,
      y = 0, yend = data$y
    )
    stem_data <- coord$transform(stem_data, panel_params)

    stems <- grid::segmentsGrob(
      x0 = stem_data$x, y0 = stem_data$y,
      x1 = stem_data$xend, y1 = stem_data$yend,
      gp = grid::gpar(col = data$colour, lwd = 0.5 * ggplot2::.pt)
    )

    # Points
    pt_data <- coord$transform(data, panel_params)
    points <- grid::pointsGrob(
      x = pt_data$x, y = pt_data$y,
      pch = data$shape, size = grid::unit(data$size, "mm"),
      gp = grid::gpar(
        col = data$colour,
        fill = scales::alpha(data$fill, data$alpha)
      )
    )

    grid::grobTree(stems, points)
  }
)
