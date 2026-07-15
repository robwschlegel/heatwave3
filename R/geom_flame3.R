#' Flame geom for MHW visualisation
#'
#' A ggplot2 geom that fills the area between two lines (temperature and
#' threshold) where the first exceeds the second, creating "flame" polygons.
#' For cold-spells, set \code{reverse = TRUE} to fill where \code{y < y2}.
#'
#' @param mapping Set of aesthetic mappings. Requires \code{x}, \code{y}, and \code{y2}.
#'   An optional \code{depth} aesthetic can be mapped for depth-resolved data
#'   (see Details).
#' @param data The data to display.
#' @param stat The statistical transformation. Default \code{"identity"}.
#' @param position Position adjustment. Default \code{"identity"}.
#' @param ... Additional arguments passed to the layer.
#' @param n Minimum number of consecutive exceeding points. Default \code{0}.
#' @param n_gap Maximum gap to bridge between exceeding segments. Default \code{0}.
#' @param reverse Logical. If \code{TRUE}, fill where \code{y < y2} (for
#'   cold-spells). Default \code{FALSE} (fill where \code{y > y2}).
#' @param na.rm Remove NAs? Default \code{FALSE}.
#' @param show.legend Show legend? Default \code{NA}.
#' @param inherit.aes Inherit aesthetics from the plot? Default \code{TRUE}.
#'
#' @return A ggplot2 layer.
#'
#' @details
#' \code{geom_flame3()} finds exceedance runs by walking each group's rows in
#' order, so a single group must be one continuous series through \code{x}.
#' For depth-resolved data (e.g. \code{hw3_export()} on a
#' \code{ts2clm3(depth_range = ...)} product), map \code{depth} as an
#' aesthetic -- \code{aes(x = t, y = temp, y2 = thresh, depth = depth)} --
#' and it is automatically folded into the row grouping, so each depth level
#' gets its own exceedance runs instead of being treated as one series with
#' depths interleaved. This happens even though \code{depth} is a plain
#' numeric column, which ggplot2's own default grouping would otherwise merge
#' into a single group (unlike a discrete aesthetic such as \code{colour}).
#'
#' @export
#'
#' @examples
#' library(ggplot2)
#' # Simple example with synthetic data
#' df <- data.frame(
#'   x = 1:100,
#'   y = sin(seq(0, 4 * pi, length.out = 100)) + 2,
#'   y2 = 1.5
#' )
#' ggplot(df, aes(x = x, y = y, y2 = y2)) +
#'   geom_flame3() +
#'   geom_line()
#'
geom_flame3 <- function(mapping = NULL, data = NULL, stat = "identity",
                        position = "identity", ..., n = 0, n_gap = 0,
                        reverse = FALSE,
                        na.rm = FALSE, show.legend = NA,
                        inherit.aes = TRUE) {
  ggplot2::layer(
    data = data, mapping = mapping, stat = stat,
    geom = GeomFlame3, position = position,
    show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(n = n, n_gap = n_gap, reverse = reverse, na.rm = na.rm, ...)
  )
}

#' @rdname geom_flame3
#' @format NULL
#' @usage NULL
#' @export
GeomFlame3 <- ggplot2::ggproto("GeomFlame3", ggplot2::Geom,

  required_aes = c("x", "y", "y2"),

  default_aes = ggplot2::aes(
    colour = NA, fill = "salmon", linewidth = 0.5,
    linetype = 1, alpha = NA
  ),

  extra_params = c("na.rm", "n", "n_gap", "reverse"),

  # A mapped `depth` column is numeric, so ggplot2's own default grouping
  # (which only folds in discrete aesthetics) leaves it out -- rows from
  # different depths would land in one group and get RLE'd together as if
  # they were one continuous series. Re-derive group to include depth
  # whenever it's present, preserving each existing group's row order.
  setup_data = function(data, params) {
    if (!is.null(data$depth)) {
      data$group <- match(
        interaction(data$group, data$depth, drop = TRUE, lex.order = TRUE),
        unique(interaction(data$group, data$depth, drop = TRUE, lex.order = TRUE))
      )
    }
    data
  },

  draw_group = function(data, panel_params, coord,
                        n = 0, n_gap = 0, reverse = FALSE, na.rm = FALSE) {
    if (reverse) {
      exceed <- data$y < data$y2
    } else {
      exceed <- data$y > data$y2
    }

    rle_exc <- rle(exceed)
    ends <- cumsum(rle_exc$lengths)
    starts <- c(1L, ends[-length(ends)] + 1L)

    keep <- rle_exc$values & rle_exc$lengths >= max(1, n + 1)

    if (n_gap > 0 && sum(keep) > 1) {
      for (i in seq_along(keep)[-1]) {
        if (!keep[i] && !rle_exc$values[i] && rle_exc$lengths[i] <= n_gap) {
          if (i > 1 && keep[i - 1] && i < length(keep) && keep[i + 1]) {
            keep[i] <- TRUE
          }
        }
      }
    }

    if (!any(keep)) return(grid::nullGrob())

    polys <- list()
    seg_id <- 0
    for (i in which(keep)) {
      s <- starts[i]; e <- ends[i]
      seg <- data[s:e, , drop = FALSE]

      if (reverse) {
        upper <- data.frame(x = seg$x, y = seg$y2)
        lower <- data.frame(x = rev(seg$x), y = rev(seg$y))
      } else {
        upper <- data.frame(x = seg$x, y = seg$y)
        lower <- data.frame(x = rev(seg$x), y = rev(seg$y2))
      }
      poly <- rbind(upper, lower)

      seg_id <- seg_id + 1
      poly$group <- seg_id
      poly$colour <- seg$colour[1]
      poly$fill <- seg$fill[1]
      poly$linewidth <- seg$linewidth[1]
      poly$linetype <- seg$linetype[1]
      poly$alpha <- seg$alpha[1]

      polys[[seg_id]] <- poly
    }

    poly_data <- do.call(rbind, polys)
    poly_data <- coord$transform(poly_data, panel_params)

    grid::polygonGrob(
      x = poly_data$x, y = poly_data$y,
      id = poly_data$group,
      gp = grid::gpar(
        col = poly_data$colour[!duplicated(poly_data$group)],
        fill = scales::alpha(
          poly_data$fill[!duplicated(poly_data$group)],
          poly_data$alpha[!duplicated(poly_data$group)]
        ),
        lwd = poly_data$linewidth[!duplicated(poly_data$group)] * ggplot2::.pt,
        lty = poly_data$linetype[!duplicated(poly_data$group)]
      )
    )
  }
)
