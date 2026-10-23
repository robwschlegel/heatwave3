# geom_flame3(): exceedance runs become polygons; missing values are not
# exceedances (they used to crash draw_group when nothing exceeded).

flame_ids <- function(y, y2 = 0, ...) {
  d <- data.frame(x = seq_along(y), y = y, y2 = y2, colour = NA, fill = "salmon",
                  linewidth = 0.5, linetype = 1, alpha = NA, group = 1)
  coord <- list(transform = function(data, panel_params) data)
  g <- GeomFlame3$draw_group(d, NULL, coord, ...)
  if (inherits(g, "null")) return(character(0))
  as.vector(tapply(as.numeric(g$x), as.integer(g$id), function(v) paste(range(v), collapse = "-")))
}

test_that("runs, n and n_gap select the polygons", {
  y <- c(-1, 1, 1, -1, 1, 1, 1, -1, -1, 1, 1, -1)
  expect_identical(flame_ids(y), c("2-3", "5-7", "10-11"))
  expect_identical(flame_ids(y, n = 2), "5-7")
  expect_identical(flame_ids(y, n_gap = 1), c("2-3", "4-4", "5-7", "10-11"))
  expect_identical(flame_ids(-y, reverse = TRUE), c("2-3", "5-7", "10-11"))
})

test_that("missing values are not exceedances and never error", {
  expect_identical(flame_ids(c(-1, 1, NA, 1, -1)), c("2-2", "4-4"))
  expect_identical(flame_ids(c(-1, NA, -1, -1)), character(0))
  expect_identical(flame_ids(c(1, 1, NA, 1, 1), n_gap = 1), c("1-2", "3-3", "4-5"))
})

test_that("without ggplot2 the geoms are NULL, not an error, and plotting says why", {
  # R CMD INSTALL serialises the whole namespace after loading, which evaluates
  # the delayed GeomFlame3/GeomLolli3; an error there blocked installation.
  local_mocked_bindings(.has_ggplot2 = function() FALSE)
  expect_null(.geom_if_ggplot2(.make_geom_flame3))
  expect_null(.geom_if_ggplot2(.make_geom_lolli3))
  expect_error(.need_ggplot2(), "Plotting needs the ggplot2 package")
})

test_that("with ggplot2 the geoms are built on first use", {
  skip_if_not_installed("ggplot2")
  expect_s3_class(.geom_if_ggplot2(.make_geom_flame3), "Geom")
  expect_s3_class(GeomLolli3, "Geom")
})

test_that("geom_lolli3 draws one stem and one point per row", {
  skip_if_not_installed("ggplot2")
  d <- data.frame(x = 1:4, y = c(1.2, 2.1, 0.8, 1.5), colour = "grey35", shape = 19,
                  size = 1, fill = NA, alpha = NA, stroke = 0.5, group = 1)
  coord <- list(transform = function(data, panel_params) data)
  g <- GeomLolli3$draw_group(d, NULL, coord)
  stems <- g$children[[1]]; points <- g$children[[2]]
  expect_s3_class(stems, "segments")
  expect_s3_class(points, "points")
  # Stems run from y = 0 up to each value, at each point's x
  expect_equal(as.numeric(stems$x0), d$x)
  expect_equal(as.numeric(stems$x1), d$x)
  expect_equal(as.numeric(stems$y0), rep(0, 4))
  expect_equal(as.numeric(stems$y1), d$y)
  expect_equal(as.numeric(points$y), d$y)
})
