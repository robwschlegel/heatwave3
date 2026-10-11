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
