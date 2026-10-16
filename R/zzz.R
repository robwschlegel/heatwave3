.onLoad <- function(libname, pkgname) {
  hw3_init_fork_safety()
  # The ggplot2 geoms are created when first used: ggplot2 is only suggested,
  # and loading it would add most of a second to every library(heatwave3).
  # Without ggplot2 they are NULL rather than an error: R CMD INSTALL's
  # final-location test serialises the whole namespace, which evaluates these
  # promises, and an error there made heatwave3 uninstallable without ggplot2.
  # geom_flame3() and geom_lolli3() call .need_ggplot2() before using them.
  ns <- topenv()
  if (!environmentIsLocked(ns)) {
    delayedAssign("GeomFlame3", .geom_if_ggplot2(.make_geom_flame3), assign.env = ns)
    delayedAssign("GeomLolli3", .geom_if_ggplot2(.make_geom_lolli3), assign.env = ns)
  }
}

.has_ggplot2 <- function() requireNamespace("ggplot2", quietly = TRUE)

# The geom made by make(), or NULL when ggplot2 is not installed
.geom_if_ggplot2 <- function(make) if (.has_ggplot2()) make() else NULL

.need_ggplot2 <- function() {
  if (!.has_ggplot2())
    stop("Plotting needs the ggplot2 package: install.packages(\"ggplot2\")",
         call. = FALSE)
}

#' Get or set the number of threads used by heatwave3
#'
#' \code{getHW3threads()} returns the current thread count.
#' \code{setHW3threads()} sets it and returns the new value.
#'
#' The default is 50 percent of available cores, which can be overridden by the
#' environment variable \code{R_HEATWAVE3_NUM_THREADS}. The per-function
#' \code{n_threads} parameter, when greater than zero, takes precedence
#' over the package-level setting for that call.
#'
#' @param threads Integer. Number of threads to use. \code{0} resets to
#'   the default (50 percent of cores).
#'
#' @return Integer: the current (or newly set) thread count.
#'
#' @examples
#' getHW3threads()
#' old <- setHW3threads(2)
#' getHW3threads()
#' setHW3threads(0)  # reset to default
#'
#' @export
getHW3threads <- function() {
  hw3_get_threads()
}

#' @rdname getHW3threads
#' @export
setHW3threads <- function(threads = 0L) {
  hw3_set_threads(as.integer(threads))
}

# n_threads argument -> C++: NULL means the session default (0 tells the C++
# core to use getHW3threads()).
.hw3_nthreads <- function(n_threads) {
  if (is.null(n_threads)) 0L else as.integer(n_threads)
}

# compression argument -> C++: one whole number from 0 to 9
.hw3_compression <- function(compression) {
  if (!is.numeric(compression) || length(compression) != 1 || is.na(compression) ||
      compression != round(compression) || compression < 0 || compression > 9)
    stop("compression must be a whole number from 0 to 9.", call. = FALSE)
  as.integer(compression)
}
