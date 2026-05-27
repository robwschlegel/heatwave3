.onLoad <- function(libname, pkgname) {
  hw3_init_fork_safety()
}

#' Get or set the number of threads used by heatwave3
#'
#' \code{getHW3threads()} returns the current thread count.
#' \code{setHW3threads()} sets it and returns the new value.
#'
#' The default is 50\% of available cores, which can be overridden by the
#' environment variable \code{R_HEATWAVE3_NUM_THREADS}. The per-function
#' \code{n_threads} parameter, when greater than zero, takes precedence
#' over the package-level setting for that call.
#'
#' @param threads Integer. Number of threads to use. \code{0} resets to
#'   the default (50\% of cores).
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
