#' Detect marine heatwaves day by day in arrays
#'
#' The in-memory counterpart of \code{\link{detect_event3}} with
#' \code{daily = "only"}: the same detection, but on arrays instead of files,
#' and with the climatology already matched to each day. It returns only
#' per-day results, because event metrics such as duration or cumulative
#' intensity need the whole event, which a short series (a forecast window,
#' for example) may cut off.
#'
#' @param temp Numeric vector, matrix or array of daily temperatures with time
#'   on the \strong{first} dimension (a vector is one series) and any number
#'   of further dimensions (pixels). The days must be consecutive, with no
#'   gaps. Missing days are \code{NA}.
#' @param seas,thresh Seasonal climatology and threshold for each day of
#'   \code{temp}, with the same dimensions. Match them to the dates
#'   beforehand, for example from a \code{\link{ts2clm3}} climatology by day
#'   of year.
#' @param minDuration Minimum duration (days) for an event. Default \code{5}.
#' @param joinAcrossGaps Logical. Join events separated by short gaps?
#'   Default \code{TRUE}.
#' @param maxGap Maximum gap length (days) to join across. Default \code{2}.
#' @param coldSpells Logical. Detect cold-spells instead of heatwaves?
#'   Default \code{FALSE}.
#' @param roundRes Decimal places for \code{intensity}. Default \code{4}.
#' @param outputs Which results to return, in this order: any of
#'   \code{"event"}, \code{"event_no"}, \code{"intensity"} and
#'   \code{"category"}.
#' @param n_threads Number of threads for the per-pixel computation. Default
#'   \code{NULL}: the session default (see \code{\link{setHW3threads}}). Use
#'   \code{1} when the calls themselves run in parallel.
#'
#' @return One array per name in \code{outputs}, each with the dimensions
#'   (and dimnames) of \code{temp}: the array itself when one output is
#'   requested, otherwise a named list in the order of \code{outputs}.
#'   \describe{
#'     \item{event}{Logical: the day belongs to an event.}
#'     \item{event_no}{Integer: the event's number within its pixel, counted
#'       from the start of this call (0 outside events). It is not an
#'       identifier across calls: the same event in two overlapping series
#'       can have different numbers.}
#'     \item{intensity}{Numeric: \code{temp - seas} on every day with data,
#'       rounded to \code{roundRes}; \code{NA} where either is missing.}
#'     \item{category}{Integer: the daily Hobday et al. (2018) category, 1 (I
#'       Moderate) to 4 (IV Extreme), on event days beyond the threshold;
#'       0 on all other days.}
#'   }
#'   A pixel with no temperature or no \code{seas} at all gets \code{FALSE},
#'   0, \code{NA} and 0. Results are identical to the per-day file of
#'   \code{\link{detect_event3}} for the same days and climatology.
#'
#' @details Whether a day belongs to an event depends on the days around it,
#'   so events are cut off at both ends of the series: a run of exceedances
#'   that started before the first day is shorter here, and one that
#'   continues past the last day may not reach \code{minDuration}. Prepend
#'   earlier days when that matters. \code{intensity} does not depend on the
#'   other days, and \code{category} only through event membership.
#'
#'   Time comes first because R stores arrays column by column, so each
#'   pixel's series is then contiguous and is read without copying. Gridded
#'   data from packages such as \pkg{terra} or \pkg{stars} can be passed as
#'   plain arrays: for example \code{t(terra::values(r))} for a SpatRaster
#'   \code{r} with one layer per day, or \code{aperm(x, c(3, 1, 2))} for a
#'   \code{[x, y, time]} array.
#'
#' @seealso \code{\link{detect_event3}}, \code{\link{ts2clm3}}
#'
#' @export
#'
#' @examples
#' # Two pixels, 30 days, a constant climatology
#' temp <- matrix(20, nrow = 30, ncol = 2)
#' temp[6:15, 1] <- 25
#' seas <- matrix(20, nrow = 30, ncol = 2)
#' thresh <- matrix(21, nrow = 30, ncol = 2)
#' res <- detect_event_array3(temp, seas, thresh)
#' which(res$event[, 1])
#' table(res$category)
detect_event_array3 <- function(temp, seas, thresh,
                                minDuration = 5,
                                joinAcrossGaps = TRUE,
                                maxGap = 2,
                                coldSpells = FALSE,
                                roundRes = 4,
                                outputs = c("event", "intensity", "category"),
                                n_threads = NULL) {

  valid <- c("event", "event_no", "intensity", "category")
  if (!is.character(outputs) || length(outputs) == 0 || anyDuplicated(outputs) ||
      !all(outputs %in% valid)) {
    stop("'outputs' must name one or more of ", paste(valid, collapse = ", "),
         ", each at most once.", call. = FALSE)
  }

  dims_of <- function(x) if (is.null(dim(x))) length(x) else dim(x)
  dims <- dims_of(temp)
  if (!identical(as.integer(dims_of(seas)), as.integer(dims)) ||
      !identical(as.integer(dims_of(thresh)), as.integer(dims))) {
    stop("temp, seas and thresh must have the same dimensions.", call. = FALSE)
  }

  res <- hw3_detect_events_array(
    temp, seas, thresh,
    ntime = as.integer(dims[1]),
    minDuration = as.integer(minDuration),
    joinAcrossGaps = as.logical(joinAcrossGaps),
    maxGap = as.integer(maxGap),
    coldSpells = as.logical(coldSpells),
    roundRes = as.integer(roundRes),
    n_threads = .hw3_nthreads(n_threads),
    event = "event" %in% outputs,
    event_no = "event_no" %in% outputs,
    intensity = "intensity" %in% outputs,
    category = "category" %in% outputs
  )

  res <- lapply(res[outputs], function(x) {
    if (!is.null(dim(temp))) {
      dim(x) <- dim(temp)
      dimnames(x) <- dimnames(temp)
    }
    x
  })
  if (length(outputs) == 1) res[[1]] else res
}
