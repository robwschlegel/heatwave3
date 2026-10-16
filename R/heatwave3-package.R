#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @useDynLib heatwave3, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats complete.cases aggregate weighted.mean
#' @importFrom utils head tail
## usethis namespace: end
NULL

# rlang's .data pronoun, used in the ggplot2 aes() calls
utils::globalVariables(".data")
