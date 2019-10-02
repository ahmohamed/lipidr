#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @return Result of `rhs(lhs, ...)`.
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' @examples
#' data(data_normalized)
#' data_normalized %>% filter_by_cv()
NULL
