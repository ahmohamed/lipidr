.quos_syms <- function(x) {
  if (rlang::is_syntactic_literal(x)) {
    return()
  } else if (is.symbol(x)) {
    return(x)
  } else if (is.call(x)) {
    return(unlist(lapply(x[TRUE][-1], .quos_syms)))
  } else if (is.pairlist(x)) {
    return(unlist(lapply(x[TRUE], .quos_syms)))
  } else if (is.list(x)) {
    return(unlist(lapply(x, .quos_syms)))
  }
}

#' @importFrom stats sd
.cv <- function(a) {
  (sd(a, na.rm = TRUE) / mean(a, na.rm = TRUE)) * 100
}

.data_internal <- function(dataset) {
  if (!exists(dataset, envir = .myDataEnv)) {
    utils::data(list = c(dataset), envir = .myDataEnv)
  }
}

#' Regex-escaping for character vectors.
#'
#' @param strings A character vector to be regex-escaped.
#' @param collapse Collapse all strings to create a single pattern (using `|`).
#'
#' @return regex-escaped string to be used for pattern matching
.as_regex <- function(strings, collapse = FALSE, prefix = "", suffix = "") {
  ret <- gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", strings)
  if (collapse) {
    ret <- paste(prefix, ret, suffix, sep = "", collapse = "|")
  } else {
    ret <- paste(prefix, ret, suffix, sep = "")
  }

  ret
}
