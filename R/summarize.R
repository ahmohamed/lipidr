#' Summarize transitions
#'
#' Calculate a single intensity for molecules with multiple transitions,
#' by determining the average or maximum intensity.
#'
#' @param data LipidomicsExperiment object.
#' @param method Choose to summarize multiple transitions by taking the average
#'   or maximum intensity. Default is `max`
#'
#' @return A LipidomicsExperiment object with single intensities per lipid molecule
#' @export
#'
#' @examples
#' datadir <- system.file("extdata", package = "lipidr")
#' filelist <- list.files(datadir, "data.csv", full.names = TRUE)
#' d <- read_skyline(filelist)
#' clinical_file <- system.file("extdata", "clin.csv", package = "lipidr")
#' d <- add_sample_annotation(d, clinical_file)
#' d_summarized <- summarize_transitions(d, method = "average")
summarize_transitions <- function(data, method = c("max", "average")) {
  stopifnot(inherits(data, "LipidomicsExperiment"))
  validObject(data)
  if (is_summarized(data)) {
    stop("data is already summarized")
  }

  method <- match.arg(method)
  multi_transitions <- to_df(data) %>%
    mutate(Molecule = fct_inorder(Molecule)) %>%
    group_by_at(vars(-1, -matches("^Product")))

  transition_gps <- split(
    multi_transitions$TransitionId,
    multi_transitions %>% group_indices()
  )
  sum_fun <- function(x) {
    if (all(is.na(x))) return (NA)
    ifelse(method == "average", mean, max)(x, na.rm = TRUE)
  }

  assay_list <- lapply(assays(data), function(m) {
    mret <- laply(transition_gps, function(x) {
      if (length(x) == 1) {
        return(m[x, ])
      } else {
        return(apply(m[x, ], 2, sum_fun))
      }
    })
    rownames(mret) <- names(transition_gps)
    return(mret)
  })

  assay_list <- as(assay_list, "SimpleList")
  mcols(assay_list) <- mcols(assays(data))

  # TODO: set Molecule(filename) as rownames
  row_data <- multi_transitions %>%
    cbind(group_idx = group_indices(.)) %>%
    summarise(rowname = first(group_idx)) %>%
    toDataFrame()

  row_data <- row_data[ row.names(assay_list[[1]]), ]

  attrs <- metadata(data)
  attrs$summarized <- TRUE
  attrs$dimnames[[1]] <- "MoleculeId"
  LipidomicsExperiment(
    assay_list = assay_list,
    metadata = attrs,
    colData = colData(data),
    rowData = row_data
  )
}

# Defined as dimname
utils::globalVariables(c("MoleculeId"))

# colnames used internally in summarize_transitions
utils::globalVariables(c("group_idx"))
