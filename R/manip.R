#' Functions to get and set attributes of LipidomicsExperiment objects
#'
#' @param data LipidomicsExperiment object.
#' @param measure Which measure to get / set attributes of.
#' @param val Value to be assigned to the attribute.
#'
#' @return Modified LipidomicsExperiment.
#' @rdname set_attr
#' @export
#' @examples
#' data(data_normalized)
#' is_logged(data_normalized, "Area")
#' is_summarized(data_normalized)
is_logged <- function(data, measure) {
  assay_annot <- mcols(assays(data), use.names = TRUE)
  !is.null(assay_annot) && TRUE %in% assay_annot[measure, "logged"]
}

#' @rdname set_attr
#' @export
set_logged <- function(data, measure, val) {
  if (is.null(mcols(assays(data)))) {
    mcols(assays(data)) <- list(logged = FALSE, normalized = FALSE)
  }
  mcols(assays(data), use.names = TRUE)[measure, "logged"] <- val
  data
}

#' @rdname set_attr
#' @export
is_normalized <- function(data, measure) {
  assay_annot <- mcols(assays(data), use.names = TRUE)
  !is.null(assay_annot) && TRUE %in% assay_annot[measure, "normalized"]
}

#' @rdname set_attr
#' @export
set_normalized <- function(data, measure, val) {
  if (is.null(mcols(assays(data)))) {
    mcols(assays(data)) <- list(logged = FALSE, normalized = FALSE)
  }
  mcols(assays(data), use.names = TRUE)[measure, "normalized"] <- val
  data
}

#' @rdname set_attr
#' @export
is_summarized <- function(data) {
  metadata(data)$summarized
}

#' @rdname set_attr
#' @export
set_summarized <- function(data, val) {
  metadata(data)$summarized <- val
  data
}

#' Get a list of molecules that couldn't be parsed by `lipidr`
#'
#' @param data LipidomicsExperiment object.
#'
#' @return A character vector of the molecule names that could not be parsed.
#'
#' @export
#' @examples
#' data(data_normalized)
#' non_parsed_molecules(data_normalized)
non_parsed_molecules <- function(data) {
  unique(rowData(data)$Molecule[rowData(data)$not_matched])
}

#' Remove molecules that couldn't be parsed by `lipidr` from the dataset
#'
#' @param data LipidomicsExperiment object.
#'
#' @return A filtered LipidomicsExperiment object.
#'
#' @export
#' @examples
#' data(data_normalized)
#' remove_non_parsed_molecules(data_normalized)
remove_non_parsed_molecules <- function(data) {
  data[!rowData(data)$not_matched, ]
}


#' Rename molecules in a dataset.
#'
#' This function enables users to rename selected molecules in the dataset,
#' so that they can be parsed correctly by `lipidr` or modify the lipid class.
#' `lipidr` automatically updates the annotation for the renamed molecules.
#'
#' @param data LipidomicsExperiment object.
#' @param old A character vector of the molecule names to be renamed.
#' @param new A character vector of the new molecule names.
#'
#' @return A LipidomicsExperiment object with molecules name and annotation
#'   updated.
#'
#' @export
#' @examples
#' data(data_normalized)
#' old_names <- rowData(data_normalized)$Molecule
#' # replace PCO with plasmenylPC
#' new_names <- sub("^LPE", "LysoPE", old_names)
#' update_molecule_names(data_normalized, old_names, new_names)
update_molecule_names <- function(data, old, new) {
  updated_names <- to_df(data, "row") %>%
    left_join(data.frame(Molecule = old, new_names = new) %>% distinct()) %>%
    mutate(Molecule = coalesce(as.character(new_names), as.character(Molecule))) %>%
    select(-new_names)

  updated_annot <- annotate_lipids(updated_names$Molecule)
  removed_cols <- colnames(updated_annot)[colnames(updated_annot) != "Molecule"]
  updated_row_data <- updated_names %>%
    select(-c(!!!rlang::syms(removed_cols))) %>%
    left_join(updated_annot)

  row_dimname <- metadata(data)$dimnames[[1]]

  rowData(data) <- toDataFrame(updated_row_data, row.names.col = row_dimname)
  data
}


#' Remove molecules with CV larger that a threshold
#'
#' @param data LipidomicsExperiment object.
#' @param cv.cutoff CV threshold (numeric).  Default is `20`.
#' @param measure Which measure used to calculate CV, usually Area (default).
#'
#' @return LipidomicsExperiment object with molecules filtered.
#' @export
#'
#' @examples
#' data(data_normalized)
#' filter_by_cv(data_normalized)
filter_by_cv <- function(data, cv.cutoff = 20, measure = "Area") {
  keep_molecules <- apply(assay(data, measure), 1, .cv) < cv.cutoff
  data[keep_molecules, ]
}


#' Used to filter out molecules that read below a density threshold.
#'
#' @param data LipidomicsExperiment object. 
#' @param fold.cutoff Fold Chain threshold (numeric). Default is '10'.
#' @param measure Which is used to measure fold chain, usually log2(Area) (default).
#'
#' @return LipidomicsExperiment object with molecules filtered
#' @export
#'
#' @examples
#' data(data_normalized)
#' filter_by_blanks(data_normalized)
filter_by_blanks <- function(data, fold.cutoff = 10, measure = "Area") {
  blanks <- lipidr:::.is_blank(data)
  blanks.molecule <- assay(data[,blanks], measure)
  blanks.molecule[is.na(blanks.molecule)] <- 0
  blanks.mean <- rowMeans(blanks.molecule)
  
  nblanks.molecule <- assay(data[,!blanks], measure)
  nblanks.mean <- rowMeans(nblanks.molecule, na.rm=T)
  
  keep.molecules <- nblanks.mean > blanks.mean*fold.cutoff
  data[keep.molecules,]
  
}
utils::globalVariables(c("new_names"))
