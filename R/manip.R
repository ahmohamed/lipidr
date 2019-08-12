#' Functions to get and set attributes of SkylineExperiment objects
#'
#' @param data SkylineExperiment object created by [read_skyline()].
#' @param measure Which measure to get / set attributes of.
#' @param val Value to be assigned to the attribute.
#'
#' @return Modified SkylineExperiment.
#' @rdname set_attr
#' @export
is_logged <- function(data, measure) {
  mcols(assays(data), use.names = TRUE)[measure, "logged"]
}

#' @rdname set_attr
#' @export
set_logged <- function(data, measure, val) {
  mcols(assays(data), use.names = TRUE)[measure, "logged"] <- val
  data
}

#' @rdname set_attr
#' @export
is_normalized <- function(data, measure) {
  mcols(assays(data), use.names = TRUE)[measure, "normalized"]
}

#' @rdname set_attr
#' @export
set_normalized <- function(data, measure, val) {
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
#' @param data SkylineExperiment object created by [read_skyline()].
#'
#' @return A character vector of the molecule namaes that could not be parsed.
#'
#' @export
non_parsed_molecules <- function(data) {
  unique(rowData(data)$Molecule[rowData(data)$not_matched])
}

#' Remove molecules that couldn't be parsed by `lipidr` from the dataset
#'
#' @param data SkylineExperiment object created by [read_skyline()].
#'
#' @return A filtered SkylineExperiment object.
#'
#' @export
remove_non_parsed_molecules <- function(data) {
  data[!rowData(data)$not_matched, ]
}


#' Rename molecules in a dataset
#' This function enables users to rename selected molecules in the dataset,
#' so that they can be parsed correctly by `lipidr` or modify the lipidclass.
#' `lipidr` automatically updates the annotation for the renamed molecules.
#'
#' @param data SkylineExperiment object created by [read_skyline()].
#' @param old A character vector of the molecule names to be renamed.
#' @param new A character vector of the new molecule names.
#'
#' @return A SkylineExperiment object with molecules name and annotation
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
    left_join(data.frame(Molecule=old, new_names=new) %>% distinct()) %>%
    mutate(Molecule=coalesce(as.character(new_names), Molecule)) %>%
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
