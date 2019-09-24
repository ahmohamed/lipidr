#' @importFrom utils read.table
.is_tabular <- function(f) {
  col_count <- function(s) {
    tryCatch(
      read.table(f, nrows = 10, sep=s) %>% ncol(),
      error = function(e){return(0)}
    )
  }
  any(
    sapply(c(",", "\t"), col_count) > 1
  )
}

.check_tabular <- function(f) {
  if(!.is_tabular(f))
    stop('Data should be tabular in CSV or tab-delimited format.')

  return(TRUE)
}

.is_skyline <- function(d){
  has_col <- function(def, exact = FALSE) {
    length(.col_exists(d, col_defs[[def]], exact_match = exact, throws = FALSE)) > 0
  }
   has_col("molecule_cols", TRUE) & has_col("intensity_cols") &
  (has_col("intensity_cols", TRUE) | has_col("replicate_cols", TRUE))
}

.check_skyline <- function(d) {
  if (.is_skyline(d)) {
    return(TRUE)
  }
  err <- 'Not a valid Skyline export.'
  has_col <- function(def, exact = FALSE) {
    length(.col_exists(d, col_defs[[def]], exact_match = exact, throws = FALSE)) > 0
  }
  if (!has_col("molecule_cols", TRUE)) {
    msg <- paste(
      err,
      'Data should have one of these columns with molecule names:',
      paste(col_defs$molecule_cols, collapse = ', ')
    )
    stop(msg)
  }
  if (!has_col("intensity_cols")) {
    msg <- paste(
      err,
      'Data should have one of these measures exported:',
      paste(col_defs$intensity_cols, collapse = ', ')
    )
    stop(msg)
  }
  if (!has_col("intensity_cols", TRUE) && !!has_col("replicate_cols", TRUE)) {
    msg <- paste(
      err,
      "In Skyline report, you should either export Replicate column",
      "or pivot by replicates"
    )
    stop(msg)
  }
}

.is_sample_annotation <- function(data, df) {
  if("Sample" %in% colnames(df)) {
    sample_col <- df$Sample
  } else {
    sample_col <- df[,1]
  }

  all(colnames(data) %in% as.character(sample_col))
}

.check_sample_annotation <- function(data, df) {
  if(!.is_sample_annotation(data, df))
    stop('All sample names must be in the first column',
      ' or a column named "Sample"')

  return(TRUE)
}

.have_lipids_molecules <- function(mols) {
  # correcting for edge case where df[,1, drop=FALSE] is passed
  mols <- unlist(mols) 
  matched <- !annotate_lipids(mols, no_match = "ignore")$not_matched
  if ((sum(matched) / length(mols)) < 0.5) {
    return(FALSE)
  }
  return(TRUE)
}
.check_lipids_molecules <- function(mols) {
  if(!.have_lipids_molecules(mols)) {
    warning('More that 50% of molecule names cannot be parsed as lipids.')
  }
}

.get_mol_dim <- function(df) {
  possible <- list(
    row_names = rownames(df),
    first_col = df[,1],
    col_names = colnames(df)
  )
  for (i in seq_along(possible)) {
    if(!is.null(possible[[i]]) && .have_lipids_molecules(possible[[i]])) {
      return(names(possible)[[i]])
    }
  }
  return('none')
}

.is_num_matrix <- function(mat) {
  if (!is.matrix(mat) || !is.numeric(c(mat))) {
    return (FALSE)
  }
  na_vals <- !is.na(mat)
  return(
    ncol(mat) > 1 &
      sum(na_vals)/length(mat) > 0.5
  )
}

col_defs <- list(
  class_cols = c("Protein Name", "Protein"),
  molecule_cols = c(
    "Peptide Name", "Peptide", "Molecule Name", "Precursor Ion Name"
  ),
  replicate_cols = c("Replicate Name", "Replicate"),
  intensity_cols = c("Area", "Height", "Area Normalized"),
  measure_cols = c(
    "Area", "Height", "Area Normalized", "Retention Time", "Background"
  )
)
