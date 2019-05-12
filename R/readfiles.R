#' Read Skyline exported files
#'
#' @param files Character vector with filepaths to 
#'   Skyline exported files in CSV format.
#' @importFrom forcats fct_inorder
#' @return SkylineExperiment object.
#' @export
#'
#' @examples
#' datadir <- system.file("extdata", package = "lipidr")
#'
#' # all csv files
#' filelist <- list.files(datadir, "data.csv", full.names = TRUE)
#' d <- read_skyline(filelist)
#'
#' # View automatically generated lipid annotations
#' rowData(d)
read_skyline <- function(files) {
  names(files) <- basename(files)
  datalist <- lapply(files, .read_skyline_file) %>% .uniform_attrs()

  original_data <- datalist %>%
    mutate(Sample = fct_inorder(Sample)) %>%
    group_by(Sample) %>%
    mutate(TransitionId = seq_len(n())) %>%
    ungroup() %>%
    .copy_attr(datalist) %>%
    .to_summarized_experiment()

  message(
    "Successfully read ", length(files), " methods.\n",
    "Your data contain ", ncol(original_data), " samples, ",
    length(unique(rowData(original_data)$Class)), " lipid classes, ",
    length(unique(rowData(original_data)$Molecule)), " lipid molecules."
  )

  original_data
}

#' Add sample annotation to Skyline data frame
#'
#' @param data SkylineExperiment object created by [read_skyline()].
#' @param annot_file CSV file with at least 2 columns, sample names & group(s).
#'
#' @return Skyline data.frame with sample group information.
#' @export
#'
#' @examples
#' datadir <- system.file("extdata", package = "lipidr")
#'
#' # all csv files
#' filelist <- list.files(datadir, "data.csv", full.names = TRUE)
#' d <- read_skyline(filelist)
#'
#' # Add clinical info to existing SkylineExperiment object
#' clinical_file <- system.file("extdata", "clin.csv", package = "lipidr")
#' d <- add_sample_annotation(d, clinical_file)
#' colData(d)
#' d$group
#'
#' # Subset samples using clinical information
#' # Note we are subsetting columns
#' d[, d$group == "QC"]
#'
#' # Subset lipids using lipid annotation
#' # Note we are subsetting rows
#' d[rowData(d)$istd, ]
add_sample_annotation <- function(data, annot_file) {
  annot <- read.csv(annot_file)
  stopifnot(ncol(annot) > 1)

  # check if any column is named "sample", otherwise take the first column
  sample_col <- grep("Sample", colnames(annot), ignore.case = TRUE)
  if (length(sample_col) > 0) {
    sample_col <- colnames(annot)[sample_col][[1]]
  } else {
    warning("No column named 'Sample', taking the first column as sample names")
    sample_col <- colnames(annot)[[1]]
  }
  annot_cols <- colnames(annot)[colnames(annot) != sample_col]
  col_data <- colData(data)

  if (any(annot_cols %in% colnames(col_data))) {
    annot_cols_exist <- annot_cols[annot_cols %in% colnames(col_data)]
    warning(
      "These annotation columns already exist in data: ",
      paste(annot_cols_exist, collapse = ", ")
    )
    annot_cols <- annot_cols[!annot_cols %in% annot_cols_exist]

    if (length(annot_cols) == 0) return(data)

    annot <- annot[, !colnames(annot) %in% annot_cols_exist]
  }

  colData(data) <- col_data %>% left_join(annot, by = c(rowname = sample_col))

  data
}

###########################################################################
#' Internal method to read skyline file
#' @param file skyline exported file in CSV format
#' @importFrom utils read.csv
#' @importFrom tidyr gather spread separate
#' @return std data.frame
.read_skyline_file <- function(file) {
  original_data <- read.csv(file, stringsAsFactors = FALSE)
  original_data[original_data == "#N/A"] <- NA

  col_defs <- list(
    class_cols = c("Protein.Name", "Protein"),
    molecule_cols = c(
      "Peptide.Name", "Peptide", "Molecule.Name", "Precursor.Ion.Name"
    ),
    replicate_cols = c("Replicate.Name", "Replicate"),
    intensity_cols = c("Area", "Height", "Area.Normalized"),
    measure_cols = c(
      "Area", "Height", "Area.Normalized", "Retention.Time", "Background"
    )
  )

  .col_exists(original_data, col_defs$class_cols)
  molecule_col <- .col_exists(original_data, col_defs$molecule_cols)[[1]]
  colnames(original_data)[[ molecule_col ]] <- "Molecule"

  .col_exists(original_data, col_defs$intensity_cols)

  if (.is_pivoted(original_data, col_defs$intensity_cols)) {
    return(.read_pivoted(original_data, col_defs))
  } else {
    return(.read_not_pivoted(original_data, col_defs))
  }
}

.read_not_pivoted <- function(original_data, col_defs) {
  intensity_cols <- .col_exists(
    original_data, col_defs$intensity_cols,
    exact_match = TRUE
  )
  replicate_cols <- .col_exists(
    original_data, col_defs$replicate_cols,
    exact_match = TRUE, throws = FALSE
  )
  if (length(replicate_cols) == 0) {
    stop(
      "In Skyline report, you should either export Replicate column ",
      "or pivot by replicates"
    )
  }

  colnames(original_data)[[ replicate_cols[[1]] ]] <- "Sample"
  measure_cols <- .col_exists(
    original_data, col_defs$measure_cols,
    exact_match = TRUE, throws = FALSE
  )

  ret <- original_data
  ret[, measure_cols] <- sapply(ret[, measure_cols], as.numeric)
  attr(ret, "skyline") <- list(
    skyline = TRUE,
    measures = colnames(original_data)[measure_cols],
    intensity_cols = colnames(original_data)[intensity_cols]
  )

  ret
}

.read_pivoted <- function(original_data, col_defs) {
  intensity_cols <- .col_exists(original_data, col_defs$intensity_cols)
  intensity_colnames <- colnames(original_data)[intensity_cols]

  # Extract sample names from intensity columns
  sample_names <- unique(sub(
    "(.*)(Area|Height|Area\\.Normalized)$", "\\1", intensity_colnames
  ))
  samples_pattern <- .as_regex(sample_names, prefix = "^", collapse = TRUE)
  sample_cols <- grep(samples_pattern, colnames(original_data))

  # Extract all measure names from sample columns
  measures <- unique(sub(
    samples_pattern, "", colnames(original_data)[sample_cols]
  ))
  measures_pattern <- .as_regex(measures, collapse = TRUE)

  colnames(original_data) <- sub(
    paste0("\\.(", measures_pattern, ")"),
    "###\\1",
    colnames(original_data)
  )

  ret <- original_data
  ret[, sample_cols] <- sapply(ret[, sample_cols], as.numeric)
  ret <- ret %>%
    gather("sample.measure", "value", sample_cols) %>%
    separate(
      sample.measure,
      into = c("Sample", "measure"),
      sep = "###", extra = "merge"
    ) %>%
    spread(measure, value)

  attr(ret, "skyline") <- list(
    skyline = TRUE,
    measures = measures,
    intensity_cols = intensity_colnames
  )
  ret
}

.col_exists <- function(d, cols, exact_match = FALSE, throws = TRUE) {
  if (exact_match) {
    col_idx <- which(colnames(d) %in% cols)
  } else {
    col_idx <- grep(.as_regex(cols, collapse = TRUE), colnames(d))
  }

  if (throws && length(col_idx) == 0) {
    stop(
      "At least one of these columns should be exported from Skyline report.",
      cols
    )
  }
  col_idx
}

.is_pivoted <- function(d, intensity_cols) {
  !any(colnames(d) %in% intensity_cols)
}

# colnames used internally in read.pivoted / not.pivoted
utils::globalVariables(c("sample.measure", "value", "measure"))
