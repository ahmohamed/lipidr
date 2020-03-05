#' Read Skyline exported files
#'
#' @param files Character vector with filepaths to
#'   Skyline exported files in CSV format.
#' @importFrom forcats fct_inorder
#' @importFrom data.table fread
#' @return LipidomicsExperiment object.
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
  if (is.data.frame(files)) {
    files <- list(dataset=files)
  }
  if (is.null(names(files))) {
    char_files <- sapply(files, is.character)
    if (any(char_files)) {
      names(files)[char_files] <- basename(unlist(files[char_files]))
    }
    if (any(!char_files)) {
      names(files)[!char_files] <- paste('dataset', seq_along(files[!char_files]))
    }
  }

  datalist <- lapply(files, .read_skyline_file) %>% .uniform_attrs()

  original_data <- datalist %>%
    mutate(Sample = fct_inorder(Sample)) %>%
    group_by(Sample) %>%
    mutate(TransitionId = seq_len(n())) %>%
    ungroup() %>%
    .copy_attr(datalist) %>%
    .to_summarized_experiment()

  .check_lipids_molecules(rowData(original_data)$Molecule)
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
#' @param data LipidomicsExperiment object.
#' @param annot_file CSV file or a data.frame with at least 2 columns,
#'   sample names & group(s).
#'
#' @return LipidomicsExperiment with sample group information.
#' @export
#'
#' @examples
#' datadir <- system.file("extdata", package = "lipidr")
#'
#' # all csv files
#' filelist <- list.files(datadir, "data.csv", full.names = TRUE)
#' d <- read_skyline(filelist)
#'
#' # Add clinical info to existing LipidomicsExperiment object
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
  if(!is.data.frame(annot_file)) {
    annot <- .read_tabular(annot_file)
  } else {
    annot <- annot_file %>% mutate_if(is.factor, as.character) %>%
      as.data.frame()
  }
  stopifnot(ncol(annot) > 1)
  .check_sample_annotation(data, annot)

  # check if any column is named "sample", otherwise take the first column
  sample_col <- grep("^Sample$", colnames(annot), ignore.case = TRUE, value = TRUE)
  if (length(sample_col) > 0) {
    sample_col <- sample_col[[1]]
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

    annot <- annot[, !colnames(annot) %in% annot_cols_exist, drop=FALSE]
  }
  annot <- distinct(annot)
  colData(data) <- col_data %>% left_join(annot, by = c(rowname = sample_col))

  data
}

###########################################################################
#' Internal method to read skyline file
#' @param file skyline exported file in CSV format
#' @importFrom tidyr gather spread separate
#' @return std data.frame
#' @noRd
.read_skyline_file <- function(file) {
  if(!is.data.frame(file)) {
    original_data <- .read_tabular(file)
  } else {
    original_data <- file
  }

  .check_skyline(original_data)

  # Check it is not an empty file
  if (!nrow(original_data)) {
    stop("file ", file, " does not have any data.")
  }
  original_data[original_data == "#N/A"] <- NA

  molecule_col <- .col_exists(original_data, col_defs$molecule_cols)[[1]]
  colnames(original_data)[[ molecule_col ]] <- "Molecule"

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
    paste0("(", measures_pattern, ")"),
    "###\\1",
    colnames(original_data)
  )

  ret <- original_data %>% mutate(TransitionId=seq_len(n()))
  ret <- ret %>%
    gather("sample.measure", "value", sample_cols) %>%
    separate(
      sample.measure,
      into = c("Sample", "measure"),
      sep = "###", extra = "merge"
    ) %>%
    mutate(Sample=trimws(Sample)) %>%
    spread(measure, value)

  if (any(colnames(ret) %in% col_defs$replicate_cols)) {
    ret <- ret[, !colnames(ret) %in% col_defs$replicate_cols]
    measures <- measures[!measures %in% col_defs$replicate_cols]
  }
  ret <- ret %>% mutate_at(vars(one_of(measures)), as.numeric)

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
      paste(cols, collapse = ", ")
    )
  }
  col_idx
}

.is_pivoted <- function(d, intensity_cols) {
  !any(colnames(d) %in% intensity_cols)
}

# colnames used internally in read.pivoted / not.pivoted
utils::globalVariables(c("sample.measure", "value", "measure"))
