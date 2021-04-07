#' Differential analysis of lipids between sample groups
#'
#' `de_analysis` and `de_design` perform differential analysis of measured
#' lipids that are associated with a sample group (annotation). `de_analysis`
#' accepts a list of contrasts, while `de_design` allows users to define a
#' design matrix, useful for complex experimental designs or for adjusting
#' possible confounding variables.
#'
#' @param data LipidomicsExperiment object,
#'   should be normalized and log2 transformed.
#' @param ... Expressions, or character strings which can be parsed to
#'   expressions, specifying contrasts. These are passed to
#'   `limma::makeContrasts`.
#' @param measure Which measure to use as intensity, usually Area (default).
#' @param group_col Name of the column containing sample groups. If not
#'   provided, defaults to first sample annotation column.
#'
#' @importFrom forcats fct_drop
#' @importFrom rlang quos
#' @importFrom stats model.matrix setNames
#' @return TopTable as returned by limma package
#' @export
#'
#' @examples
#' # type ?normalize_pqn to see how to normalize and log2-transform your data
#' data(data_normalized)
#'
#' # Specifying contrasts
#' de_results <- de_analysis(
#'   data_normalized,
#'   HighFat_water - NormalDiet_water,
#'   measure = "Area"
#' )
de_analysis <- function(data, ..., measure = "Area", group_col = NULL) {
  if (is.null(group_col)) {
    if (ncol(colData(data)) > 0) {
      group_col <- names(colData(data))[[1]]
    } else {
      stop("Please add clinical data or specify a group column")
    }
  }

  symbols <- as.character(.quos_syms(quos(...)))

  group <- colData(data)[[group_col]]
  if (!all(symbols %in% as.character(group))) {
    stop(
      "These constrast variables are not present in ", group_col, " column: ",
      paste(symbols[!symbols %in% as.character(group)], collapse=", ")
    )
  }
  data <- data[, group %in% symbols]

  group <- fct_drop(colData(data)[[group_col]])
  design <- model.matrix(~ 0 + group)
  colnames(design) <- gsub("group", "", colnames(design))
  return(de_design(data = data, design = design, ..., measure = measure))
}

#' @param design Design matrix generated from [model.matrix()],
#'   or a design formula.
#' @param coef Column number or column name specifying which coefficient of
#'   the linear model is of interest.
#'
#' @importFrom rlang is_formula
#' @importFrom limma topTable lmFit makeContrasts contrasts.fit eBayes
#' @export
#' @rdname de_analysis
#' @examples
#' # Using formula
#' de_results_formula <- de_design(
#'   data = data_normalized,
#'   design = ~group,
#'   coef = "groupHighFat_water",
#'   measure = "Area"
#' )
#'
#' # Using design matrix
#' design <- model.matrix(~group, data = colData(data_normalized))
#' de_results_design <- de_design(
#'   data = data_normalized,
#'   design = design,
#'   coef = "groupHighFat_water",
#'   measure = "Area"
#' )
de_design <- function(data, design, ..., coef = NULL, measure = "Area") {
  if (is_formula(design)) {
    design <- model.matrix(design, data = colData(data))
    if (!identical(colnames(data), rownames(design))) {
      warning(
        'These samples are not present in the design matrix: ',
        paste(colnames(data) [!colnames(data) %in% rownames(design)], collapse = ", "),
        '. Possibly because the grouping columns have missing values.'
      )
      data <- data[, rownames(design)]
    }
  }
  if (!is.matrix(design)) {
    stop("design should be a matrix or formula")
  }
  if (!limma::is.fullrank(design)) {
    stop("Tested variables are redundant (Design matrix is not full rank).")
  }

  vfit <- lmFit(assay(data, measure), design)

  if (is.null(coef)) {
    if (length(quos(...)) == 0) {
      warning(
        "No contrasts or coefficients are provided. ",
        "ANOVA-style analysis will be performed using all group."
      )
      # Exclude the first column (intercept)
      coef <- list("ANOVA" = seq(2, ncol(design)))
    } else {
      contr.matrix <- limma::makeContrasts(..., levels = colnames(design))
      vfit <- limma::contrasts.fit(vfit, contrasts = contr.matrix)
      coef <- setNames(seq_len(ncol(contr.matrix)), colnames(contr.matrix))
    }
  } else {
    if (!coef %in% colnames(design)) {
      stop(
        "One or more coefficients is not in the design matrix.",
        " Allowed values are ", colnames(design)
      )
    }
    names(coef) <- coef
  }

  efit <- eBayes(vfit)
  dimname_x <- metadata(data)$dimnames[[1]]

  top <- lapply(
    coef, function(x)
      topTable(efit, number = Inf, coef = x) %>% rownames_to_column(dimname_x)
  ) %>%
    bind_rows(.id = "contrast")

  top <- to_df(data, dim = "row") %>%
    select(
      one_of("Molecule", "Class", "total_cl", "total_cs", "istd", dimname_x)
    ) %>%
    .left_join_silent(top)
  attr(top, 'measure') <- measure
  return(top)
}

#' @describeIn de_analysis gets a list of significantly changed lipids for
#' each contrast.
#'
#' @param de.results Output of [de_analysis()].
#' @param p.cutoff Significance threshold.  Default is `0.05`.
#' @param logFC.cutoff Cutoff limit for log2 fold change.  Default is `1`.
#'   Ignored in multi-group (ANOVA-style) comparisons.
#'
#' @return `significant_molecules` returns a character vector with names of
#'   significantly differentially changed lipids.
#' @export
#'
#' @examples
#' significant_molecules(de_results)
significant_molecules <- function(de.results, p.cutoff = 0.05,
                                  logFC.cutoff = 1) {
  if (!"logFC" %in% colnames(de.results)) {
    message(
      "de.results contains ANOVA-style comparison.",
      " LogFC cutoff will be ignored"
    )
    ret <- de.results %>%
      filter(adj.P.Val < p.cutoff) %>%
      (function(x) split(x$Molecule, x$contrast))
    return(ret)
  }

  de.results %>%
    filter(adj.P.Val < p.cutoff, abs(logFC) > logFC.cutoff) %>%
    (function(x) split(x$Molecule, x$contrast))
}

#' @describeIn de_analysis plots a volcano chart for differential analysis
#' results.
#'
#' @param show.labels Whether labels should be displayed for
#'   significant lipids.  Default is `TRUE`.
#'
#' @return `plot_results_volcano` returns a ggplot object.
#' @export
#' @examples
#' plot_results_volcano(de_results, show.labels = FALSE)
plot_results_volcano <- function(de.results, show.labels = TRUE) {
  if (!"logFC" %in% colnames(de.results)) {
    message(
      "de.results contains ANOVA-style comparison.",
      " Average Experssion will be plotted instead of logFC."
    )
    p <- ggplot(de.results, aes(`Average Intensity`, -log10(adj.P.Val), color = Class, label = Molecule)) +
      geom_point() + xlab("Average Intensity")
  } else {
    p <- ggplot(de.results, aes(logFC, -log10(adj.P.Val), color = Class, label = Molecule)) +
      geom_point() +
      geom_vline(xintercept = c(1, -1), lty = 2)
  }

  p <- p +
    geom_hline(yintercept = -log10(0.05), lty = 2) +
    facet_wrap(~contrast)
  if (show.labels) {
    sig <- de.results$adj.P.Val < 0.05
    if ("logFC" %in% colnames(de.results)) {
      sig <- sig & abs(de.results$logFC) > 1
    }
    p + geom_text(
      aes(label = ifelse(sig, Molecule, "")),
      vjust = -.5, size = 3, color = "black"
    )
  }
  .display_plot(p)
}

# colnames used in topTable
utils::globalVariables(c("logFC", "Average Intensity", "P.Value", "adj.P.Val", "contrast"))
