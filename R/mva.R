# colnames used internally in top.lipids
utils::globalVariables(c("molrank"))


#' Perform multivariate analyses to investigate sample clustering
#'
#' `mva` performs multivariate analysis using several possible methods.
#' The available methods are PCA, PCoA, OPLS and OPLS-DA. The OPLS method
#' requires a numeric y-variable, whilst OPLS-DA requires two groups for
#' comparison. By default, for OPLS and OPLS-DA the number of predictive and
#' orthogonal components are set to 1.
#' Blank samples are automatically detected (using TIC) and excluded.
#' Missing data are imputed using average lipid intensity across all samples.
#'
#' @param data SkylineExperiment object created by [read_skyline()].
#' @param measure Which measure to use as intensity, usually Area (default).
#'   The measure should be already summarized and normalized.
#' @param method Either PCA, PCoA, OPLS or OPLS-DA.  Default is `PCA`.
#' @param group_col Sample annotation to use as grouping column. If not
#'   provided, samples are treated independently.
#' @param groups A numeric grouping (OPLS) or two groups to be used for
#'   supervised analysis (OPLS-DA), ignored in other methods.
#' @param ... Extra arguments to be passed to [opls()] for OPLS-DA,
#'   ignored in other methods.
#'
#' @return Multivariate analysis results in `mvaresults` object.
#'   The object contains the following:\itemize{
#'     \item scores    Sample scores
#'     \item loadings   Feature or component loadings (not for PCoA)
#'     \item method   Multivariate method that was used
#'     \item row_data   Lipid molecule annotations
#'     \item col_data   Sample annotations
#'     \item original_object   Original output object as returned by
#'         corresponding analysis methods
#'   }
#'
#' @importFrom stats cmdscale prcomp
#' @importFrom ropls opls
#' @importFrom forcats fct_drop
#'
#' @export
mva <- function(data, measure = "Area",
                method = c("PCA", "PCoA", "OPLS", "OPLS-DA"),
                group_col = NULL, groups = NULL, ...) {
  stopifnot(inherits(data, "SkylineExperiment"))
  method <- match.arg(method)
  data_f <- data[!rowData(data)$itsd, !.is_blank(data)]
  d <- data_f %>%
    assay(measure) %>%
    .replace_na_rowmean() %>%
    t()

  if (is.null(group_col)) {
    if (ncol(colData(data)) > 0) {
      group_col <- names(colData(data))[[1]]
    } else {
      group_col <- NULL
    }
  }

  if (method == "PCA") {
    object <- prcomp(d, scale = TRUE)
    return(structure(list(
      scores = object$x,
      loadings = object$rotation,
      variance = object$stdev,
      var.pct = round((object$sdev^2) / sum(object$sdev^2), 3) * 100,
      method = "PCA",
      row_data = rowData(data_f),
      col_data = colData(data_f),
      group_col = group_col
    ),
    class = c("mvaResults", "pca"),
    original_object = object
    ))
  }
  else if (method == "PCoA") {
    return(structure(list(
      scores = cmdscale(dist(d)),
      method = "PCoA",
      row_data = rowData(data_f),
      col_data = colData(data_f),
      group_col = group_col
    ),
    class = c("mvaResults", "pcoa")
    ))
  }

  # method either "OPLS" or "OPLS-DA"
  if (is.null(group_col)) {
    stop("Please add clinical data or specify a group column")
  }

  if (length(group_col) == 1) {
    # group_col is the name of the grouping column
    group_vector <- colData(data_f)[[group_col]]
  } else {
    # group_col is a vector with grouping
    group_vector <- group_col
  }

  if (method == "OPLS-DA") {
    # group_vector should either have 2 values,
    # or a subset should be provided in groups
    if (length(unique(group_vector)) != 2) {
      if (length(unique(groups)) != 2) {
        stop("Please provide 2 groups for comparison in OPLS-DA")
      }
      # all groups in the subset should be present in the vector
      if (!all(groups %in% group_vector)) {
        stop("Provided groups are not in the grouping column.")
      }
    }
  } else {
    # This is OPLS
    # group_vector should be numeric
    if (!is.numeric(group_vector)) {
      stop("Please provide a numeric y-variable for comparison in OPLS")
    }
    # if group subset is provided, it should be numeric,
    # and values should be present in group_vector
    if (!is.null(groups)) {
      if (!is.numeric(groups)) {
        stop(
          "Please provide a numeric groups variable for comparison in OPLS"
        )
      }
    }
    if (!all(groups %in% group_vector)) {
      stop("Provided groups are not in the grouping column.")
    }
  }
  # By now we know that group_vector is of correct type
  # groups, if provided, have correct type and values.
  if (!is.null(groups)) {
    data_f <- data_f[, group_vector %in% groups]
    d <- d[group_vector %in% groups, ]
    group_vector <- fct_drop(group_vector[group_vector %in% groups])
  }
  object <- run_opls(d, y = group_vector, ...)

  return(structure(list(
    scores = data.frame(object@scoreMN[, 1], object@orthoScoreMN[, 1]),
    loadings = data.frame(object@loadingMN[, 1], object@orthoLoadingMN[, 1]),
    summary = object@modelDF,
    method = method,
    row_data = rowData(data_f),
    col_data = colData(data_f),
    group_col = group_col
  ),
  class = c("mvaResults", "opls"),
  original_object = object
  ))

}

run_opls <- function(data, y,
                     predI = 1, orthoI = 1,
                     scaleC = "standard",
                     fig.pdfC = NULL, ...) {
  opls(
    data,
    y = y,
    predI = predI, orthoI = orthoI,
    scaleC = scaleC,
    fig.pdfC = fig.pdfC, ...
  )
}

#' @importFrom stats var qf median dist
plot_opls <- function(mvaresults, components,
                      color_by, ellipse = TRUE, hotelling = TRUE) {
  ret <- .get_mds_matrix(mvaresults, components, color_by)
  d <- ret$mds_matrix
  color_by <- ret$color_by

  N <- nrow(d)
  pscores <- d[, 2]
  oscores <- d[, 3]

  hotFisN <- (N - 1) * 2 * (N^2 - 1) / (N^2 * (N - 2)) * qf(0.95, 2, N - 2)


  p <- ggplot(d, aes_string(
    colnames(d)[[2]], colnames(d)[[3]],
    label = "Sample", color = color_by
  ))

  if (ellipse) {
    p <- p + stat_ellipse(
      geom = "polygon", alpha = 0.3, linetype = "blank",
      aes_string(fill = color_by), type = "norm"
    )
  }
  if (hotelling) {
    p <- p + gg_circle(
      rx = sqrt(var(pscores) * hotFisN),
      ry = sqrt(var(oscores) * hotFisN),
      xc = 0, yc = 0
    )
  }
  sm <- mvaresults$summary
  p <- p + geom_hline(yintercept = 0, color = "gray") +
    geom_vline(xintercept = 0, color = "gray") +
    geom_point(size = 3) +
    xlab(paste("p1:", sm["p1", "R2X"] * 100, "%")) +
    ylab(paste("o1:", sm["o1", "R2X"] * 100, "%")) +
    labs(color = "Group", fill = "Group") +
    theme_grey(base_size = 10) +
    annotate(
      "text",
      x = Inf, y = Inf,
      label = paste("R2X:", sm["sum", "R2X(cum)"]),
      vjust = 1, hjust = 1, size = 3
    ) +
    annotate(
      "text",
      x = Inf, y = Inf,
      label = paste("R2Y:", sm["sum", "R2Y(cum)"]),
      vjust = 2.5, hjust = 1, size = 3
    ) +
    annotate(
      "text",
      x = Inf, y = Inf,
      label = paste("Q2:", sm["sum", "Q2(cum)"]),
      vjust = 4, hjust = 1, size = 3
    )

  .display_plot(p)
}

#' @describeIn mva plots a multivariate scatterplot of sample scores to investigate
#' sample clustering.
#'
#' @param mvaresults Results obtained from [mva()].
#' @param components Which components to plot. Ignored for PCoA, OPLS and
#'   OPLS-DA results. Default is first 2 components.
#' @param color_by Sample annotation (or lipid annotation in case of
#'   `plot_mva_loadings`) to use as color. Defaults to individual samples /
#'   lipids
#'
#' @return `plot_mva` returns a ggplot of the sample scores.
#' @export
#' @examples
#' data(data_normalized)
#'
#' # PCA
#' mvaresults <- mva(data_normalized, measure = "Area", method = "PCA")
#' plot_mva(mvaresults, color_by = "group")
#' # NOT RUN
#' # plot_mva(mvaresults, color_by = "Diet", components = c(2, 3))
#'
#' # PCoA
#' mvaresults <- mva(data_normalized, measure = "Area", method = "PCoA")
#' # NOT RUN
#' # plot_mva(mvaresults, color_by = "group")
#'
#' # OPLS-DA
#' mvaresults <- mva(
#'   data_normalized,
#'   method = "OPLS-DA", group_col = "Diet", groups=c("HighFat", "Normal")
#' )
#' plot_mva(mvaresults, color_by = "group")
plot_mva <- function(mvaresults, components = c(1, 2), color_by = NULL) {
  stopifnot(inherits(mvaresults, "mvaResults"))
  if (inherits(mvaresults, "opls")) {
    return(plot_opls(mvaresults, components, color_by))
  }

  ret <- .get_mds_matrix(mvaresults, components, color_by)
  mds_matrix <- ret$mds_matrix
  color_by <- ret$color_by
  cols <- colnames(mds_matrix)

  p <- ggplot(mds_matrix, aes_string(
    cols[[2]], cols[[3]],
    label = "Sample", color = color_by
  )) + geom_point(size = 3, pch = 16) +
    geom_text(vjust = -.5, size = 3, color = "black")

  if (inherits(mvaresults, "pca") && !is.null(mvaresults$var.pct)) {
    xlabel <- paste(
      cols[[2]],
      "(", mvaresults$var.pct[[ components[[1]] ]], "%)"
    )
    ylabel <- paste(
      cols[[3]],
      "(", mvaresults$var.pct[[ components[[2]] ]], "%)"
    )
    p <- p + xlab(xlabel) + ylab(ylabel)
  }
  .display_plot(p)
}


#' @describeIn mva Plot a multivariate scatterplot of feature loadings
#' to investigate feature importance.
#'
#' @param top.n Number of top ranked features to highlight in the plot.
#'   If omitted, returns top 10 lipids.
#'
#' @return `plot_mva_loadings` returns a ggplot of the loadings.
#' @export
#'
#' @examples
#' plot_mva_loadings(mvaresults, color_by = "Class", top.n = 10)
plot_mva_loadings <- function(mvaresults, components = c(1, 2),
                              color_by = NULL,
                              top.n = nrow(mvaresults$loadings)) {
  stopifnot(inherits(mvaresults, "mvaResults"))
  stopifnot(inherits(mvaresults, "opls"))
  ret <- .get_loading_matrix(mvaresults, components, color_by)
  mds_matrix <- ret$mds_matrix
  mds_matrix$molrank <- rank(-abs(mds_matrix[, 2]))

  color_by <- ret$color_by
  sm <- mvaresults$summary

  p <- ggplot(mds_matrix, aes_string(
    colnames(mds_matrix)[[2]], colnames(mds_matrix)[[3]],
    color = color_by
  )) +
    geom_hline(yintercept = 0, color = "gray") +
    geom_vline(xintercept = 0, color = "gray") +
    xlab(paste("p1:", sm["p1", "R2X"] * 100, "%")) +
    ylab(paste("o1:", sm["o1", "R2X"] * 100, "%")) +
    theme_grey(base_size = 10) +
    geom_point(size = 3, pch = 16, aes(alpha = molrank > top.n)) +
    scale_alpha_manual(values = c(1, 0.5))

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    xlimits <- max(abs(mds_matrix[, 2])) * 1.25
    ylimits <- max(abs(mds_matrix[, 3])) * 1.1

    p <- p + ggrepel::geom_label_repel(
      aes(label = ifelse(molrank > top.n, "", Molecule)),
      size = 2.4, direction = "both", segment.alpha = 0.6,
      label.padding = 0.15, force = 0.5
    ) + xlim(-xlimits, xlimits) + ylim(-ylimits, ylimits)
  } else {
    p <- p + geom_text(
      vjust = -.5, size = 3, color = "black",
      aes(label = ifelse(molrank > top.n, "", Molecule))
    )
  }

  .display_plot(p)
}

#' @describeIn mva extracts top lipids from OPLS-DA results
#'
#' @return `top_lipids` returns s dataframe of `top.n` lipids with
#'   their annotations.
#' @export
#'
#' @examples
#' top_lipids(mvaresults, top.n = 10)
top_lipids <- function(mvaresults, top.n = 10) {
  stopifnot(inherits(mvaresults, "mvaResults"))
  stopifnot(inherits(mvaresults, "opls"))

  ret <- .get_loading_matrix(mvaresults, c(1, 2), "Molecule")
  mds_matrix <- ret$mds_matrix
  mds_matrix$molrank <- rank(-abs(mds_matrix[, 2]))
  mds_matrix <- mds_matrix[, -c(1, 2, 3)]
  mds_matrix %>%
    filter(molrank <= top.n) %>%
    arrange(molrank)
}

# Function to plot Hotelling's T-squared ellipse
# Adapted from https://github.com/tyrannomark/bldR/blob/master/R/L2017.R
# GPL-3 license
gg_circle <- function(rx, ry, xc, yc, color = "black", fill = NA, ...) {
  x <- xc + rx * cos(seq(0, pi, length.out = 100))
  ymax <- yc + ry * sin(seq(0, pi, length.out = 100))
  ymin <- yc + ry * sin(seq(0, -pi, length.out = 100))
  annotate(
    "ribbon",
    x = x, ymin = ymin, ymax = ymax,
    color = color, fill = fill, ...
  )
}

.process_prcomp <- function(x, choices = c(1, 2), scale = 1) {
  if (length(choices) != 2L) {
    stop("length of choices must be 2")
  }
  scores <- x$x
  if (!length(scores)) {
    stop(gettextf("object '%s' has no scores", deparse(substitute(x))),
      domain = NA
    )
  }

  lam <- x$sdev[choices]
  n <- NROW(scores)
  lam <- lam * sqrt(n)
  if (scale < 0 || scale > 1) {
    warning("'scale' is outside [0, 1]")
  }
  if (scale != 0) {
    lam <- lam^scale
  } else {
    lam <- 1
  }

  list(t(t(scores[, choices]) / lam), t(t(x$rotation[, choices]) * lam))
  return(lam)
}

.get_loading_matrix <- function(mvaresults, components = c(1, 2),
                                color_by = NULL) {
  stopifnot(inherits(mvaresults, "mvaResults"))
  mds_matrix <- mvaresults$loadings[, components]
  mds_matrix <- mds_matrix %>%
    as.data.frame() %>%
    rownames_to_column("LipidID")

  if (!is.null(color_by)) {
    row_data <- mvaresults$row_data %>%
      as.data.frame() %>%
      rownames_to_column("LipidID")

    mds_matrix <- mds_matrix %>%
      .left_join_silent(row_data)
  }
  return(list(mds_matrix = mds_matrix, color_by = color_by))
}

.get_mds_matrix <- function(mvaresults, components = c(1, 2), color_by = NULL) {
  stopifnot(inherits(mvaresults, "mvaResults"))
  mds_matrix <- mvaresults$scores
  if (mvaresults$method == "PCA") {
    mds_matrix <- mds_matrix[, components]
  }

  mds_matrix <- mds_matrix %>%
    as.data.frame() %>%
    rownames_to_column("Sample")

  if (is.null(color_by)) {
    if (!is.null(mvaresults$group_col)) {
      color_by <- mvaresults$group_col
    } else {
      color_by <- "Sample"
    }
  }
  if (color_by != "Sample") {
    col_data <- mvaresults$col_data %>%
      as.data.frame() %>%
      rownames_to_column("Sample")

    mds_matrix <- mds_matrix %>%
      .left_join_silent(col_data)
  }
  return(list(mds_matrix = mds_matrix, color_by = color_by))
}
