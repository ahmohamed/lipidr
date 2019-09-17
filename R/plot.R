#' Informative plots to investigate samples
#'
#' `lipidr` supports two types of plots for sample quality checking.\cr\cr
#' `tic` plots a bar chart for total sample intensity.\cr\cr
#' `boxplot` plots a boxplot chart to examine the distribution of values
#' per sample.
#' @param data LipidomicsExperiment object.
#' @param type plot type, either `tic` or `boxplot`. Default is `tic`.
#' @param measure Which measure to use as intensity, usually Area,
#'   Area Normalized or Height. Default is `Area`
#' @param log Whether values should be log2 transformed. Default is `TRUE`
#'
#' @return A ggplot object.
#' @export
#' @examples
#' data(data_normalized)
#'
#' plot_samples(data_normalized, type = "tic", "Area", log = TRUE)
#' plot_samples(data_normalized, type = "tic", "Background", log = FALSE)
#' plot_samples(
#'   data_normalized[, data_normalized$group == "QC"],
#'   type = "boxplot",
#'   measure = "Retention Time", log = FALSE
#' )
plot_samples <- function(data, type = c("tic", "boxplot"),
  measure = "Area", log = TRUE) {
  stopifnot(inherits(data, "LipidomicsExperiment"))
  validObject(data)
  type <- match.arg(type)
  dlong <- to_long_format(data, measure)
  measure <- sym(measure)
  if (log) {
    measure <- .check_log(data, measure)
  }
  if (type == "tic") {
    return(.display_plot(.plot_sample_tic(dlong, measure)))
  }

  .display_plot(.plot_sample_boxplot(dlong, measure))
}

.plot_sample_tic <- function(dlong, measure) {
  ggplot(dlong, aes(Sample, !!measure)) +
    stat_summary(fun.y = mean, geom = "bar") +
    facet_wrap(~filename, ncol = 1, scales = "free_y") +
    theme(axis.text.x = element_text(angle = -90, vjust = 0.5)) +
    guides(size = FALSE)
}

.plot_sample_boxplot <- function(dlong, measure) {
  ggplot(dlong, aes(Sample, !!measure)) + geom_boxplot() +
    facet_wrap(~filename, ncol = 1, scales = "free_y") +
    theme(axis.text.x = element_text(angle = -90, vjust = 0.5)) +
    guides(size = FALSE)
}

#' Informative plots to investigate lipid classes
#'
#' `lipidr` supports two types of plots for to visualize at lipid classes.\cr\cr
#' `sd` plots a bar chart for standard deviation of a certain measure in each
#' class. This plot type is usually used to look at standard deviations of
#' intensity in each class, but can also be used to look at different measures
#' such as `Retention Time`, to ensure all lipids are eluted within the expected
#' range. To assess instrumental variation apply the function to technical
#' quality control samples. \cr\cr
#' `boxplot` Plots a boxplot chart to examine the distribution of values per
#' class. This plot type is usually used to look at the intensity distribution
#' in each class, but can also be used to look at different measures, such as
#' `Retention Time` or `Background`.
#'
#' @param data LipidomicsExperiment object.
#' @param type plot type, either `boxplot` or `sd`. Default is `boxplot`.
#' @param measure Which measure to plot the distribution of: usually Area,
#'   Area Normalized, Height or Retention Time. Default is `Area`
#' @param log Whether values should be log2 transformed. Default is `TRUE`
#'   (Set FALSE for retention time).
#'
#' @return A ggplot object.
#' @export
#' @examples
#' data(data_normalized)
#'
#' d_qc <- data_normalized[, data_normalized$group == "QC"]
#' plot_lipidclass(d_qc, "sd", "Area", log = TRUE)
#' plot_lipidclass(d_qc, "sd", "Retention Time", log = FALSE)
#' plot_lipidclass(d_qc, "boxplot", "Area", log = TRUE)
#' plot_lipidclass(d_qc, "boxplot", "Retention Time", log = FALSE)
plot_lipidclass <- function(data, type = c("boxplot", "sd"),
  measure = "Area", log = TRUE) {
  stopifnot(inherits(data, "LipidomicsExperiment"))
  validObject(data)
  type <- match.arg(type)
  dlong <- to_long_format(data, measure)
  measure <- sym(measure)
  if (log) {
    measure <- .check_log(data, measure)
  }
  if (type == "sd") {
    return(.display_plot(.plot_class_sd(dlong, measure)))
  }

  .display_plot(.plot_class_boxplot(dlong, measure))
}

.plot_class_sd <- function(dlong, measure) {
  ggplot(dlong, aes(Class, !!measure, fill = Class)) +
    stat_summary(fun.y = sd, geom = "bar") +
    facet_wrap(~filename, scales = "free_x") +
    theme(axis.text.x = element_text(angle = -90, vjust = 0.5)) +
    ylab(paste("SD of", as_label(measure)))
}

.plot_class_boxplot <- function(dlong, measure) {
  ggplot(dlong, aes(Class, !!measure, fill = Class)) +
    geom_boxplot() +
    facet_wrap(~filename, scales = "free_x") +
    theme(axis.text.x = element_text(angle = -90, vjust = 0.5))
}

#' Plot logFC of lipids per class showing chain information
#'
#' Plot a chart of (log2) fold changes of lipids per class showing chain
#' lengths and saturations. If multiple molecules with the same total chain
#' length and saturation are present in the dataset, the `measure` is averaged,
#' and the number of molecules is indicated on the plot.
#'
#' @param de_results Output of [de_analysis()].
#' @param measure Which measure to plot the distribution of: logFC, P.Value,
#'   Adj.P.Val. Default is `logFC`
#' @param contrast Which comparison to plot. if not provided, defaults to the
#'   the first comparison.
#'
#' @return A ggplot object.
#' @export
#' @examples
#' data(data_normalized)
#' de_results <- de_analysis(
#'   data_normalized,
#'   HighFat_water - NormalDiet_water,
#'   measure = "Area"
#' )
#' plot_chain_distribution(de_results)
plot_chain_distribution <- function(de_results, contrast = NULL,
                                    measure = "logFC") {
  if (is.null(contrast)) {
    use_contrast <- de_results$contrast[[1]]
  }
  measure <- sym(measure)

  de_results <- de_results$Molecule %>%
    annotate_lipids() %>%
    filter(!istd) %>%
    .left_join_silent(de_results)

  de_results <- de_results %>%
    filter(contrast == use_contrast) %>%
    select(total_cl, total_cs, !!measure, Class) %>%
    mutate_at(vars(total_cl, total_cs), factor) %>%
    group_by(Class, total_cl, total_cs) %>%
    summarise(!!measure := mean(!!measure), nmolecules = n())

  p <- ggplot(de_results, aes(total_cs, total_cl, fill = logFC)) + geom_tile() +
    facet_wrap(~Class) +
    xlab("Total chain unsaturation") + ylab("Total chain length") +
    scale_fill_gradient2(midpoint = 0) +
    geom_text(aes(label = ifelse(nmolecules > 1, nmolecules, "")))

  .display_plot(p)
}

#' Plot a regulation trend line  between logFC and chain annotation
#'
#' Fit and plot a regression line of (log2) fold changes and chain
#' lengths or unsaturations. If multiple comparisons are included, one
#' regression is plotted for each.
#'
#' @param de_results Output of [de_analysis]().
#' @param annotation Whether to fit trend line against chain `length` or `unsat`.
#'
#' @return A ggplot object.
#' @export
#' @examples
#' data(data_normalized)
#' de_results <- de_analysis(
#'   data_normalized,
#'   HighFat_water - NormalDiet_water,
#'   NormalDiet_DCA - NormalDiet_water,
#'   measure = "Area"
#' )
#' plot_trend(de_results, "length")
plot_trend <- function(de_results, annotation = c("length", "unsat")) {
  annotation <- match.arg(annotation)
  x <- rlang::sym(ifelse(annotation == "length", "total_cl", "total_cs"))
  x_lab <- ifelse(annotation == "length", "Chain length", "Unsaturated bonds")
  x_breaks <- if (annotation == "length") seq(10, 80, 2) else seq(0, 12, 1)


  ggplot(de_results, aes(!!x, logFC, color = contrast)) +
    geom_hline(color = "black", yintercept = 0, lty = 2) +
    geom_smooth(alpha = 0.2) +
    labs(x = x_lab, y = "logFC") +
    scale_x_continuous(breaks = x_breaks) +
    theme(legend.position = "bottom")
}

#' Informative plots to investigate individual lipid molecules
#'
#' `lipidr` supports three types of plots for to visualize at lipid molecules.
#' \cr\cr
#' `cv` plots a bar chart for coefficient of variation of lipid molecules. This
#' plot type is usually used to investigate the CV in lipid intensity or
#' retention time, in QC samples. \cr\cr
#' `sd` plots a bar chart for standard deviations of a certain measure in each
#' lipid. This plot type is usually used to look at standard deviation of
#' intensity for each lipid, but can also be used to look at different
#' measures such as `Retention Time`, to ensure all lipids elute within
#' expected range. \cr\cr
#' `boxplot` plots a boxplot chart to examine the distribution of values per
#' lipid. This plot type is usually used to look at intensity distribution
#' for each lipid, but can also be used to look at different measures, such as
#' `Retention Time` or `Background`.
#'
#' @param data LipidomicsExperiment object.
#' @param type plot type, either `cv`, `sd` or `boxplot`. Default is `cv`.
#' @param measure Which measure to plot the distribution of: usually Area,
#'   Area Normalized or Height. Default is `Area`
#' @param log Whether values should be log2 transformed
#'   (Set FALSE for retention time). Default is `TRUE`
#'
#' @return A ggplot object.
#' @export
#' @examples
#' data(data_normalized)
#' d_qc <- data_normalized[, data_normalized$group == "QC"]
#'
#' # plot the variation in intensity and retention time of all measured
#' #   lipids in QC samples
#' plot_molecules(d_qc, "cv", "Area")
#' plot_molecules(d_qc, "cv", "Retention Time", log = FALSE)
#'
#' # plot the variation in intensity, RT of ISTD (internal standards)
#' #   in QC samples
#' d_istd_qc <- data_normalized[
#'   rowData(data_normalized)$istd,
#'   data_normalized$group == "QC"
#' ]
#' plot_molecules(d_istd_qc, "sd", "Area")
#' plot_molecules(d_istd_qc, "sd", "Retention Time", log = FALSE)
#'
#' plot_molecules(d_istd_qc, "boxplot")
#' plot_molecules(d_istd_qc, "boxplot", "Retention Time", log = FALSE)
plot_molecules <- function(data, type = c("cv", "sd", "boxplot"),
  measure = "Area", log = TRUE) {
  stopifnot(inherits(data, "LipidomicsExperiment"))
  validObject(data)
  type <- match.arg(type)
  dlong <- to_long_format(data, measure)
  measure <- sym(measure)
  if (log) {
    measure <- .check_log(data, measure)
  }
  if (type == "cv") {
    return(.display_plot(.plot_molecule_cv(dlong, measure)))
  }
  else if (type == "sd") {
    return(.display_plot(.plot_molecule_sd(dlong, measure)))
  }

  .display_plot(.plot_molecule_boxplot(dlong, measure))
}

.plot_molecule_sd <- function(dlong, measure) {
  ggplot(
    dlong,
    aes(Molecule, !!measure, fill = Class, color = Class)
  ) +
    stat_summary(fun.y = sd, geom = "bar") +
    facet_wrap(~filename, scales = "free_y") + coord_flip() +
    theme(axis.text.x = element_text(angle = -90, vjust = 0.5)) +
    ylab(paste("SD of", as_label(measure)))
}

.plot_molecule_cv <- function(dlong, measure) {
  ggplot(
    dlong,
    aes(Molecule, !!measure, fill = Class, color = Class)
  ) +
    stat_summary(fun.y = .cv, geom = "bar") + coord_flip() +
    facet_wrap(~filename, scales = "free_y") +
    theme(axis.text.x = element_text(angle = -90, vjust = 0.5)) +
    ylab(paste("CV of", as_label(measure)))
}

.plot_molecule_boxplot <- function(dlong, measure) {
  ggplot(
    dlong,
    aes(Molecule, !!measure, fill = Class, color = Class)
  ) +
    geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3) + coord_flip() +
    facet_wrap(~filename, scales = "free_y") +
    theme(axis.text.x = element_text(angle = -90, vjust = 0.5))
}

#' Plot an annotated heatmap
#' Plots a hierarchically clustered heatmap showing selected sample and
#' lipid molecule annotations.
#'
#' @param data LipidomicsExperiment object.
#' @param measure Which measure to plot the distribution of: usually Area,
#'   Area Normalized, Height or Retention Time. Default is `Area`.
#' @param molecule_annotation The column name for lipid annotation, default to `Class`.
#' @param sample_annotation The column name for sample annotation, default to `all`.
#' @param cluster_cols "none","hclust", or "k-means" for no clustering,
#'   hierarchical clustering, and k-means clustering of rows respectively.
#'   Default is "hclust".
#' @param cluster_rows "none","hclust", or "k-means" for no clustering,
#'   hierarchical clustering, and k-means clustering of rows respectively.
#'   Default is "hclust".
#' @param scale character indicating if the values should be centered and
#'   scaled in either the row direction or the column direction, or none.
#'   Corresponding values are "row", "cols" and "none".
#' @param ... Additional parameters passed to [iheatmapr::iheatmap()].
#'
#' @return A heatmap plot
#' @export
#'
#' @examples
#' data(data_normalized)
#' plot_heatmap(data_normalized, sample_annotation = "group")
plot_heatmap <- function(data, measure = "Area",
  molecule_annotation = "Class", sample_annotation = "all",
  cluster_cols = "hclust", cluster_rows = "hclust",
  scale = "rows", ...) {
  if (!requireNamespace("iheatmapr", quietly = TRUE)) {
    stop("Package 'iheatmapr' must be installed for heatmap plots")
  }
  col_annotation <- colData(data) %>% as.data.frame()
  row_annotation <- rowData(data) %>% as.data.frame()
  if (!"all" %in% sample_annotation) {
    if (is.null(sample_annotation) || FALSE %in% sample_annotation) {
      col_annotation <- NULL
    } else {
      col_annotation <- col_annotation[, sample_annotation, drop=FALSE]
    }
  }
  if (!"all" %in% molecule_annotation) {
    if (is.null(molecule_annotation) || FALSE %in% molecule_annotation) {
      row_annotation <- NULL
    } else {
      row_annotation <- row_annotation[, molecule_annotation, drop=FALSE]
    }
  }
  dim_names <- metadata(data)$dimnames
  iheatmapr::iheatmap(assay(data, measure),
    row_title = dim_names[[1]], col_title = dim_names[[2]],
    cluster_cols = cluster_cols, cluster_rows = cluster_rows,
    row_annotation = row_annotation, col_annotation = col_annotation,
    scale = scale, ...
  )
}

.display_plot <- function(p) {
  if (.myDataEnv$interactive) {
    p <- plotly::ggplotly(p)
  }
  return(p)
}



# colnames used in plot_chain_distribution
utils::globalVariables(c("nmolecules"))
