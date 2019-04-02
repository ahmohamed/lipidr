#' @import ggplot2
{}


#' Plot a bar chart for total sample intensity
#'
#' @param data Skyline data.frame created by [read_skyline()].
#' @param measure Which measure to use as intensity, usually Area,
#'   Area.Normalized or Height.
#' @param log Whether values should be log2 transformed.
#'
#' @return A ggplot object.
#' @export
#' @examples
#' data(data_normalized)
#'
#' plot_sample_tic(data_normalized, "Area", log = TRUE)
#' plot_sample_tic(data_normalized, "Background", log = FALSE)
plot_sample_tic <- function(data, measure = "Area", log = TRUE) {
  stopifnot(inherits(data, "SkylineExperiment"))
  dlong <- to_long_format(data, measure)
  if (log == TRUE) {
    measure <- .check_log(data, measure)
  }
  p <- ggplot(dlong, aes_string("Sample", measure)) + stat_sum(geom = "bar") +
    facet_wrap(~filename, ncol = 1, scales = "free_y") +
    theme(axis.text.x = element_text(angle = -90, vjust = 0.5)) +
    guides(size = FALSE)

  .display_plot(p)
}

#' Plot a boxplot chart to examine the distribution of values per sample
#'
#' The function should usually be used to look at intensity distribution in
#' each sample ensuring they are normalized. It can also be used to look at
#' different measures such as  `Retention.Time` or `Background`.
#'
#' @param data Skyline data.frame created by [read_skyline()].
#' @param measure Which measure to plot the distribution of: usually Area,
#'   Area.Normalized or Height.
#' @param log Whether values should be log2 transformed.
#'
#' @return A ggplot object.
#' @export
#' @examples
#' data(data_normalized)
#'
#' plot_sample_boxplot(
#'   data_normalized[, data_normalized$group == "QC"],
#'   measure = "Retention.Time", log = FALSE
#' )
#' # NOT RUN
#' # plot_sample_boxplot(data_normalized, "Area", log = TRUE)
plot_sample_boxplot <- function(data, measure = "Area", log = TRUE) {
  stopifnot(inherits(data, "SkylineExperiment"))
  dlong <- to_long_format(data, measure)
  if (log == TRUE) {
    measure <- .check_log(data, measure)
  }

  p <- ggplot(dlong, aes_string("Sample", measure)) + geom_boxplot() +
    facet_wrap(~filename, ncol = 1, scales = "free_y") +
    theme(axis.text.x = element_text(angle = -90, vjust = 0.5)) +
    guides(size = FALSE)

  .display_plot(p)
}



#' Plot a bar chart for standard deviation of a certain measure in each class
#'
#' The function is usually used to look at standard deviations of intensity in
#' each class, but can also be used to look at different measures such as
#' `Retention.Time`, to ensure all lipids are eluted within the expected range.
#' To assess instrumental variation apply the function to technical quality
#' control samples.
#'
#' @param data Skyline data.frame created by [read_skyline()].
#' @param measure Which measure to plot the distribution of: usually Area,
#'   Area.Normalized, Height or Retention.Time
#' @param log Whether values should be log2 transformed
#'   (Set FALSE for retention time).
#'
#' @return A ggplot object.
#' @export
#' @examples
#' data(data_normalized)
#'
#' plot_class_sd(data_normalized, "Area", log = TRUE)
#' plot_class_sd(data_normalized, "Retention.Time", log = FALSE)
plot_class_sd <- function(data, measure = "Area", log = TRUE) {
  stopifnot(inherits(data, "SkylineExperiment"))
  dlong <- to_long_format(data, measure)
  if (log == TRUE) {
    measure <- .check_log(data, measure)
  }

  p <- ggplot(dlong, aes_string("Class", measure, fill = "Class")) +
    stat_summary(fun.y = sd, geom = "bar") +
    facet_wrap(~filename, scales = "free_x") +
    theme(axis.text.x = element_text(angle = -90, vjust = 0.5)) +
    ylab(paste("SD of", measure))

  .display_plot(p)
}

#' Plot a boxplot chart to examine the distribution of values per class
#'
#' The function is usually used to look at the intensity distribution in each
#' class, but can also be used to look at different measures, such as
#' `Retention.Time` or `Background`.
#'
#' @param data Skyline data.frame created by [read_skyline()].
#' @param measure Which measure to plot the distribution of: usually Area,
#'   Area.Normalized or Height.
#' @param log Whether values should be log2 transformed.
#'
#' @return A ggplot object.
#' @export
#' @examples
#' data(data_normalized)
#' d_qc <- data_normalized[, data_normalized$group == "QC"]
#'
#' plot_class_boxplot(d_qc, "Area", log = TRUE)
#' plot_class_boxplot(d_qc, "Retention.Time", log = FALSE)
plot_class_boxplot <- function(data, measure = "Area", log = TRUE) {
  stopifnot(inherits(data, "SkylineExperiment"))
  dlong <- to_long_format(data, measure)
  if (log == TRUE) {
    measure <- .check_log(data, measure)
  }

  p <- ggplot(dlong, aes_string("Class", measure, fill = "Class")) +
    geom_boxplot() +
    facet_wrap(~filename, scales = "free_x") +
    theme(axis.text.x = element_text(angle = -90, vjust = 0.5))

  .display_plot(p)
}


#' Plot logFC of lipids per class showing chain information
#' Plot a chart of (log2) fold changes of lipids per class showing chain
#' lengths and saturations. If multiple molecules with the same total chain
#' length and saturation are present in the dataset, the `measure` is averaged,
#' and the number of molecules is indicated on the plot.
#'
#' @param de_results Output of [de_analysis()].
#' @param measure Which measure to plot the distribution of: logFC, P.Value,
#'   Adj.P.Val.
#' @param contrast Which comparison to plot.
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
    filter(!itsd) %>%
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

#' Plot a bar chart for standard deviations of a certain measure in each lipid
#'
#' The function is usually used to look at standard deviation of intensity for
#' each lipid, but can also be used to look at different measures such as
#' `Retention.Time`, to ensure all lipids elute within expected range.
#'
#' @param data Skyline data.frame created by [read_skyline()].
#' @param measure Which measure to plot the distribution of: usually Area,
#'   Area.Normalized or Height.
#' @param log Whether values should be log2 transformed
#'   (Set FALSE for retention time).
#'
#' @return A ggplot object.
#' @export
#' @examples
#' data(data_normalized)
#'
#' # plot the variation in intensity of ITSD (internal standards) in QC samples
#' d_itsd_qc <- data_normalized[
#'   rowData(data_normalized)$itsd,
#'   data_normalized$group == "QC"
#' ]
#' plot_molecule_sd(d_itsd_qc, "Area")
#'
#' # plot the variation in intensity and retention time of all measured
#' #   lipids in QC samples
#' plot_molecule_sd(data_normalized[, data_normalized$group == "QC"], "Area")
#' plot_molecule_sd(
#'   data_normalized[, data_normalized$group == "QC"],
#'   measure = "Retention.Time", log = FALSE
#' )
plot_molecule_sd <- function(data, measure = "Area", log = TRUE) {
  stopifnot(inherits(data, "SkylineExperiment"))
  dlong <- to_long_format(data, measure)
  if (log == TRUE) {
    measure <- .check_log(data, measure)
  }
  p <- ggplot(
    dlong,
    aes_string("Molecule", measure, fill = "Class", color = "Class")
  ) +
    stat_summary(fun.y = sd, geom = "bar") +
    facet_wrap(~filename, scales = "free_y") + coord_flip() +
    theme(axis.text.x = element_text(angle = -90, vjust = 0.5)) +
    ylab(paste("SD of", measure))

  .display_plot(p)
}

#' Bar chart for coefficient of variation of lipid molecules
#'
#' Plot a bar chart for coefficient of variation (CV) of a certain measure
#' in each lipid. The functions is usually used to investigate the CV in lipid
#' intensity or retention time, in QC samples.
#'
#' @param data Skyline data.frame created by [read_skyline()].
#' @param measure Which measure to plot the distribution of: usually Area,
#'   Area.Normalized or Height.
#' @param log whether values should be log2 transformed
#'   (Set FALSE for retention time).
#'
#' @return A ggplot object.
#' @export
#' @examples
#' data(data_normalized)
#'
#' # plot the variation in intensity and retention time of all measured
#' #   lipids in QC samples
#' d_qc <- data_normalized[, data_normalized$group == "QC"]
#' plot_molecule_cv(d_qc, "Area")
#' plot_molecule_cv(d_qc, "Retention.Time", log = FALSE)
plot_molecule_cv <- function(data, measure = "Area", log = TRUE) {
  stopifnot(inherits(data, "SkylineExperiment"))
  dlong <- to_long_format(data, measure)
  if (log == TRUE) {
    measure <- .check_log(data, measure)
  }
  p <- ggplot(
    dlong,
    aes_string("Molecule", measure, fill = "Class", color = "Class")
  ) +
    stat_summary(fun.y = .cv, geom = "bar") + coord_flip() +
    facet_wrap(~filename, scales = "free_y") +
    theme(axis.text.x = element_text(angle = -90, vjust = 0.5)) +
    ylab(paste("SD of", measure))

  .display_plot(p)
}

#' Plot a boxplot chart to examine the distribution of values per lipid
#'
#' The function is usually used to look at intensity distribution for each
#' lipid, but can also be used to look at different measures, such as
#' `Retention.Time` or `Background`.
#'
#' @param data Skyline data.frame created by [read_skyline()].
#' @param measure Which measure to plot the distribution of: usually Area,
#'   Area.Normalized or Height.
#' @param log Whether values should be log2 transformed.
#'
#' @return A ggplot object.
#' @export
#' @examples
#' data(data_normalized)
#'
#' # plot the variation in intensity of ITSD (internal standards) in QC samples
#' d_itsd_qc <- data_normalized[
#'   rowData(data_normalized)$itsd,
#'   data_normalized$group == "QC"
#' ]
#' plot_molecule_boxplot(d_itsd_qc)
#' plot_molecule_boxplot(d_itsd_qc, "Retention.Time", log = FALSE)
plot_molecule_boxplot <- function(data, measure = "Area", log = TRUE) {
  stopifnot(inherits(data, "SkylineExperiment"))
  dlong <- to_long_format(data, measure)
  if (log == TRUE) {
    measure <- .check_log(data, measure)
  }
  p <- ggplot(
    dlong,
    aes_string("Molecule", measure, fill = "Class", color = "Class")
  ) +
    geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3) + coord_flip() +
    facet_wrap(~filename, scales = "free_y") +
    theme(axis.text.x = element_text(angle = -90, vjust = 0.5))

  .display_plot(p)
}

.display_plot <- function(p) {
  if (.myDataEnv$interactive) {
    p <- plotly::ggplotly(p)
  }
  return(p)
}
