#' Lipid set enrichment analysis (LSEA)
#'
#' @param de.results Output of [de_analysis()].
#' @param rank.by Statistic used to rank the lipid list.  Default is `logFC`.
#' @param ... Extra parameters passed to [fgsea::fgsea()].
#'
#' @return `lsea` returns enrichment results (data.frame) as returned from
#'   [fgsea::fgsea()].
#'   The results also contain the following attributes: \itemize{
#'     \item de.results    Original de.results input.
#'     \item rank.by   Measure used to rank lipid molecules.
#'     \item sets   Lipid sets tested, with their member molecules.
#'   }
#'
#' @importFrom rlang :=
#' @export
#'
#' @examples
#' data(data_normalized)
#' de_results <- de_analysis(
#'   data_normalized,
#'   HighFat_water - NormalDiet_water,
#'   measure = "Area"
#' )
#' enrich_results <- lsea(
#'   de_results,
#'   rank.by = "logFC", minSize = 4, nperm = 1000
#' )
lsea <- function(de.results,
                 rank.by = c("logFC", "P.Value", "Adj.P.Val"), ...) {
  rank.by <- match.arg(rank.by)
  rank_by_sym <- sym(rank.by)

  sets <- gen_lipidsets(de.results$Molecule)

  de.results <- de.results %>%
    arrange(-!!rank_by_sym) %>%
    group_by(contrast, Molecule) %>%
    summarise(!!rank_by_sym := first(!!rank_by_sym))

  contrast.list <- split.data.frame(de.results, de.results$contrast)


  res <- lapply(contrast.list, function(x)
    .gsea_fun(
      pathways = sets,
      stats = to_named_list(x, "Molecule", rank.by),
      ...
  )) %>%
    bind_rows(.id = "contrast")

  res <- res %>% arrange(padj) %>% rename(set = pathway)
  attr(res, "de.results") <- de.results
  attr(res, "rank.by") <- rank.by
  attr(res, "sets") <- sets

  return(res)
}

.gsea_fun <- function(pathways, stats, minSize=2, nperm=10000, ...) {
  fgsea::fgsea(
    pathways = pathways, stats = stats,
    minSize = minSize, nperm = nperm, ...
  )
}

#' @describeIn lsea gets a list of significantly changed lipid sets
#'
#' @param enrich.results Output of [lsea()].
#' @param p.cutoff Significance threshold.  Default is `0.05`.
#' @param size.cutoff Minimum number of lipids in a set tested for enrichment.
#'    Default is `2`.
#'
#' @return `significant_lipidsets` returns a list of character vectors of
#'   significantly enriched sets for each contrast.
#' @export
#' @examples
#' sig_lipidsets <- significant_lipidsets(enrich_results)
significant_lipidsets <- function(enrich.results, p.cutoff = 0.05,
                                  size.cutoff = 2) {
  enrich.results %>%
    filter(padj < p.cutoff, size > size.cutoff) %>%
    (function(x) split(x$set, x$contrast))
}


#' @describeIn lsea is usually used to look at log2 fold change
#' distribution of lipids in each class, marking significantly enriched classes.
#' Can also be used to plot `P.Value` or `Adj.P.Val`.
#'
#' @param significant.sets List of significantly changed lipid sets
#'   (output of [significant_lipidsets()]).
#' @param measure Which measure to plot the distribution of: logFC, P.Value,
#'   Adj.P.Val. Default is `logFC`.
#'
#' @return `plot_class_enrichment` returns a ggplot object.
#' @export
#' @examples
#' plot_class_enrichment(de_results, sig_lipidsets)
plot_class_enrichment <- function(de.results, significant.sets,
                                  measure = "logFC") {
  significant.sets <- lapply(
    significant.sets,
    function(c) sub("^Class_", "", c[grep("^Class_", c)])
  )
  de_results <- de.results$Molecule %>%
    annotate_lipids() %>%
    .left_join_silent(de_results) %>%
    group_by(contrast) %>%
    mutate(Significant = Class %in% significant.sets[[ contrast[[1]] ]]) %>%
    ungroup()

  p <- ggplot(de_results, aes_string("Class", measure, color = "Significant")) +
    geom_boxplot() + geom_hline(yintercept = 0, lty = 2) +
    facet_wrap(~contrast, scales = "free_x") +
    scale_color_manual(values = c("black", "red")) +
    theme(axis.text.x = element_text(angle = -90, vjust = 0.5))

  .display_plot(p)
}

#' Generate lipid sets from lipid molecule names
#'
#' @param molecules A character vector containing lipid molecule names.
#'
#' @return List of lipid sets
#' @export
#' @importFrom tidyr gather unite
#' @examples
#' data(data_normalized)
#' molecules <- rowData(data_normalized)$Molecule
#' gen_lipidsets(molecules)
gen_lipidsets <- function(molecules) {
  data_ <- annotate_lipids(molecules) %>%
    filter(!itsd) %>%
    distinct(Molecule, Class, total_cl, total_cs) %>%
    gather("collection", "value", Class, total_cl, total_cs) %>%
    unite("set", collection, value, sep = "_")

  split(as.character(data_$Molecule), data_$set)
}

# colnames used internally in gen.lipidset
utils::globalVariables(c("collection", "value"))

# colnames output by fgsea
utils::globalVariables(c("padj", "size", "pathway"))
