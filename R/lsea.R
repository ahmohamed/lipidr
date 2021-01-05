#' Lipid set enrichment analysis (LSEA)
#'
#' @param de.results Output of [de_analysis()].
#' @param rank.by Statistic used to rank the lipid list.  Default is `logFC`.
#' @param min_size Minimum number of molecules in a set to be included
#'   in enrichment.
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
#'   rank.by = "logFC", min_size = 4, nperm = 1000
#' )
lsea <- function(de.results,
  rank.by = c("logFC", "P.Value", "adj.P.Val"), min_size = 2, ...) {
  rank.by <- match.arg(rank.by)
  rank_by_sym <- sym(rank.by)

  sets <- gen_lipidsets(de.results, min_size)
  if (length(unlist(sets)) == 0) {
    stop('Unable to generate lipid sets, possibly because of missing annotations.')
  }

  de.results <- de.results %>%
    arrange(-!!rank_by_sym) %>%
    group_by(contrast, Molecule) %>%
    summarise(!!rank_by_sym := first(!!rank_by_sym))

  contrast.list <- split.data.frame(de.results, de.results$contrast)


  res <- lapply(contrast.list, function(x)
    .gsea_fun(
      pathways = sets,
      stats = to_named_list(x, "Molecule", rank.by),
      minSize = min_size,
      ...
    )) %>%
    bind_rows(.id = "contrast")

  res <- res %>% arrange(padj) %>% rename(set = pathway)
  attr(res, "de.results") <- de.results
  attr(res, "rank.by") <- rank.by
  attr(res, "sets") <- sets

  return(res)
}

.gsea_fun <- function(pathways, stats, minSize = 2, nperm = 10000, ...) {
  fgsea::fgsea(
    pathways = pathways, stats = stats,
    minSize = minSize, ...
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


#' @export
#' @rdname lsea
plot_class_enrichment <- function(de.results, significant.sets,
  measure = "logFC") {
  .Deprecated("plot_enrichment")
  plot_enrichment(de.results, significant.sets, measure=measure)
}

#' @describeIn lsea is usually used to look at log2 fold change
#' distribution of lipids in each class, chain length or unsaturation,
#' marking significantly enriched sets. It can also be used to plot `P.Value`
#' or `Adj.P.Val`.
#'
#' @param significant.sets List of significantly changed lipid sets
#'   (output of [significant_lipidsets()]).
#' @param measure Which measure to plot the distribution of: logFC, P.Value,
#'   Adj.P.Val. Default is `logFC`.
#' @param annotation Which lipid set collection to plot.
#'
#' @return `plot_enrichment` returns a ggplot object.
#' @importFrom forcats fct_recode
#' @export
#' @examples
#' plot_enrichment(de_results, sig_lipidsets, annotation="class")
#' plot_enrichment(de_results, sig_lipidsets, annotation="length")
plot_enrichment <- function(de.results, significant.sets,
  annotation=c("class", "length", "unsat"), measure = "logFC") {

  annotation <- match.arg(annotation)
  collection <- c(class="Class", length="total_cl", unsat="total_cs")[[annotation]]
  prefix = paste0("^", collection, "_")

  significant.sets <- lapply(
    significant.sets,
    function(c) sub(prefix, "", c[grep(prefix, c)])
  )
  de_results <- de.results %>% group_by(contrast) %>%
    mutate(Significant = factor(
      as.character(!!sym(collection)) %in% significant.sets[[ contrast[[1]] ]],
      levels = c(TRUE, FALSE)
    )) %>%
    mutate(Enrichment = fct_recode(
      Significant, "Significant"="TRUE", "Not significant"="FALSE")
    ) %>%
    ungroup() %>%
    mutate(
      collection = as.character(!!sym(collection)),
      measure = !!sym(measure),
    )

  p <- ggplot(de_results, aes(collection, measure, color = Enrichment)) +
    geom_boxplot() + geom_hline(yintercept = 0, lty = 2) +
    facet_wrap(~contrast, scales = "free_x") +
    scale_color_manual(values = c(`Not significant`="black", `Significant`="red"), drop=FALSE) +
    theme(axis.text.x = element_text(angle = -90, vjust = 0.5))

  .display_plot(p)
}

#' Generate lipid sets from lipid molecule names
#'
#' @param molecules A character vector containing lipid molecule names.
#' @param min_size Minimum number of molecules in a set to be included
#'   in enrichment.
#'
#' @return List of lipid sets
#' @export
#' @importFrom tidyr gather unite
#' @examples
#' data(data_normalized)
#' molecules <- rowData(data_normalized)$Molecule
#' gen_lipidsets(molecules)
gen_lipidsets <- function(molecules, min_size=2) {
  if (is.data.frame(molecules) &&
    all(c( "Molecule", "Class", "total_cl", "total_cs") %in% colnames(molecules))
    ) {
    df <- molecules
  } else {
    df <- annotate_lipids(molecules)
  }
  clean_df <- df %>%
    filter(!istd) %>%
    distinct(Molecule, Class, total_cl, total_cs)

  data_  <- clean_df%>%
    gather("collection", "value", Class, total_cl, total_cs) %>%
    filter(!is.na(value)) %>% # Remove NA sets
    unite("set", collection, value, sep = "_") %>%
    group_by(set) %>% filter(n() >= min_size)

  sets <- split(as.character(data_$Molecule), data_$set)
  lens <- lapply(sets, length) == length(unique(unlist(sets)))
  if (any(lens)) {
    warning(
      'These sets contained all molecules, and were excluded: ',
      paste(names(sets[lens]), collapse = ", ")
    )
    sets <- sets[!lens]
  }
  sets
}

# colnames used internally in gen.lipidset
utils::globalVariables(c("collection", "value", "set"))

# colnames output by fgsea
utils::globalVariables(c("padj", "size", "pathway", "Significant", "Enrichment"))
