#' @importFrom ggplot2 ggplot aes_string aes facet_wrap theme theme_dark theme_grey guides coord_flip
#' @importFrom ggplot2 stat_sum stat_summary scale_fill_gradient2 scale_color_gradient2 element_text
#' @importFrom ggplot2 geom_point geom_boxplot geom_bar geom_text geom_tile geom_vline geom_hline xlab ylab labs
#' @importFrom ggplot2 scale_color_manual stat_ellipse annotate
{} 


#' Plot a bar chart for total sample intensity
#'
#' @param data Skyline data.frame created by \code{\link{read_skyline}}
#' @param measure which measure to use as intensity, usually Area, Area.Normalized or Height
#' @param log whether values should be log2 transformed
#'
#' @export
#'
#' @examples
#' datadir = system.file("extdata", package="lipidr")
#' filelist = list.files(datadir, "data.csv", full.names = TRUE)
#' d = read_skyline(filelist)
#' 
#' plot_sample_tic(d, "Area", log=TRUE)
#' plot_sample_tic(d, "Background", log=FALSE)
plot_sample_tic <- function(data, measure="Area", log=TRUE){
  stopifnot(inherits(data, "SkylineExperiment"))
  dlong = to_long_format(data, measure)
  if (log == TRUE) {
    measure = .check_log(data, measure)
  }
  ggplot(dlong, aes_string("Sample", measure)) + stat_sum(geom="bar") + 
    facet_wrap(~filename, ncol=1, scales="free_y") + 
    theme(axis.text.x = element_text(angle = -90, vjust = 0.5)) + guides(size=FALSE)
}

#' Plot a boxplot chart to examine distribution of values per sample
#'
#' The function is usually used to look at intensity distribution in each sample
#' ensuring they are normalized. It can also be used to look at different measures 
#' such as  `Retention.Time` or `Background`.
#'
#' @param data Skyline data.frame created by \code{\link{read_skyline}}
#' @param measure which measure to plot the distribution of: usually Area, Area.Normalized or Height
#' @param log whether values should be log2 transformed
#'
#' @export
#' @examples
#' datadir = system.file("extdata", package="lipidr")
#' filelist = list.files(datadir, "data.csv", full.names = TRUE)
#' d = read_skyline(filelist)
#' 
#' plot_sample_boxplot(d, "Area", log=TRUE)
#' plot_sample_boxplot(d, "Retention.Time", log=FALSE)
plot_sample_boxplot <- function(data, measure="Area", log=TRUE){
  stopifnot(inherits(data, "SkylineExperiment"))
  dlong = to_long_format(data, measure)
  if (log == TRUE) {
    measure = .check_log(data, measure)
  }
  
  ggplot(dlong, aes_string("Sample", measure)) + geom_boxplot() + 
    facet_wrap(~filename, ncol=1, scales="free_y") + 
    theme(axis.text.x = element_text(angle = -90, vjust = 0.5)) + guides(size=FALSE)
}



#' Plot a bar chart for standard deviation of a certain measure in each class
#'
#' The function is usually used to look at standard deviations of intensity in each class,
#' but can also be used to look at different measures such as  `Retention.Time`,
#' to ensure all lipids are eluted within expected range. 
#' To assess instrumental variation apply the function to technical quality control samples.
#'
#' @param data Skyline data.frame created by \code{\link{read_skyline}}
#' @param measure which measure to plot the distribution of: usually Area, Area.Normalized, Height or Retention.Time
#' @param log whether values should be log2 transformed (Set FALSE for retention time).
#'
#' @export
#' @examples 
#' datadir = system.file("extdata", package="lipidr")
#' filelist = list.files(datadir, "data.csv", full.names = TRUE)
#' d = read_skyline(filelist)
#' 
#' plot_class_sd(d, "Area", log=TRUE)
#' plot_class_sd(d, "Retention.Time", log=FALSE)
plot_class_sd <- function(data, measure="Area", log=TRUE){
  stopifnot(inherits(data, "SkylineExperiment"))
  dlong = to_long_format(data, measure)
  if (log == TRUE) {
    measure = .check_log(data, measure)
  }
  
  ggplot(dlong, aes_string("Class", measure, fill="Class")) + 
    stat_summary(fun.y = sd, geom="bar")  + 
    facet_wrap(~filename, scales="free_x") + 
    theme(axis.text.x = element_text(angle = -90, vjust = 0.5)) + 
    ylab(paste("SD of", measure))
}

#' Plot a boxplot chart to examine distribution of values per class
#'
#' The function is usually used to look at intensity distribution in each class,
#' but can also be used to look at different measures, such as  `Retention.Time` or `Background`.
#'
#' @param data Skyline data.frame created by \code{\link{read_skyline}}
#' @param measure which measure to plot the distribution of: usually Area, Area.Normalized or Height
#' @param log whether values should be log2 transformed
#'
#' @export
#' @examples 
#' datadir = system.file("extdata", package="lipidr")
#' filelist = list.files(datadir, "data.csv", full.names = TRUE)
#' d = read_skyline(filelist)
#' 
#' plot_class_boxplot(d, "Area", log=TRUE)
#' plot_class_boxplot(d, "Retention.Time", log=FALSE)
plot_class_boxplot <- function(data, measure="Area", log=TRUE){
  stopifnot(inherits(data, "SkylineExperiment"))
  dlong = to_long_format(data, measure)
  if (log == TRUE) {
    measure = .check_log(data, measure)
  }
  
  ggplot(dlong, aes_string("Class", measure, fill="Class")) + 
    geom_boxplot()  + 
    facet_wrap(~filename, scales="free_x") + 
    theme(axis.text.x = element_text(angle = -90, vjust = 0.5))
}


#' Plot a boxplot chart to log fold changes of lipids per class
#'
#' The function is usually used to look at logFC distribution of lipids in each class,
#' marking significantly enriched classes. Can also be used to plot `P.Value` or `Adj.P.Val`
#'
#' @param de_results output of \code{\link{de_analysis}}
#' @param significant.sets list of significantly changed lipid sets (output of \code{\link{significant_lipidsets}})
#' @param measure which measure to plot the distribution of: logFC, P.Value, Adj.P.Val
#'
#' @export
plot_class_enrichment <- function(de_results, significant.sets, measure="logFC"){
  significant.sets = lapply(significant.sets,
    function(c) sub("^Class_", "", c[grep("^Class_", c)] )
  )
  de_results = de_results %>% annotate_lipids() %>% 
    .left_join.silent(de_results) %>%
    group_by(contrast) %>%
    mutate(Significant=Class %in% significant.sets[[ contrast[[1]] ]]) %>%
    ungroup()
  
  ggplot(de_results, aes_string("Class", measure, color="Significant")) + 
    geom_boxplot()  + geom_hline(yintercept=0, lty=2) + 
    facet_wrap(~contrast, scales="free_x") + 
    scale_color_manual(values=c("black", "red"))+
    theme(axis.text.x = element_text(angle = -90, vjust = 0.5))
}

#' Plot a chart for log fold changes of lipids per class showing chain lengths and saturations
#'
#' @param de_results output of \code{\link{de_analysis}}
#' @param measure which measure to plot the distribution of: logFC, P.Value, Adj.P.Val
#' @param contrast which comparison to plot.
#'
#' @export
plot_chain_distribution <- function(de_results, contrast=NULL, measure="logFC"){
  if(is.null(contrast)) {
    contrast = de_results$contrast[[1]]
  }
  measure = sym(measure)
  
  de_results = de_results %>% 
    select(Molecule, !!measure) %>%
    annotate_lipids() %>% 
    filter(!itsd) %>%
    .left_join.silent(de_results) %>%
    group_by(clean_name) %>%
    ungroup
  
  de_results = de_results[de_results$contrast==contrast, ] %>%
    mutate_at(vars(total_cl:total_cs), factor)
  
  ggplot(de_results, aes(total_cs, total_cl, fill=logFC)) + geom_tile() + 
    facet_wrap(~Class) + 
    xlab("Total chain unsaturation") + ylab("Total chain length") +
    scale_fill_gradient2(midpoint = 0)
}

#' Plot a bar chart for standard deviation if a certain measure in each lipid
#'
#' The function is usually used to look at standard deviation of intensity for each lipid,
#' but can also be used to look at different measures such as  `Retention.Time`,
#' to ensure all lipids elute within expected range.
#'
#' @param data Skyline data.frame created by \code{\link{read_skyline}}
#' @param measure Which measure to plot the distribution of: usually Area, Area.Normalized or Height
#' @param log Whether values should be log2 transformed (Set FALSE for retention time).
#'
#' @export
#' @examples 
#' datadir = system.file("extdata", package="lipidr")
#' filelist = list.files(datadir, "data.csv", full.names = TRUE)
#' d = read_skyline(filelist)
#' 
#' plot_molecule_sd(d)
#' plot_molecule_sd(d, "Retention.Time", log=FALSE)
plot_molecule_sd <- function(data, measure="Area", log=TRUE){
  stopifnot(inherits(data, "SkylineExperiment"))
  dlong = to_long_format(data, measure)
  if (log == TRUE) {
    measure = .check_log(data, measure)
  }
  ggplot(dlong, aes_string("Molecule", measure, fill="Class", color="Class")) + 
    stat_summary(fun.y = sd, geom="bar") +
    facet_wrap(~filename, scales="free_y") + coord_flip() +
    theme(axis.text.x = element_text(angle = -90, vjust = 0.5)) + 
    ylab(paste("SD of", measure))
}

plot_molecule_cv <- function(data, measure="Area", log=TRUE){
  stopifnot(inherits(data, "SkylineExperiment"))
  dlong = to_long_format(data, measure)
  if (log == TRUE) {
    measure = .check_log(data, measure)
  }
  ggplot(dlong, aes_string("Molecule", measure, fill="Class", color="Class")) + 
    stat_summary(fun.y = .cv, geom="bar") + coord_flip() +
    facet_wrap(~filename, scales="free_y") + 
    theme(axis.text.x = element_text(angle = -90, vjust = 0.5)) + 
    ylab(paste("SD of", measure))
}

#' Plot a boxplot chart to examine distribution of values per lipid
#'
#' The function is usually used to look at intensity distribution for each lipid,
#' but can also be used to look at different measures, such as  `Retention.Time` or `Background`.
#'
#' @param data Skyline data.frame created by \code{\link{read_skyline}}
#' @param measure which measure to plot the distribution of: usually Area, Area.Normalized or Height
#' @param log whether values should be log2 transformed
#'
#' @export
#' @examples 
#' datadir = system.file("extdata", package="lipidr")
#' filelist = list.files(datadir, "data.csv", full.names = TRUE)
#' d = read_skyline(filelist)
#' 
#' plot_molecule_boxplot(d)
#' plot_molecule_boxplot(d, "Retention.Time", log=FALSE)
plot_molecule_boxplot <- function(data, measure="Area", log=TRUE){
  stopifnot(inherits(data, "SkylineExperiment"))
  dlong = to_long_format(data, measure)
  if (log == TRUE) {
    measure = .check_log(data, measure)
  }
  ggplot(dlong, aes_string("Molecule", measure, fill="Class", color="Class")) + 
    geom_boxplot(outlier.size=0.5, outlier.alpha=0.3)  + coord_flip() + 
    facet_wrap(~filename, scales="free_y") + 
    theme(axis.text.x = element_text(angle = -90, vjust = 0.5)) 
}

#' Plot a volcano chart for differential analysis results.
#'
#' @param de_results output of \code{\link{de_analysis}}
#' @param show.labels whether labels show be displayed for significant lipids
#' 
#' @export
plot_results_volcano <- function(de_results, show.labels=TRUE) {
  de_results %>%
    mutate_at(vars(matches("P.Val")), log10) %>%
    (function(.){ 
      p = ggplot(., aes(logFC, -adj.P.Val, color=Class, label=Molecule)) + 
       geom_point() + 
       geom_hline(yintercept = -log10(0.05), lty=2) + 
       geom_vline(xintercept = c(1, -1), lty=2) + 
       facet_wrap(~contrast)
     if (show.labels) {
       p + geom_text(aes(label=ifelse(adj.P.Val < log10(0.05) & abs(logFC) > 1, Molecule, "")), vjust=-.5, size=3, color="black")
     }
     return(p)
    })
}

