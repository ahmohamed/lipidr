test_workflow <- function(d, clin_file=NULL, measure="Area", group_col=NULL, groups=NULL, ...) {

  if (!is.null(clin_file)) {
    d <- add_sample_annotation(d, clin_file)
  }
  p <- plot_samples(d, "tic", measure, color=group_col)
  p <- plot_samples(d, "boxplot", measure, color=group_col)
  p <- plot_molecules(d, "sd", measure, color="Class")
  p <- plot_molecules(d, "cv", measure, color="Class")
  p <- plot_molecules(d, "boxplot", measure, color="Class")
  p <- plot_lipidclass(d, "boxplot", measure)
  p <- plot_lipidclass(d, "sd", measure)
  p <- plot_heatmap(d, measure, sample_annotation = group_col, molecule_annotation = "Class")
  p <- plot_heatmap(d, measure, sample_annotation = "all", molecule_annotation = "Class")

  if (!is_summarized(d)) {
    summarize_transitions(d, method = "max")
    d <- summarize_transitions(d, method = "average")
  }
  if (any(rowData(d)$istd)) {
    normalize_istd(d, measure = measure)
  }
  d <- normalize_pqn(d, measure = measure)

  pca <- mva(d, measure, method = "PCA", group_col = group_col)
  p <- plot_mva(pca, color_by = group_col)

  pcoa <- mva(d, measure, method = "PCoA", group_col = group_col)
  p <- plot_mva(pcoa, color_by = group_col)

  if (!is.null(groups)) {
    oplsda <- mva(d, measure, method = "OPLS-DA", group_col = group_col, groups = groups)
    p <- plot_mva(oplsda, color_by = group_col)

    d$num_group <- factor(d[[group_col]]) %>% as.numeric()
    d_nona <- d[, !is.na(d$num_group)]
    opls <- mva(d_nona, measure, method = "OPLS-DA", group_col = group_col, groups = groups)
    p <- plot_mva(opls, color_by = group_col)

    comparison <- paste(groups, collapse = " - ")
    de_results <- de_analysis(d, ..., measure=measure, group_col = group_col)
    significant_molecules(de_results)
    p <- plot_results_volcano(de_results, TRUE)
    p <- plot_trend(de_results)

    sets <- gen_lipidsets(de_results)
    if (length(sets) > 0) {
      lsea(de_results, rank.by = "logFC")
      lsea(de_results, rank.by = "P.Value")
      enrich <- lsea(de_results, rank.by = "adj.P.Val")
      p <- plot_enrichment(de_results, significant_lipidsets(enrich), annotation="class")
      p <- plot_enrichment(de_results, significant_lipidsets(enrich), annotation="length")
      p <- plot_enrichment(de_results, significant_lipidsets(enrich), annotation="unsat")
      p <- plot_chain_distribution(de_results, contrast = comparison, measure = "logFC")
      p <- plot_chain_distribution(de_results, contrast = comparison, measure = "adj.P.Val")
    }
  }
  if (!is.null(group_col)) {
    de_results <- de_design(d, design = as.formula(paste0("~+", group_col)), measure=measure)
    significant_molecules(de_results)
    p <- plot_results_volcano(de_results, TRUE)
    # p <- plot_trend(de_results) TODO: use coefs for groups
  }
}
context("test-workflow")
test_that("Can run workflow with skyline non-pivot", {
  skip_on_bioc()
  file = "A1.csv"
  d <- read_skyline(file)
  test_workflow(d)
})

test_that("Can run workflow with multiple skyline non-pivot", {
  skip_on_bioc()
  file = list("A1.csv", "F2.csv")
  d <- read_skyline(file)
  test_workflow(d)
})

test_that("Can run workflow with skyline non-pivot and sample annotations", {
  skip_on_bioc()
  file = "A1.csv"
  d <- read_skyline(file)
  df <- gen_sample_annot(d)
  test_workflow(d, clin_file=df, measure="Area", group_col = "Group", groups = c("A", "B"), A-B)
})

test_that("Can run workflow with skyline pivot and sample annotations", {
  skip_on_bioc()
  file = "A1_pivot.csv"
  d <- read_skyline(file)
  df <- gen_sample_annot(d)
  test_workflow(d, clin_file=df, measure="Area", group_col = "Group", groups = c("A", "B"), A-B)
})

test_that("Can run workflow with with measure other than Area", {
  skip_on_bioc()
  file = .read_tabular("A1.csv") %>% rename("Normalized Area"=Area) %>% save_temp_csv()
  d <- read_skyline(file)
  df <- gen_sample_annot(d)
  test_workflow(d, clin_file=df, measure="Normalized Area", group_col = "Group", groups = c("A", "B"), A-B)
})

# num matrix
test_that("Can run workflow with with num matrix", {
  skip_on_bioc()
  file = "A1.csv"
  d <- read_skyline(file)
  d <- assay(d, "Area") %>% `row.names<-`(rowData(d)$Molecule) %>% as_lipidomics_experiment()
  test_workflow(d)
})

# num matrix with annot
test_that("Can run workflow with with num matrix with sample annot", {
  skip_on_bioc()
  file = .read_tabular("A1.csv") %>% rename("Normalized Area"=Area) %>% save_temp_csv()
  d <- read_skyline(file)
  df <- gen_sample_annot(d)
  d <- assay(d, "Normalized Area") %>% `row.names<-`(rowData(d)$Molecule) %>% as_lipidomics_experiment()
  test_workflow(d, clin_file=df, measure="Area", group_col = "Group", groups = c("A", "B"), A-B)
})

# group_col missing vals
test_that("Can run workflow with where some sample annot are missing", {
  skip_on_bioc()
  file = .read_tabular("A1.csv") %>% rename("Normalized Area"=Area) %>%
    mutate(Peptide=ifelse(1:nrow(.) %% 2 == 0,  NA, Peptide)) %>%
    save_temp_csv(quote=FALSE, na = "")

  d <- read_skyline(file)
  df <- gen_sample_annot(d) %>% mutate(Group=ifelse(1:nrow(.) %% 2 == 0,  NA, Group))
  d <- assay(d, "Normalized Area") %>% `row.names<-`(rowData(d)$Molecule) %>% as_lipidomics_experiment()
  test_workflow(d, clin_file=df, measure="Area", group_col = "Group", groups = c("A", "B"), A-B)
})

# Mols missing vals
test_that("Can run workflow with where some molecules are missing", {
  skip_on_bioc()
  file = .read_tabular("A1.csv") %>% rename("Normalized Area"=Area) %>%
    mutate(Peptide=ifelse(1:nrow(.) %% 10 == 0,  NA, Peptide)) %>%
    save_temp_csv(quote=FALSE, na = "")

  d <- read_skyline(file)
  df <- gen_sample_annot(d)
  d <- assay(d, "Normalized Area") %>% `row.names<-`(rowData(d)$Molecule) %>% as_lipidomics_experiment()
  test_workflow(d, clin_file=df, measure="Area", group_col = "Group", groups = c("A", "B"), A-B)
})

# all non-lipid
test_that("Can run workflow with all non-lipid molecules in Skyline format", {
  skip_on_bioc()
  file = .read_tabular("A1.csv") %>% rename("Normalized Area"=Area) %>%
    mutate(Peptide=paste0("any", Peptide, "any")) %>%
    save_temp_csv(quote=FALSE, na = "")

  d <- read_skyline(file)
  df <- gen_sample_annot(d)
  test_workflow(d, clin_file=df, measure="Normalized Area",group_col = "Group", groups = c("A", "B"), A-B)
})
