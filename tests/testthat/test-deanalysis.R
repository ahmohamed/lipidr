test_that("Multiple group comparison works", {
  data("data_normalized")
  design <- model.matrix(~group, data = colData(data_normalized))
  de_results <- de_design(data_normalized, design = design)
  # de_results = de_analysis(data_normalized, group_col = "group")
  expect_true(all(
    c("AveExpr", "adj.P.Val", "Class", "Molecule") %in% colnames(de_results)
  ))
  expect_false("logFC" %in% colnames(de_results))
  
  p = plot_results_volcano(de.results = de_results)
  eval(p)  # we need to explicitly force-eval the plot
})

test_that("Pairwise will fail without contrasts", {
  data("data_normalized")
  expect_error(
    de_analysis(data_normalized, group_col = "group"),
    "No contrasts provided"
  )
})
