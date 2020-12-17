test_that("Interactive plots work", {
  use_interactive_graphics()
  data("data_normalized")
  plot_molecules(data_normalized, 'cv')
  plot_molecules(data_normalized, 'sd')
  plot_molecules(data_normalized, 'boxplot')
  
  plot_samples(data_normalized, 'tic')
  plot_samples(data_normalized, 'boxplot')
  
  plot_lipidclass(data_normalized, 'sd')
  plot_lipidclass(data_normalized, 'boxplot')
  expect_true(TRUE) # Needed so devtools doesn't think it's an empty test.
})
