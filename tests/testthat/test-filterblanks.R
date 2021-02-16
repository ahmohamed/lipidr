library(lipidr)

d <- read_skyline(c("A1.csv", "F2.csv"))



expect_assay_equal <- function(e1, e2, measure) {
  mat <- assay(e1, measure)
  rownames(mat) <- rowData(e1)$Sample
  
  sum_mat <- assay(e2, measure)
  rownames(sum_mat) <- rowData(e2)$Molecule
  expect_equal(mat, sum_mat[rownames(mat), ] )
}



context("test-filterblanks")

test_that("filter_by_blanks does not work when there are no blanks", {
  
  d <- read_skyline(c("A1.csv", "F2.csv"))
  
  expect_valid_lipidex(d, c(22,11))
  
  blanks <- colnames(d) %in% c("Blank_1", "Blank_2")
  d_no_blanks <- d[, !blanks]
  expect_error(filter_by_blanks(d_no_blanks), "there are no blanks")
  
})
  
  
  
test_that("removes molecules below intensity threshold", {
  
  d <- read_skyline(c("A1.csv", "F2.csv"))
  
  expect_valid_lipidex(d, c(22,11))
  
  blanks <- colnames(d) %in% c("Blank_1", "Blank_2")
  feature <- rownames(d)[[1]]
 
  area <- assay(d, "Area")
  area[feature, blanks] <- mean(area[feature, !blanks], na.rm = TRUE)
  assay(d, "Area") <- area
  
  area <- assay(d, "Area")
  area[feature, blanks] <- area[feature, blanks] / 5
  assay(d, "Area") <- area
  
  d_filtered <- filter_by_blanks(d)
  
  expect_valid_lipidex(d_filtered, c(21,11))
  
})


test_that("does not effect molecules above the intensity threshold", {
  
  d <- read_skyline(c("A1.csv", "F2.csv"))
  
  expect_valid_lipidex(d, c(22,11))
  blanks <- colnames(d) %in% c("Blank_1", "Blank_2")
  feature <- rownames(d)[[1]]
  
  area <- assay(d, "Area")
  area[feature, blanks] <- mean(area[feature, !blanks], na.rm = TRUE)
  assay(d, "Area") <- area
  
  area <- assay(d, "Area")
  area[feature, blanks] <- area[feature, blanks] / 20
  assay(d, "Area") <- area
  
  d_filtered <- filter_by_blanks(d)
  
  expect_valid_lipidex(d_filtered, c(22,11))
  
})

test_that("filter_by_blanks works with 'Retention Time' as measure", {
  
  d <- read_skyline(c("A1.csv", "F2.csv"))
  
  expect_valid_lipidex(d, c(22,11))
  
  blanks <- colnames(d) %in% c("Blank_1", "Blank_2")
  feature <- rownames(d)[[1]]
  
  
  
  r_time <- assay(d, "Retention Time")
  r_time[feature, blanks] <- mean(r_time[feature, !blanks], na.rm = TRUE)
  assay(d, "Retention Time") <- r_time 
  
  r_time <- assay(d, "Retention Time")
  r_time[feature, blanks] <- r_time[feature, blanks] / 20
  assay(d, "Retention Time") <- r_time
  
  d_filtered <- filter_by_blanks(d, measure = "Retention Time")
  
  expect_valid_lipidex(d_filtered, c(2, 11))
  
})
