library(lipidr)
de <- read_skyline(c("A1.csv", "F2.csv")) %>%
  add_sample_annotation(gen_sample_annot(.)) %>%
  summarize_transitions() %>%
  normalize_pqn() %>%
  de_analysis(A-B) %>%
  group_by_at(vars(logFC:B)) %>%
  summarise_all(first)

context("test-lsea")
test_that("should work with multiple sets and give correct results", {
  lsea(de, rank.by = "logFC")
  lsea(de, rank.by = "P.Value")

  de2 <- de
  de2[de2$Class == "PE", ]$logFC = seq_along(de2[de2$Class == "PE", ]$logFC) + 2
  res <- lsea(de2, rank.by = "logFC")
  expect_true("Class_PE" %in% significant_lipidsets(res))
  expect_true(res[res$set == "Class_PE", "NES"] > 0)
  de2[de2$Class == "PE", ]$logFC = seq_along(de2[de2$Class == "PE", ]$logFC) * -1
  res <- lsea(de2, rank.by = "logFC")
  expect_true("Class_PE" %in% significant_lipidsets(res))
  expect_true(res[res$set == "Class_PE", "NES"] < 0)

  de2[de2$Class == "PE", ]$adj.P.Val = seq_along(de2[de2$Class == "PE", ]$adj.P.Val) + 2
  res <- lsea(de2, rank.by = "adj.P.Val")
  expect_true("Class_PE" %in% significant_lipidsets(res))
  expect_true(res[res$set == "Class_PE", "NES"] > 0)
  de2[de2$Class == "PE", ]$adj.P.Val = seq_along(de2[de2$Class == "PE", ]$adj.P.Val) * -1
  res <- lsea(de2, rank.by = "logFC")
  expect_true("Class_PE" %in% significant_lipidsets(res))
  expect_true(res[res$set == "Class_PE", "NES"] < 0)
})

test_that("gen_lipidsets can accept de_results or a molecule list", {
  sets <- gen_lipidsets(de)
  sets2 <- gen_lipidsets(de$Molecule)
  lapply(names(sets), function(i) expect_setequal(sets[[i]], sets2[[i]]) )
})

test_that("Can filter sets based on min size", {
  set_size <- length(gen_lipidsets(de)$Class_LPE)
  expect_true("Class_LPE" %in% names(gen_lipidsets(de)))
  expect_true("Class_LPE" %in% names(gen_lipidsets(de, min_size = set_size)))
  expect_false("Class_LPE" %in% names(gen_lipidsets(de, min_size = set_size + 1)))

  expect_true("Class_LPE" %in% lsea(de)$set)
  expect_true("Class_LPE" %in% lsea(de, min_size = set_size)$set)
  expect_false("Class_LPE" %in% lsea(de, min_size = set_size + 1)$set)
})

test_that("Can filter sets with NA values", {
  de2 <- de
  de2$Class[de2$Class == "LPE"] <- NA
  expect_false("Class_LPE" %in% names(gen_lipidsets(de2)))
  expect_false("Class_LPE" %in% lsea(de2)$set)
})

test_that("Can filter sets that have all molecules", {
  de2 <- de
  de2$total_cs = 0
  expect_warning(
    sets <- gen_lipidsets(de2),
    "These sets contained all molecules, and were excluded: total_cs_0"
  )
  expect_false(
    any(grepl("total_cs", names(sets)))
  )
  expect_warning(
    res <- lsea(de2),
    "These sets contained all molecules, and were excluded: total_cs_0"
  )
  expect_false(
    any(grepl("total_cs", res$set))
  )
})

test_that("Use updated annotations when provided", {
  de2 <- de
  de2$Class[de2$Class == "LPE"] <- "CLS"
  expect_true("Class_CLS" %in% names(gen_lipidsets(de2)))
  expect_true("Class_CLS" %in% lsea(de2)$set)

  d <- read_skyline(c("A1.csv", "F2.csv")) %>%
    add_sample_annotation(gen_sample_annot(.)) %>%
    summarize_transitions() %>%
    normalize_pqn()

  rowData(d)$Class[1:3] <- "CLS"
  de2 <- d %>% de_analysis(A-B) %>%
    group_by_at(vars(logFC:B)) %>%
    summarise_all(first)

  expect_true("Class_CLS" %in% names(gen_lipidsets(de2)))
  expect_true("Class_CLS" %in% lsea(de2)$set)
})

test_that("Will generate error when no sets pass the filters", {
  expect_equal(length(gen_lipidsets(de, min_size = 22)), 0)
  expect_error(
    lsea(de, min_size = 22),
    "Unable to generate lipid sets, possibly because of missing annotations"
  )
  de2 <- de
  de2$Class <- "CLS"
  de2$total_cl <- 16
  de2$total_cs <- 0
  expect_warning(
    sets <- gen_lipidsets(de2),
    "These sets contained all molecules, and were excluded"
  )
  expect_equal(length(sets), 0)
  expect_error(
    suppressWarnings(lsea(de2)),
    "Unable to generate lipid sets, possibly because of missing annotations"
  )

  de2 <- de
  de2$Class <- NA
  de2$total_cl <- NA
  de2$total_cs <- NA
  sets <- gen_lipidsets(de2)
  expect_equal(length(sets), 0)
  expect_error(
    suppressWarnings(lsea(de2)),
    "Unable to generate lipid sets, possibly because of missing annotations"
  )
})

