mols <- c("PC(O-32:0)", "PC(O-16:0/16:0)", "So d18:0", "18:1(d7) LysoPC")

context("test-annotation")
test_that("Can match lipids that follow supported patterns", {
  a <- annotate_lipids(mols)
  expect_true(is.data.frame(a))
  expect_equal(nrow(a), 4)
  expect_false(any(a$not_matched))
})

test_that("Can work with factors", {
  a <- annotate_lipids(factor(mols))
  expect_true(is.data.frame(a))
  expect_equal(nrow(a), 4)
  expect_false(any(a$not_matched))
})

test_that("Cannot match lipids that do not follow supported patterns", {
  expect_warning(a <- annotate_lipids(c("a", "b", "cde", "any 15")), "Some lipid names couldn't be parsed")
  expect_true(is.data.frame(a))
  expect_equal(nrow(a), 4)
  expect_true(all(a$not_matched))
})

test_that("Can handle mixed matching and non-matching lipids", {
  expect_warning(a <- annotate_lipids(c(mols, "any 15")), "Some lipid names couldn't be parsed")
  expect_true(is.data.frame(a))
  expect_equal(nrow(a), 5)
  expect_equal(sum(a$not_matched), 1)
})

test_that("should not give warning on non-matching lipids if no_match=ignore", {
  expect_silent(a <- annotate_lipids(c(mols, "any 15"), no_match = "ignore"))
  expect_true(is.data.frame(a))
  expect_equal(nrow(a), 5)
  expect_equal(sum(a$not_matched), 1)
})

test_that("should not give warning on non-matching lipids if no_match=remove", {
  expect_silent(a <- annotate_lipids(c(mols, "any 15"), no_match = "remove"))
  expect_true(is.data.frame(a))
  expect_equal(nrow(a), 4)
  expect_equal(sum(a$not_matched), 0)
})
