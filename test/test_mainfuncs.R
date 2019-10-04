library(testthat)
library("UbiquitinAnalysis")

test_that("To test GetSimilarPercentage", {
  protein1 = "Q9UHB7"
  protein2 = "Q9UKV5"

  value <- GetSimilarPercentage(protein1, protein2)

  expect_that(value, is_a("percentage"))
  expect_that(trunc(value), equals("56.8%"))
})

test_check("UbiquitinAnalysis")

