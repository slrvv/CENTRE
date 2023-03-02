
test_that("compute generic feature function runs as expected", {
  ##Defining inputs for compute features function

  cat("Check ComputeGenericFeatures")
  input <- readRDS(file = system.file("extdata",
                                      "input_generic_features.rds",
                                      package = "CENTRE"))
  pred <- computeGenericFeatures(input)

  pred$pair <- paste(pred$enhancer_id,
                                 pred$gene_id2,
                                 sep = "_")
  expected <- readRDS(file = system.file("extdata",
                                         "expected_generic_features.rds",
                                         package = "CENTRE"))

  testthat::expect_equal(length(unique(pred$pair)), nrow(pred))
  testthat::expect_equal(dim(pred), c(19,9))
  testthat::expect_equal(pred$crup_cor, expected$cor_CRUP[,1])
  testthat::expect_equal(pred$combined_tests, expected$combined_tests, tolerance =1e-5)
})
