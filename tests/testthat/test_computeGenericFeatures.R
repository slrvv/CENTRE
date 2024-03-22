
test_that("compute generic feature function runs as expected", {
<<<<<<< Updated upstream
  ##Defining inputs for compute features function

  cat("Check ComputeGenericFeatures")
  input <- readRDS(file = system.file("extdata",
                                      "input_generic_features.rds",
                                      package = "CENTRE"))
=======
  #Defining inputs for compute features function
  input <- readRDS(file = system.file("extdata",
                                      "input_generic_features.rds",
                                      package = "CENTRE"))

  ##Expected result for this function
  expected <- readRDS(file = system.file("extdata",
                                         "expected_generic_features.rds",
                                         package = "CENTRE"))

  ## gen features with all of the input
>>>>>>> Stashed changes
  pred <- computeGenericFeatures(input)

  pred$pair <- paste(pred$enhancer_id,
                                 pred$gene_id2,
                                 sep = "_")
<<<<<<< Updated upstream
  expected <- readRDS(file = system.file("extdata",
                                         "expected_generic_features.rds",
                                         package = "CENTRE"))

  testthat::expect_equal(length(unique(pred$pair)), nrow(pred))
  testthat::expect_equal(dim(pred), c(19,9))
  testthat::expect_equal(pred$crup_cor, expected$cor_CRUP[,1])
  testthat::expect_equal(pred$combined_tests, expected$combined_tests, tolerance =1e-5)
=======

  ## gen features with only one pair
  input1 <- as.data.frame(matrix(c("ENSG00000059728.6", "EH38E3350767"),
                          nrow = 1,
                          ncol = 2,
                          byrow = T))

  pred1 <- computeGenericFeatures(input1)

  ## gen features with two pairs
  input2 <- as.data.frame(matrix(c("ENSG00000059728.6", "EH38E3350767",
                                   "ENSG00000071677.1" ,"EH38E3410785"),
                                 nrow = 2,
                                 ncol = 2,
                                 byrow = T))

  pred2 <- computeGenericFeatures(input2)

  ## pairs are unique and the dimensions ob the returned dataset are the expected ones
  testthat::expect_equal(length(unique(pred$pair)), nrow(pred))
  testthat::expect_equal(dim(pred), c(19,6))

  testthat::expect_equal(pred$crup_cor, expected$crup_cor[,1])
  testthat::expect_equal(pred$distance, expected$distance)
  testthat::expect_equal(pred$combined_tests,
                         expected$combined_tests,
                         tolerance =1e-5)
  ## different sized inputs still give correct result for combined_tests and
  ## cor_crup
  testthat::expect_equal(pred1[pred1$gene_id2 == "ENSG00000059728", 4],
                         expected[expected$gene_id2 == "ENSG00000059728", 7][,1])

  testthat::expect_equal(pred1[pred1$gene_id2 == "ENSG00000059728", 5],
                         expected[expected$gene_id2 == "ENSG00000059728", 8],
                         tolerance =1e-5)

  testthat::expect_equal(pred2[pred2$gene_id2 == "ENSG00000059728", 4],
                         expected[expected$gene_id2 == "ENSG00000059728", 7][,1])

  testthat::expect_equal(pred2[pred2$gene_id2 == "ENSG00000059728", 5],
                         expected[expected$gene_id2 == "ENSG00000059728", 8],
                         tolerance =1e-5)
>>>>>>> Stashed changes
})
