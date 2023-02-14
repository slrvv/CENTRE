test_that("compute centrePrediction function runs as expected", {
  cat("Check centrePrediction")
  generic_features <- readRDS(file = system.file("extdata",
                                                 "expected_generic_features.rds",
                                                 package = "CENTRE"))
  celltype_features <- readRDS(file = system.file("extdata",
                                                  "expected_celltype_features.rds",
                                                  package = "CENTRE"))
  predictions_expected <- readRDS(file = system.file("extdata",
                                                     "expected_predictions.rds",
                                                     package = "CENTRE"))
  predictions <- centrePrediction(celltype_features,
                                generic_features)


  expect_equal(predictions$label, predictions_expected$label)
  expect_equal(predictions$predictions,
             predictions_expected$predictions,
             tolerance =1e-5)
})
