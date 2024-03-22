test_that("compute centrePrediction function runs as expected", {
  startPart("Check centrePrediction full input")
  generic_features <- readRDS(file = system.file("extdata",
                                                 "expected_generic_features.rds",
                                                 package = "CENTRE"))
  celltype_features <- readRDS(file = system.file("extdata",
                                                  "expected_celltype_features.rds",
                                                  package = "CENTRE"))
  ## add some pairs that are positive

  predictions_expected <- readRDS(file = system.file("extdata",
                                                     "expected_predictions.rds",
                                                     package = "CENTRE"))
  predictions <- centrePrediction(celltype_features,
                                generic_features)

  expect_equal(predictions$label, predictions_expected$label)
  expect_equal(predictions$score,
             predictions_expected$score,
             tolerance =1e-5) # changed name of feature

  ## only one pair
  startPart("Check centrePrediction only one pair")

  generic_features1 <- generic_features[generic_features$gene_id2 == "ENSG00000202276", ]
  celltype_features1 <- celltype_features[celltype_features$gene_id2 == "ENSG00000202276", ]
  predictions_expected1 <- predictions_expected[predictions_expected$pairs == "EH38E1519134_ENSG00000202276",]
  predictions1 <- centrePrediction(celltype_features1,
                                  generic_features1)
  expect_equal(predictions1$label, predictions_expected1$label)
  expect_equal(predictions1$score,
               predictions_expected1$score,
               tolerance =1e-5)
})


