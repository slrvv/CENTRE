test_that("compute centrePrediction function runs as expected", {
    message("Check centrePrediction full input")
    message("Get generic features...")

    pairs <- readRDS(file = system.file("extdata",
        "input_cellType_pairs.rds",
        package = "CENTRE"
    ))

    generic_features <- computeGenericFeatures(pairs)

    celltype_features <- readRDS(file = system.file("extdata",
        "expected_cellType_HeLa_reduced.rds",
        package = "CENTRE"
    ))

    predictions_expected <- readRDS(file = system.file("extdata",
        "expected_predictions_hela.rds",
        package = "CENTRE"
    ))
    predictions <- centrePrediction(
        celltype_features,
        generic_features
    )

    expect_equal(predictions$label, predictions_expected$label)
    expect_equal(predictions$score,
        predictions_expected$score,
        tolerance = 1e-5
    ) # changed name of feature

})
