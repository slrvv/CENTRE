test_that("compute generic feature function runs as expected", {
    message("Test that errors are raised as expected")
    ## test errors are raised
    testthat::expect_error(computeGenericFeatures())

    error_dat <- data.frame(
        gene_id = c("ENSG00000059728.6"),
        enh_id = c("EH38E3350767")
    )
    testthat::expect_error(computeGenericFeatures(error_dat))


    # rename the same dataframe used for createPairs for the sake of storing less
    # data


    ## Expected result for this function
    expected <- readRDS(file = system.file("extdata",
        "expected_generic_features.rds",
        package = "CENTRE"
    ))

    input <- expected[, c("gene_id1", "enhancer_id")]
    ## generate features with all of the input
    message("Test computeGenericFeatures with full input")
    pred <- computeGenericFeatures(input)

    pred$pair <- paste(pred$enhancer_id,
        pred$gene_id1,
        sep = "_"
    )
    expected$pair <- paste(expected$enhancer_id,
        expected$gene_id1,
        sep = "_"
    )
    ## ordered for comparison
    expected <- expected[order(expected$pair), ]
    pred <- pred[order(pred$pair), ]

    ## pairs are unique and the dimensions ob the returned dataset are the expected ones
    testthat::expect_equal(length(unique(pred$pair)), nrow(pred))
    testthat::expect_equal(dim(pred), c(19, 6))
    # columns return expected values

    testthat::expect_equal(pred$crup_cor, expected$crup_cor)
    testthat::expect_equal(pred$distance, abs(expected$distance))
    testthat::expect_equal(pred$combined_tests,
        expected$combined_tests,
        tolerance = 1e-5
    )
})
