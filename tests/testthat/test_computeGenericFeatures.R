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
    message("Test function with one pair and two pairs")
    ## gen features with only one pair
    input1 <- data.frame(
        gene_id1 = c("ENSG00000059728.6"),
        enhancer_id = c("EH38E3350767")
    )

    pred1 <- computeGenericFeatures(input1)

    ## gen features with two pairs
    input2 <- data.frame(
        gene_id1 = c("ENSG00000059728.6", "ENSG00000071677.1"),
        enhancer_id = c("EH38E3350767", "EH38E3410785")
    )
    colnames(input2) <- c("gene_id1", "enhancer_id")
    pred2 <- computeGenericFeatures(input2)


    ## different sized inputs still give correct result for combined_tests and
    ## cor_crup
    crup_cor_pred1 <- pred1 %>%
        dplyr::filter(gene_id1 == "ENSG00000059728") %>%
        dplyr::select(crup_cor)
    crup_cor_pred2 <- pred2 %>%
        dplyr::filter(gene_id1 == "ENSG00000059728") %>%
        dplyr::select(crup_cor)
    crup_cor_expected <- expected %>%
        dplyr::filter(gene_id1 == "ENSG00000059728") %>%
        dplyr::select(crup_cor)
    combined_test_pred1 <- pred1 %>%
        dplyr::filter(gene_id1 == "ENSG00000059728") %>%
        dplyr::select(combined_tests)
    combined_test_pred2 <- pred2 %>%
        dplyr::filter(gene_id1 == "ENSG00000059728") %>%
        dplyr::select(combined_tests)
    combined_test_expected <- expected %>%
        dplyr::filter(gene_id1 == "ENSG00000059728") %>%
        dplyr::select(combined_tests)

    testthat::expect_equal(
        crup_cor_pred1$crup_cor,
        crup_cor_expected$crup_cor
    )

    testthat::expect_equal(
        combined_test_pred1,
        combined_test_expected
    )

    testthat::expect_equal(
        crup_cor_pred2$crup_cor,
        crup_cor_expected$crup_cor
    )

    testthat::expect_equal(
        combined_test_pred2,
        combined_test_expected
    )
})
