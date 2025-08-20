# Test the createPairs function
test_that("createPairs() runs as expected", {
    benchmarkPairsGenes <- readRDS(file = system.file("extdata",
        "expected_gene_pairs.rds",
        package = "CENTRE"
    ))

    benchmarkPairsEnh <- readRDS(file = system.file("extdata",
        "expected_enh_pairs.rds",
        package = "CENTRE"
    ))

    message("Test that createPairs() from genes works as expected.\n")
    genes <- c("ENSG00000105281", "ENSG00000180458")

    pairComp <- createPairs(genes)

    pairComp$pair <- paste(pairComp$enhancer_id, pairComp$gene_id1, sep = "_")
    benchmarkPairsGenes$pair <- paste(benchmarkPairsGenes$enhancer_id,
        benchmarkPairsGenes$gene_id,
        sep = "_"
    )
    benchmarkPairsGenes <- benchmarkPairsGenes[order(benchmarkPairsGenes$enhancer_id), ]
    pairComp <- pairComp[order(pairComp$enhancer_id), ]
    # checking there are no duplicate pairs being returned
    testthat::expect_equal(length(unique(pairComp$pair)), nrow(pairComp))
    testthat::expect_equal(pairComp, benchmarkPairsGenes)

    # check that the number of columns is correct
    testthat::expect_equal(ncol(pairComp), 3)

    message("Test that createPairs() from enhancers works as expected.\n")
    enhancer <- c("EH38E3750708", "EH38E2776554")
    pairCompEnh <- createPairs(enhancer, enhancerCentered = TRUE)
    # ordering for comparison
    benchmarkPairsEnh <- benchmarkPairsEnh[order(benchmarkPairsEnh$gene_id1), ]
    pairCompEnh <- pairCompEnh[order(pairCompEnh$gene_id1), ]
    testthat::expect_equal(pairCompEnh, benchmarkPairsEnh)
    # test that it throws error when it should
    testthat::expect_error(createPairs())
})
