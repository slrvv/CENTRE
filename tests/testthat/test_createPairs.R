#Test the createPairs function
test_that("createPairs() runs as expected", {
  benchmark <- readRDS(file= system.file("extdata",
                                         "input_generic_features.rds",
                                         package = "CENTRE"))
  genes <- as.data.frame(benchmark[, 1])

  pair_data <- createPairs(genes)

  pair_data$pair <- paste(pair_data$enhancer_id, pair_data$gene_id1, sep = "_")
  testthat::expect_equal(length(unique(pair_data$pair)), nrow(pair_data))
  testthat::expect_equal(ncol(pair_data), 3)
})
