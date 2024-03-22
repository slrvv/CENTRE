#Test the createPairs function
test_that("createPairs() runs as expected", {
  benchmark <- readRDS(file= system.file("extdata",
                                         "input_generic_features.rds",
                                         package = "CENTRE"))
  genes <- as.data.frame(benchmark[, 1])

  testthat::expect_error(createPairs(genes))

  colnames(genes) <- c("gene_id")
  pair_data <- createPairs(genes)


  pair_data$pair <- paste(pair_data$enhancer_id, pair_data$gene_id1, sep = "_")
  benchmark$pair <- paste(benchmark[,2], benchmark[,1], sep = "_")

  #checking there are no duplicate pairs being returned
  testthat::expect_equal(length(unique(pair_data$pair)), nrow(pair_data))

  #check that the number of columns is correct
  testthat::expect_equal(ncol(pair_data), 3)

})
