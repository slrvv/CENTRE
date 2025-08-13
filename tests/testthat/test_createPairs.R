#Test the createPairs function
test_that("createPairs() runs as expected", {
  benchmark <- readRDS(file= system.file("extdata",
                                         "input_generic_features.rds",
                                         package = "CENTRE"))
  benchmark2 <- readRDS(file= system.file("extdata",
                                          "output_enh_pairs.rds",
                                          package = "CENTRE"))
  genes <- benchmark[, 1]
  pair_data <- createPairs(genes)

  pair_data$pair <- paste(pair_data$enhancer_id, pair_data$gene_id1, sep = "_")
  benchmark$pair <- paste(benchmark[,2], benchmark[,1], sep = "_")

  #checking there are no duplicate pairs being returned
  testthat::expect_equal(length(unique(pair_data$pair)), nrow(pair_data))

  #check that the number of columns is correct
  testthat::expect_equal(ncol(pair_data), 3)
  
  enhancer <- c("EH38E3750708", "EH38E2776554")
  pair_enh <- createPairs(enhancer, enhancerCentered = TRUE)
  testthat::expect_equal(pair_enh, benchmark2)
  testthat::expect_error(createPairs())
})
