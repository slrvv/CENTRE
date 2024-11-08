#Test the createPairs function
test_that("createPairs() runs as expected", {
  benchmark <- readRDS(file= system.file("extdata",
                                         "input_generic_features.rds",
                                         package = "CENTRE"))
  benchmark2 <- readRDS(file= system.file("extdata",
                                          "output_enh_pairs.rds",
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
  
  enhancer <- data.frame(c("EH38E3750708", "EH38E2776554"))
  testthat::expect_error(createPairs(enhancer, enhancerCentered = TRUE))
  colnames(enhancer) <- "enhancer_id"
  pair_enh <- createPairs(enhancer, enhancerCentered = TRUE)
  testthat::expect_equal(pair_enh, benchmark2)
})
