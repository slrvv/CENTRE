

file1 <- "/project/CRUP_scores/toSara/Thyroid.GTEx-Benchmark.V38.edited.txt"
benchmark <- read.table(file = file1, header = T, sep = "\t")
genes <- benchmark[,1]
genes <- as.data.frame(genes[1:50])

test_that("createPairs() runs as expected",{
  pair_data <- createPairs(genes)
  pair_data$pair <- paste(pair_data$enhancer_id, pair_data$gene_id1, sep= "_")
  expect_equal(length(unique(pair_data$pair)), nrow(pair_data))
})
