context("Genomic distances")


test_that("distances are correct for the input with gene and enhancer",{
  result1 <- read.table("~/Documentos/EPI/tests/testthat/test_result1_genomicdistances.txt", header=T, stringsAsFactors = F)
  example1 <- read.table("~/Documentos/EPI_v2/example1.txt")
  expect_identical(distances(example1), as.data.frame(result1))
})
