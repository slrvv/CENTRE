### Testing the genomic distances function

test_that("distances are correct for the input with gene and enhancer", {
  result1 <- read.table("~/Documentos/EPI/tests/testthat/test_result1_genomicdistances.txt",
                        header = T, stringsAsFactors = F,
                        colClasses = c("character", "character", "character"))
  example1 <- read.table("~/Documentos/EPI_v2/example1.txt")
  result_actual <- as.data.frame(distances(example1))
  expect_identical(result_actual, as.data.frame(result1))
})

test_that("distances are correct for the input with just genes", {
  result1 <- read.table("~/Documentos/EPI/tests/testthat/result2_test_genomic_distances.txt",
                        header = T, stringsAsFactors = F,
                        colClasses = c("character", "character", "character"))
  example1 <- read.table("~/Documentos/EPI_v2/example2.txt")
  result_actual <- as.data.frame(distances(example1))
  expect_identical(result_actual, as.data.frame(result1))
})


test_that("error when computing distances between differet chromosomes", {
   example <- as.data.frame(c("ENSG00000286832.1", "EH38E2776520"))
   expect_error(distances(example))
})
