### Testing the genomic distances function

test_that("distances are correct for the input with gene and enhancer", {
  correct_result <- 750312
  example1 <- as.data.frame(cbind(c("ENSG00000227232.5"), c("EH38E2776542")))
  result_function <- as.data.frame(distances(example1))
  expect_identical(result_function$distance, correct_result)
})

# test_that("distances are correct for the input with just genes", {
#   result1 <- read.table("~/Documentos/EPI/tests/testthat/result2_test_genomic_distances.txt",
#                         header = T, stringsAsFactors = F,
#                         colClasses = c("character", "character", "character"))
#   example1 <- read.table("~/Documentos/EPI_v2/example2.txt")
#   result_actual <- as.data.frame(distances(example1))
#   expect_identical(result_actual, as.data.frame(result1))
# })


test_that("warning when computing distances between differet chromosomes", {
   example <- as.data.frame(cbind(c("ENSG00000286832.1"), c("EH38E3565295")))
   expect_warning(distances(example))
})
