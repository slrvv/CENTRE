context("Compute features function")

##Defining inputs for compute features function

files <- c("/project/CRUP_scores/CRUP_scores/ENCODE/Single_ended/Tissues/body_pancreas/H3K4me1/H3K4me1.bam",
           "/project/CRUP_scores/CRUP_scores/ENCODE/Single_ended/Tissues/body_pancreas/H3K4me3/H3K4me3.bam",
           "/project/CRUP_scores/CRUP_scores/ENCODE/Single_ended/Tissues/body_pancreas/H3K27ac/H3K27ac.bam" )

inputs <- "/project/CRUP_scores/CRUP_scores/ENCODE/Single_ended/Tissues/body_pancreas/Controls/input.bam"

metaData <- data.frame(HM = c("H3K4me1","H3K4me3","H3K27ac"),
                         condition = c(1,1,1), replicate = c(1,1,1),
                         bamFile = files, inputFile = rep(inputs,3))
tpmfile <- "/project/CRUP_scores/total-RNA-seq/body-of-pancreas/ENCSR586SYA/body-of-pancreas.tsv"
expected.out1 <- read.table("/project/CRUP_scores/EPI/tests/output1.txt", header = T)
#expected.out2 <-

##Testing the compute features function

test_that("The warnings and errors in compute_features() work",{
  file <- "/home/lopez_s/filedoesntexist"
  expect_condition(compute_features(file, metaData, 1, tpmfile))
  example <- as.data.frame(cbind(c("ENSG00000286832.1"), c("EH38E3565295")))
  expect_warning(distances_gene_enhancer(example))
})

test_that("compute_features() runs as expected",{
  file1 <- "/project/CRUP_scores/EPI/example1.txt"
  file2 <- "/project/CRUP_scores/EPI/example2.txt"
  expect_equal(compute_features(file1, metaData, 1, tpmfile), expected.out1, tolerance = 1e-5)
  #expect_equal(compute_features(file2, metaData, 3), expected.out1, tolerance = 1e-5)
})



### Testing the genomic distances function

test_that("distances are correct for the input with gene and enhancer", {
  correct_result <- 750312
  example1 <- as.data.frame(cbind(c("ENSG00000227232.5"), c("EH38E2776542")))
  result_function <- as.data.frame(distances_gene_enhancer(example1))
  expect_identical(result_function$distance, correct_result)
})

