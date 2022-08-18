################################################################################
# Testing compute generic features and compute cell type specific function
################################################################################

##Defining inputs for compute features function

files <- c("/project/CRUP_scores/CRUP_scores/ENCODE/Single_ended/Tissues/thyroid_gland/H3K4me1/H3K4me1.bam",
           "/project/CRUP_scores/CRUP_scores/ENCODE/Single_ended/Tissues/thyroid_gland/H3K4me3/H3K4me3.bam",
           "/project/CRUP_scores/CRUP_scores/ENCODE/Single_ended/Tissues/thyroid_gland/H3K27ac/H3K27ac.bam" )

inputs <- "/project/CRUP_scores/CRUP_scores/ENCODE/Single_ended/Tissues/thyroid_gland/Controls/input.bam"

metaData <- data.frame(HM = c("H3K4me1","H3K4me3","H3K27ac"),
                         condition = c(1,1,1), replicate = c(1,1,1),
                         bamFile = files, inputFile = rep(inputs,3))

tpmpath <- "/project/CRUP_scores/total-RNA-seq/thyroid-gland/ENCSR000AFK/thyroid-gland.tsv"
tpmfile <-  read.table(tpmpath, sep = "", stringsAsFactors = F, header = T)

###ASK TRIS FOR THE EXPECTED OUTPUT FILES BUT FOR NOW RUN HER CODE AND TAKE IT AS EXPECTED OUTPUT

expected.out1 <- read.table("/project/CRUP_scores/EPI/tests/output1.txt", header = T)


##Testing the compute generic features function

test_that("The warnings and errors in computeGenericFeatures() work",{
  file <- "/home/lopez_s/filedoesntexist"
  expect_condition(computeGenericFeatures(file))
  ## Warning for pairs in diff chromosomes
})

test_that("compute features functions runs as expected",{
  ###Checking that output has correct number of rows
  cat("Check ComputeGenericFeatures")
  file1 <- "/project/CRUP_scores/toSara/Thyroid.GTEx-Benchmark.V38.edited.txt"
  file1 <- read.table(file1, sep = "\t", stringsAsFactors = F)
  output <- computeGenericFeatures(file1)
  output$pair <- paste(output$enhancer_id, output$gene_id2, sep= "_")
  expect_equal(length(unique(output$pair)), nrow(output))
  cat("Check ComputeCellTypeFeatures")
  features_generic <- output
  output <- computeCellTypeFeatures(metaData, 1, "single", tpmfile, features_generic)
  output$pair <- paste(output$enhancer_id, output$gene_id2, sep= "_")
  expect_equal(length(unique(output$pair)), nrow(output))
  print(output$EP_prob_enh.1 !=  output$PP_prob_enh.1 )
  #### Checking that output is the expected

})





### Testing the genomic distances function



