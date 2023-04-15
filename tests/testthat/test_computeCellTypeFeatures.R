################################################################################
# Testing that the pipeline works
################################################################################




##Testing the compute generic features function

test_that("computeCellTypeFeatures functions runs as expected for thyroid", {
  ##Defining inputs for compute features function

  files <- c("/project/CRUP_scores/CRUP_scores/ENCODE/Single_ended/Tissues/thyroid_gland/H3K4me1/H3K4me1.bam",
             "/project/CRUP_scores/CRUP_scores/ENCODE/Single_ended/Tissues/thyroid_gland/H3K4me3/H3K4me3.bam",
             "/project/CRUP_scores/CRUP_scores/ENCODE/Single_ended/Tissues/thyroid_gland/H3K27ac/H3K27ac.bam")

  inputs <- "/project/CRUP_scores/CRUP_scores/ENCODE/Single_ended/Tissues/thyroid_gland/Controls/input.bam"

  metaData <- data.frame(HM = c("H3K4me1", "H3K4me3", "H3K27ac"),
                         condition = c(1, 1, 1), replicate = c(1, 1, 1),
                         bamFile = files, inputFile = rep(inputs, 3))

  tpmpath <- "/project/CRUP_scores/total-RNA-seq/thyroid-gland/ENCSR000AFK/thyroid-gland.tsv"
  tpmfile <-  read.table(tpmpath, sep = "", stringsAsFactors = F, header = T)

  ###expected_predictions
  generic_features<- readRDS(file = system.file("extdata",
                                         "expected_generic_features.rds",
                                         package = "CENTRE"))
  expcelltype_features <- readRDS(file = system.file("extdata",
                                                  "expected_celltype_features.rds",
                                                  package = "CENTRE"))


  cat("Check ComputeCellTypeFeatures")
  celltype_features <- computeCellTypeFeatures(metaData,
                                               replicate = 1,
                                               input.free = FALSE,
                                               cores = 1,
                                               sequencing = "single",
                                               tpmfile = tpmfile,
                                               featuresGeneric = generic_features)
  celltype_features$pair <- paste(celltype_features$enhancer_id,
                                  celltype_features$gene_id2,
                                  sep = "_")
  expect_equal(length(unique(celltype_features$pair)), nrow(celltype_features))

  testthat::expect_equal(celltype_features$TPM,
                         expcelltype_features$TPM,
                         tolerance = 1e-8)
  testthat::expect_equal(celltype_features$reg_dist_enh,
                         expcelltype_features$reg_dist_enh)
  testthat::expect_equal(celltype_features$norm_reg_dist_enh,
                         expcelltype_features$norm_reg_dist_enh,
                         tolerance = 1e-5)
  testthat::expect_equal(celltype_features$reg_dist_enh,
                         expcelltype_features$reg_dist_enh)
  testthat::expect_equal(celltype_features$norm_reg_dist_enh,
                         expcelltype_features$norm_reg_dist_enh,
                         tolerance = 1e-5)



  #### Checking that output is the expected

})

# files <- c("/project/CRUP_scores/CRUP_scores/ENCODE/Single_ended/Cell_lines/GM12878/H3K4me1/H3K4me1.bam",
#            "/project/CRUP_scores/CRUP_scores/ENCODE/Single_ended/Cell_lines/GM12878/H3K4me3/H3K4me3.bam",
#            "/project/CRUP_scores/CRUP_scores/ENCODE/Single_ended/Cell_lines/GM12878/H3K27ac/H3K27ac.bam")
# inputs <- "/project/CRUP_scores/CRUP_scores/ENCODE/Single_ended/Cell_lines/GM12878/Controls/input.bam"
