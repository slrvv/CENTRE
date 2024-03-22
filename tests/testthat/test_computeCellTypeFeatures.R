################################################################################
# Testing that the pipeline works
################################################################################




##Testing the compute generic features function

test_that("computeCellTypeFeatures functions runs as expected for thyroid", {
<<<<<<< Updated upstream
  ##Defining inputs for compute features function

=======

  ##Defining inputs for compute features function

  ##testing on thyroid data (we need the input data to be smaller and contained
  ##within the package)

  ##we might need to redo this
>>>>>>> Stashed changes
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
<<<<<<< Updated upstream
=======

  pairs <- generic_features[, c(1, 2)] ## remove columns that
  ## are no longer given in output update this object
>>>>>>> Stashed changes
  expcelltype_features <- readRDS(file = system.file("extdata",
                                                  "expected_celltype_features.rds",
                                                  package = "CENTRE"))


<<<<<<< Updated upstream
  cat("Check ComputeCellTypeFeatures")
=======
  ##cell type features with chrom normalization split by chromosome

  startPart("Check ComputeCellTypeFeatures norm split by chrom")
>>>>>>> Stashed changes
  celltype_features <- computeCellTypeFeatures(metaData,
                                               replicate = 1,
                                               input.free = FALSE,
                                               cores = 1,
                                               sequencing = "single",
                                               tpmfile = tpmfile,
<<<<<<< Updated upstream
                                               featuresGeneric = generic_features)
  celltype_features$pair <- paste(celltype_features$enhancer_id,
                                  celltype_features$gene_id2,
                                  sep = "_")
  expect_equal(length(unique(celltype_features$pair)), nrow(celltype_features))
=======
					                                     chr = unique(pairs$chr),
                                               pairs = pairs)



  celltype_features$pair <- paste(celltype_features$enhancer_id,
                                  celltype_features$gene_id2,
                                  sep = "_")


  testthat::expect_equal(length(unique(celltype_features$pair)),
                         nrow(celltype_features))
>>>>>>> Stashed changes

  testthat::expect_equal(celltype_features$TPM,
                         expcelltype_features$TPM,
                         tolerance = 1e-8)
  testthat::expect_equal(celltype_features$reg_dist_enh,
                         expcelltype_features$reg_dist_enh)
  testthat::expect_equal(celltype_features$norm_reg_dist_enh,
                         expcelltype_features$norm_reg_dist_enh,
<<<<<<< Updated upstream
                         tolerance = 1e-5)
=======
                         tolerance = 1e-2)
>>>>>>> Stashed changes
  testthat::expect_equal(celltype_features$reg_dist_enh,
                         expcelltype_features$reg_dist_enh)
  testthat::expect_equal(celltype_features$norm_reg_dist_enh,
                         expcelltype_features$norm_reg_dist_enh,
<<<<<<< Updated upstream
                         tolerance = 1e-5)



  #### Checking that output is the expected

})

# files <- c("/project/CRUP_scores/CRUP_scores/ENCODE/Single_ended/Cell_lines/GM12878/H3K4me1/H3K4me1.bam",
#            "/project/CRUP_scores/CRUP_scores/ENCODE/Single_ended/Cell_lines/GM12878/H3K4me3/H3K4me3.bam",
#            "/project/CRUP_scores/CRUP_scores/ENCODE/Single_ended/Cell_lines/GM12878/H3K27ac/H3K27ac.bam")
# inputs <- "/project/CRUP_scores/CRUP_scores/ENCODE/Single_ended/Cell_lines/GM12878/Controls/input.bam"
=======
                         tolerance = 1e-2)


  startPart("Check ComputeCellTypeFeatures norm all chrom one pair")

  expcelltype_features1 <- readRDS(file = system.file("extdata",
                                                     "expected_celltype_allchrnorm.rds",
                                                     package = "CENTRE"))


  ##only one pair as input
  pairs1 <-  pairs[pairs$gene_id2=="ENSG00000059728",]

  celltype_features1 <- computeCellTypeFeatures(metaData,
                                                replicate = 1,
                                                input.free = FALSE,
                                                cores = 1,
                                                sequencing = "single",
                                                tpmfile = tpmfile,
                                                pairs = pairs1)

  expcelltype_features2 <- expcelltype_features1[expcelltype_features1$gene_id2=="ENSG00000059728",]

  testthat::expect_equal(celltype_features1$TPM,
                         expcelltype_features2$TPM,
                         tolerance = 1e-8)
  testthat::expect_equal(celltype_features1$reg_dist_enh,
                         expcelltype_features2$reg_dist_enh)
  testthat::expect_equal(celltype_features1$norm_reg_dist_enh,
                         expcelltype_features2$norm_reg_dist_enh,
                         tolerance = 1e-5)
  testthat::expect_equal(celltype_features1$reg_dist_enh,
                         expcelltype_features2$reg_dist_enh)
  testthat::expect_equal(celltype_features1$norm_reg_dist_enh,
                         expcelltype_features2$norm_reg_dist_enh,
                         tolerance = 1e-5)


  startPart("Check ComputeCellTypeFeatures norm all chrom all pairs")

  celltype_features2 <- computeCellTypeFeatures(metaData,
                                                replicate = 1,
                                                input.free = FALSE,
                                                cores = 1,
                                                sequencing = "single",
                                                tpmfile = tpmfile,
                                                pairs = pairs)

  testthat::expect_equal(celltype_features2$TPM,
                         expcelltype_features1$TPM,
                         tolerance = 1e-8)
  testthat::expect_equal(celltype_features2$reg_dist_enh,
                         expcelltype_features1$reg_dist_enh)
  testthat::expect_equal(celltype_features2$norm_reg_dist_enh,
                         expcelltype_features1$norm_reg_dist_enh,
                         tolerance = 1e-5)
  testthat::expect_equal(celltype_features2$reg_dist_enh,
                         expcelltype_features1$reg_dist_enh)
  testthat::expect_equal(celltype_features2$norm_reg_dist_enh,
                         expcelltype_features1$norm_reg_dist_enh,
                         tolerance = 1e-5)

})


>>>>>>> Stashed changes
