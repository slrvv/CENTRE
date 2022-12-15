################################################################################
# Testing that the pipeline works
################################################################################




##Testing the compute generic features function

test_that("compute features functions runs as expected for thyroid", {
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
  predictions_expected <- read.table("/project/CRUP_scores/EPI/tests/testthat/Thyroid.GTEx-Benchmark.V38.predictions.txt"
                            , sep = "\t", header = T)

  ###Checking that output has correct number of rows
  cat("Check ComputeGenericFeatures")
  file1 <- "/project/CRUP_scores/toSara/Thyroid.GTEx-Benchmark.V38.txt"
  file1 <- read.table(file1, sep = "\t", stringsAsFactors = F, header = T)
  file1 <- file1[,c("V2", "symbol38")]
  file1 <- file1[sample(nrow(file1), 20), ]
  row.names(file1) <- NULL
  generic_features <- computeGenericFeatures(file1)

  generic_features$pair <- paste(generic_features$enhancer_id,
                                 generic_features$gene_id2,
                                 sep = "_")

  expect_equal(length(unique(generic_features$pair)), nrow(generic_features))


  cat("Check ComputeCellTypeFeatures")
  celltype_features <- computeCellTypeFeatures(metaData,
                                    cores = 1,
                                    "single",
                                    tpmfile,
                                    generic_features)
  celltype_features$pair <- paste(celltype_features$enhancer_id,
                                  celltype_features$gene_id2,
                                  sep = "_")
  expect_equal(length(unique(celltype_features$pair)), nrow(celltype_features))

  #### Checking that output is the expected
  cat("Check centrePrediction")
  predictions <- centrePrediction(celltype_features,
                                  generic_features,
                                  model = "/project/CRUP_scores/EPI/model/centre2_final_model.txt" )
  predictions$pair <- paste(predictions$enhancer_id,
                                 predictions$gene_id2,
                                 sep = "_")
  predictions_expected$pair <- paste(predictions_expected$enhancer_id,
                                 predictions_expected$gene_id2,
                                 sep = "_")
  predictions_expected <- predictions_expected[predictions_expected$pair == predictions$pair,]
  expect_equal(predictions$label, predictions_expected$label)
  # expect_equal(predictions$predictions,
  #              predictions_expected$predictions,
  #              tolerance =1e-5)
})

# files <- c("/project/CRUP_scores/CRUP_scores/ENCODE/Single_ended/Cell_lines/GM12878/H3K4me1/H3K4me1.bam",
#            "/project/CRUP_scores/CRUP_scores/ENCODE/Single_ended/Cell_lines/GM12878/H3K4me3/H3K4me3.bam",
#            "/project/CRUP_scores/CRUP_scores/ENCODE/Single_ended/Cell_lines/GM12878/H3K27ac/H3K27ac.bam")
# inputs <- "/project/CRUP_scores/CRUP_scores/ENCODE/Single_ended/Cell_lines/GM12878/Controls/input.bam"
