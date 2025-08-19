test_that("computeCellTypeFeatures functions runs as expected for thyroid", {
  ##Defining inputs for compute features function

  ##testing on thyroid data (we need the input data to be smaller and contained
  ##within the package)
  message("Loading the data from ExperimentHub and the package")
  eh <- ExperimentHub::ExperimentHub()
  # ##we might need to redo this
  files <- c(system.file("extdata",
                         "example/HeLa_H3K4me1.REF_chr19_reduced.bam",
                         package = "CENTRE"),
             system.file("extdata",
                         "example/HeLa_H3K4me3.REF_chr19_reduced.bam",
                         package = "CENTRE"), 
             system.file("extdata",
                         "example/HeLa_H3K4me3.REF_chr19_reduced.bam",
                         package = "CENTRE"))
  inputs <- system.file("extdata",
                        "example/HeLa_input.REF_chr19_reduced.bam",
                        package = "CENTRE")
  metaData <- data.frame(HM = c("H3K4me1", "H3K4me3", "H3K27ac"),
                         condition = c(1, 1, 1), replicate = c(1, 1, 1),
                         bamFile = files, inputFile = rep(inputs, 3))

  tpmpath <- suppressMessages(unname(eh[["EH9545"]]))
  tpmfile <-  read.table(tpmpath, sep = "", stringsAsFactors = F, header = T)
  tpmfile <- tpmfile[grep("E", tpmfile$gene_id), ]
  
  ##selected pairs chr19
  pairs <- readRDS(file = system.file("extdata",
                                         "input_cellType_pairs.rds",
                                         package = "CENTRE"))
  ## are no longer given in output update this object
  expcelltype_features <- readRDS(file = system.file("extdata",
                                                  "expected_cellType_HeLa_reduced.rds",
                                                  package = "CENTRE"))
  ##cell type features with chrom normalization split by chromosome

  message("Check ComputeCellTypeFeatures all pairs")
  celltype_features <- computeCellTypeFeatures(metaData,
                                               replicate = 1,
                                               input.free = FALSE,
                                               cores = 1,
                                               sequencing = "single",
                                               tpmData = tpmfile,
					                                     chr = "chr19",
                                               pairs = pairs)


  celltype_features$pair <- paste(celltype_features$enhancer_id,
                                  celltype_features$gene_id1,
                                  sep = "_")


  testthat::expect_equal(length(unique(celltype_features$pair)),
                         nrow(celltype_features))

  testthat::expect_equal(celltype_features$TPM,
                         expcelltype_features$TPM,
                         tolerance = 1e-8)
  testthat::expect_equal(celltype_features$reg_dist_enh,
                         expcelltype_features$reg_dist_enh)
  testthat::expect_equal(celltype_features$norm_reg_dist_enh,
                         expcelltype_features$norm_reg_dist_enh,
                         tolerance = 1e-2)
  testthat::expect_equal(celltype_features$reg_dist_enh,
                         expcelltype_features$reg_dist_enh)
  testthat::expect_equal(celltype_features$norm_reg_dist_enh,
                         expcelltype_features$norm_reg_dist_enh,
                         tolerance = 1e-2)

  ##only one pair as input
  message("Check ComputeCellTypeFeatures one pair")
  pairs1 <-  pairs[pairs$enhancer_id=="EH38E1958626",]

  celltype_features1 <- computeCellTypeFeatures(metaData,
                                                replicate = 1,
                                                input.free = FALSE,
                                                cores = 1,
                                                sequencing = "single",
                                                tpmData = tpmfile,
                                                chr = "chr19",
                                                pairs = pairs1)

  expcelltype_features1 <- expcelltype_features[expcelltype_features$enhancer_id=="EH38E1958626",]

  testthat::expect_equal(celltype_features1$TPM,
                         expcelltype_features1$TPM,
                         tolerance = 1e-8)
  testthat::expect_equal(celltype_features1$reg_dist_enh,
                         expcelltype_features1$reg_dist_enh)
  testthat::expect_equal(celltype_features1$norm_reg_dist_enh,
                         expcelltype_features1$norm_reg_dist_enh,
                         tolerance = 1e-5)
  testthat::expect_equal(celltype_features1$reg_dist_enh,
                         expcelltype_features1$reg_dist_enh)
  testthat::expect_equal(celltype_features1$norm_reg_dist_enh,
                         expcelltype_features1$norm_reg_dist_enh,
                         tolerance = 1e-5)
  
  
  message("Check ComputeCellTypeFeatures all pairs with 2 cores")
  celltype_features <- computeCellTypeFeatures(metaData,
                                               replicate = 1,
                                               input.free = FALSE,
                                               cores = 2,
                                               sequencing = "single",
                                               tpmData = tpmfile,
                                               chr = "chr19",
                                               pairs = pairs)
  celltype_features$pair <- paste(celltype_features$enhancer_id,
                                  celltype_features$gene_id1,
                                  sep = "_")
  
  
  testthat::expect_equal(length(unique(celltype_features$pair)),
                         nrow(celltype_features))
  
  testthat::expect_equal(celltype_features$TPM,
                         expcelltype_features$TPM,
                         tolerance = 1e-8)
  testthat::expect_equal(celltype_features$reg_dist_enh,
                         expcelltype_features$reg_dist_enh)
  testthat::expect_equal(celltype_features$norm_reg_dist_enh,
                         expcelltype_features$norm_reg_dist_enh,
                         tolerance = 1e-2)
  testthat::expect_equal(celltype_features$reg_dist_enh,
                         expcelltype_features$reg_dist_enh)
  testthat::expect_equal(celltype_features$norm_reg_dist_enh,
                         expcelltype_features$norm_reg_dist_enh,
                         tolerance = 1e-2)
  
  message("Check that error messages are raised")
  testthat::expect_error(computeCellTypeFeatures(metaData,
                          replicate = 1,
                          input.free = FALSE,
                          cores = 2,
                          sequencing = "single",
                          tpmData = tpmfile,
                          chr = "chr19"))
  testthat::expect_error(computeCellTypeFeatures(metaData,
                                                 replicate = 1,
                                                 input.free = FALSE,
                                                 cores = 2,
                                                 sequencing = "single",
                                                 chr = "chr19",
                                                 pairs = pairs))
  colnames(pairs) <- c("gene","enh")
  testthat::expect_error(computeCellTypeFeatures(metaData,
                                                 replicate = 1,
                                                 input.free = FALSE,
                                                 cores = 2,
                                                 sequencing = "single",
                                                 tpmData = tpmfile,
                                                 chr = "chr19",
                                                 pairs = pairs))

})
