#' Compute cell type specific features
#'
#' Computes the cell type specific features needed for the CENTRE classification
#' step.
#'
#'
#' @param metaData Dataframe indicating the paths to the ChIP-seq experiments.
#' More information on the format here `crupR::normalize()`
#' @param cores Number of cores to compute the CRUP score features
#' @param sequencing Type of sequencing of the ChIP-seq experiments "paired" or
#' "single". The parameter takes single as default
#' @param tpmfile Dataframe of two columns one with the RNA-seq TPM values,
#' one with the names of the genes given as ENSEMBLE ID's
#' @param featuresGeneric The output of `CENTRE::computeGenericFeatures()`
#'
#'
#' @return
#' A table containting the following computed features :
#'* CRUP enhancer score for enhancer region, promoter region and the region
#'between the enhancer and the promoter
#'* CRUP promoter score for enhancer region, promoter region and the region
#'between the enhancer and the promoter
#'* RNA-seq TPM values
#'
#'
#' @examples
#' candidates <- read.table(system.file("extdata",
#' "exampleids.txt", package = "CENTRE"), header = T)
#'
#' #Remember to give the columns the name "gene_id"
#' colnames(candidates) <- c("gene_id")
#'
#' #Generate the candidate pairs
#' candidate_pairs <- createPairs(candidates)
#'
#' #Compute the generic features for given names
#' generic_features <- computeGenericFeatures(candidate_pairs)
#'
#' ## Prepare the data needed for computing cell type
#' featuresfiles <- c(system.file("extdata","HeLa_H3K4me1.bam", package = "CENTRE"),
#'                   system.file("extdata","HeLa_H3K4me3.bam", package = "CENTRE"),
#'                   system.file("extdata","HeLa_H3K27ac.bam", package = "CENTRE"))
#'
#'inputs <- system.file("extdata", "HeLa_input.bam", package = "CENTRE")
#' metaData <- data.frame(HM = c("H3K4me1", "H3K4me3", "H3K27ac"),
#'                        condition = c(1, 1, 1), replicate = c(1, 1, 1),
#'                        bamFile = files, inputFile = rep(inputs, 3))
#'#More information on this step is found in the crupR documentation
#'tpmfile <- read.table(system.file("extdata", "HeLa.tsv", package = "CENTRE"),
#'                      sep = "", stringsAsFactors = F, header = T)
#'
#'celltype_features <- computeCellTypeFeatures(metaData,
#'                                            cores = 1,
#'                                            "single",
#'                                            tpmfile,
#'                                            generic_features)
#'@export
#'@importFrom crupR normalize getEnhancers
#'@import utils
#'@importFrom GenomicRanges GRanges findOverlaps elementMetadata
#'@importFrom IRanges IRanges
#'@importFrom stats reshape

computeCellTypeFeatures <- function(metaData,
                                    condition,
                                    replicate,
                                    mapq = 10,
                                    input.free = FALSE,
                                    cores,
                                    sequencing = "single",
                                    tpmfile,
                                    featuresGeneric) {
  start_time <- Sys.time()
  ## Pre-eliminary checks and computations


  ## Computing the crup scores
  startPart("Computing CRUP score features")

  ## Calling normalization step only on the chromosomes we have
  chr <- unique(featuresGeneric$chr)
  normalized <- crupR::normalize(metaData = metaData,
                                 condition = condition,
                                 replicate = replicate,
                                 mapq = mapq,
                                 input.free = input.free,
                                 genome = "hg38",
                                 sequencing = sequencing,
                                 chroms = chr,
                                 C = cores)
  #Get CRUP enhancer probabilities
  crupScores <- crupR::getEnhancers(data = normalized, C = cores, all = TRUE)
  crupScores <- crupScores$D
  ##Making the ranges for the enhancers
  list_enh <- as.data.frame(unique(featuresGeneric$enhancer_id))
  colnames(list_enh) <- c("enhancer_id")

  regions_enhancer <-  merge(list_enh,
                    ccres_enhancer[, c("V1", "V5", "new_start", "new_end")],
                    by.x = "enhancer_id",
                    by.y = "V5")


  ##Making the ranges for the genes

  list_prom <- as.data.frame(unique(featuresGeneric$gene_id2))
  colnames(list_prom) <- c("gene_id2")
  regions_prom <- merge(list_prom,
                   gencode[, c("chr", "gene_id1", "new_start", "new_end")],
                   by.x = "gene_id2",
                   by.y = "gene_id1")


  cat("Getting the CRUP-EP scores for enhancer, promoter and the regulatory
      distance")

  #Crup enhancer scores for enhancer
  crup_EP_enh <- compute_crup_enhancer(regions_enhancer,
                                       list_enh,
                                       crupScores)

  crup_features <- merge(featuresGeneric,
                        crup_EP_enh,
                        by.x = "enhancer_id",
                        by.y = "cres_name",
                        all.x = TRUE)
  #CRUP enhancer scores for promoter
  crup_EP_prom <- compute_crup_promoter(regions_prom,
                                        list_prom,
                                        crupScores)

  crup_features <- merge(crup_features,
                        crup_EP_prom,
                        by.x = "gene_id2",
                        by.y = "gene_name",
                        all.x = TRUE)
  ##crup enhancer scores for distance
  crup_features <- compute_crup_reg_distance_enh(crup_features, crupScores)


  ##Get CRUP promoter probabilities

  ### Compute the promoter probability from probA and probE
  crupScores$probP <- crupScores$probA *(1 - crupScores$probE)

  cat("Getting the CRUP-PP scores for enhancer")

  #Crup promoter scores for enhancer

  crup_PP_enh <- compute_crup_enhancer(regions_enhancer,
                                       list_enh,
                                       crupScores,
                                       promprob = T)

  crup_features <- merge(crup_features,
                        crup_PP_enh,
                        by.x = "enhancer_id",
                        by.y = "cres_name",
                        all.x = TRUE)

  #Crup promoter scores for promoter

  crup_PP_prom <- compute_crup_promoter(regions_prom,
                                        list_prom,
                                        crupScores,
                                        promprob = T)
  crup_features <- merge(crup_features,
                        crup_PP_prom,
                        by.x = "gene_id2",
                        by.y = "gene_name",
                        all.x = TRUE)

  #Crup promoter scores for distance
  crup_features <- compute_crup_reg_distance_prom(crup_features,
                                             crupScores)


  endPart()

  startPart("Getting the TPM values")
  features_table_all <- get_rnaseq(crup_features, tpmfile)


  ###Some renaming and so on
  features_table_all[is.na(features_table_all)] <- 0
  features_table_all <- features_table_all[, c(1, 2, 9, 10, 11, 12, 13, 14, 15,
                                              16, 17, 18, 22, 23, 24, 25,
                                              26, 27, 28, 29, 30, 31, 32, 33,
                                              34, 35, 36)]


  colnames(features_table_all) <- c("gene_id2",
                                    "enhancer_id",
                                    "EP_prob_enh.1",
                                    "EP_prob_enh.2",
                                    "EP_prob_enh.3",
                                    "EP_prob_enh.4",
                                    "EP_prob_enh.5",
                                    "EP_prob_gene.1",
                                    "EP_prob_gene.2",
                                    "EP_prob_gene.3",
                                    "EP_prob_gene.4",
                                    "EP_prob_gene.5",
                                    "reg_dist_enh",
                                    "norm_reg_dist_enh",
                                    "PP_prob_enh.1",
                                    "PP_prob_enh.2",
                                    "PP_prob_enh.3",
                                    "PP_prob_enh.4",
                                    "PP_prob_enh.5",
                                    "PP_prob_gene.1",
                                    "PP_prob_gene.2",
                                    "PP_prob_gene.3",
                                    "PP_prob_gene.4",
                                    "PP_prob_gene.5",
                                    "reg_dist_prom",
                                    "norm_reg_dist_prom",
                                    "RNA_seq")

  cat(paste0('time: ', format(Sys.time() - start_time), "\n"))
  endPart()
  return(features_table_all)

}
