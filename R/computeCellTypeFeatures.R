#' Compute cell type specific features
#'
#' Computes the cell type specific features needed for the CENTRE classification
#' step.
#'
#'
#' @param metaData Dataframe indicating the paths to the ChIP-seq experiments.
#' More information on the format here `crupR::normalize
#' @param replicate: The number of replicates of the ChIP-seq experiments
#' that need to be normalized.
#' @param input.free: Boolean value indicating whether a Control/Input ChIP-seq
#' experiment is provided to go with the Histone Modification ChIP-seq experiments.
#' If the parameter is set to FALSE the normalization of ChIP-seq experiments
#' will be run in input.free mode.
#' @param cores Number of cores to compute the CRUP score features
#' @param sequencing Type of sequencing of the ChIP-seq experiments "paired" or
#' "single". The parameter takes single as default
#' @param tpmfile Dataframe of two columns one with the RNA-seq TPM values,
#' one with the names of the genes given as ENSEMBLE ID's
#' @param chr NULL or a vector of chromosomes. Use only if the crupR
#' normalization should be done for certain chromosomes. If this parameter is
#' not used crupR normalization is done for all chromosomes.
#' Using only certain chromosomes for normalization might change results
#' and is not the intented used of crupR or CENTRE.
#' @param pairs The output of `CENTRE::createPairs()`
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
#' #Create gene enhancer pairs
#' genes <- as.data.frame(c("ENSG00000130203.10",
#' "ENSG00000171119.3"))
#' colnames(genes) <- c("gene_id") #It is important to name the column gene_id
#' pairs <- CENTRE::createPairs(genes)
#'
#' #Compute generic features
#' colnames(pairs) <- c("gene_id", "enhancer_id")
#' generic_features <- CENTRE::computeGenericFeatures(pairs)
#'
#' #Compute Cell-type features
#' files <- c(system.file("extdata/example","HeLa_H3K4me1.REF_chr19.bam", package = "CENTRE"),
#' system.file("extdata/example","HeLa_H3K4me3.REF_chr19.bam", package = "CENTRE"),
#' system.file("extdata/example","HeLa_H3K27ac.REF_chr19.bam", package = "CENTRE"))
#' # Control ChIP-seq experiment to go with the rest of ChIP-seqs
#' inputs <- system.file("extdata/example", "HeLa_input.REF_chr19.bam", package = "CENTRE")
#' metaData <- data.frame(HM = c("H3K4me1", "H3K4me3", "H3K27ac"),
#'                      condition = c(1, 1, 1), replicate = c(1, 1, 1),
#'                       bamFile = files, inputFile = rep(inputs, 3))
#'tpmfile <- read.table(system.file("extdata/example", "HeLa-S3.tsv", package = "CENTRE"),
#'                       sep = "", stringsAsFactors = F, header = T)
#'celltype_features <- CENTRE::computeCellTypeFeatures(metaData,
#'                                                     replicate = 1,
#'                                                     input.free = FALSE,
#'                                                     cores = 1,
#'                                                     sequencing = "single",
#'                                                     tpmfile = tpmfile,
#'                                                     pairs = generic_features)
#'
#'@export
#'@importFrom crupR normalize getEnhancers
#'@import utils
#'@importFrom GenomicRanges GRanges findOverlaps elementMetadata
#'@importFrom IRanges IRanges
#'@importFrom stats reshape

computeCellTypeFeatures <- function(metaData,
                                    replicate,
                                    input.free = FALSE,
                                    cores,
                                    sequencing = "single",
                                    tpmfile,
                                    chr = NULL,
                                    pairs) {
  startTime <- Sys.time()
  ## Computing the crup scores
  startPart("Computing CRUP score features")
  ## Calling normalization step only on the chromosomes we have
  normalized <- crupR::normalize(metaData = metaData,
                                 condition = 1,
                                 replicate = replicate,
                                 mapq = 10,
                                 input.free = input.free,
                                 genome = "hg38",
                                 sequencing = sequencing,
                                 chroms = chr,
                                 C = cores)
  #Get CRUP enhancer probabilities
  crupScores <- crupR::getEnhancers(data = normalized, C = cores, all = TRUE)
  crupScores <- crupScores$D
  ## check what parts of this are necessary
  listEnh <- as.data.frame(unique(pairs$enhancer_id))
  colnames(listEnh) <- c("enhancer_id")
  listProm <- as.data.frame(unique(pairs$gene_id2))
  colnames(listProm) <- c("gene_id2")
  #Get Gencode and CCRes anntotations for the input genes and enhancers
  regions <- createRegionsDf(listProm, listEnh, pairs)
  pairs$pair <- paste(pairs$enhancer_id, pairs$gene_id2, sep = "_")

  cat("Getting the CRUP-EP scores for enhancer, promoter and the regulatory
      distance\n")
  #Crup enhancer scores for enhancer
  crupEPenh <- compute_crup_enhancer(regions,
                                     crupScores)
  crupFeatures <- merge(pairs,
                        crupEPenh,
                        by.x = "enhancer_id",
                        by.y = "enhancer_id",
                        all.x = TRUE)

  #CRUP enhancer scores for promoter
  crupEPprom <- compute_crup_promoter(regions,
                                      crupScores)
  crupFeatures <- merge(crupFeatures,
                        crupEPprom,
                        by.x = "gene_id2",
                        by.y = "gene_id2",
                        all.x = TRUE)

  #create the between_ranges objects that is used for the distance calculations
  betweenRanges <- createBetweenRanges(regions)
  crupFeatures <- compute_crup_reg_distance_enh(crupFeatures,
                                                crupScores,
                                                betweenRanges)
  ##Get CRUP promoter probabilities

  # Compute the promoter probability from probA and probE
  # In CRUP probA is the probability of a region being an active reg. element
  # probE is the probability of a region being an active enhancer
  crupScores$probP <- crupScores$probA *(1 - crupScores$probE)

  cat("Getting the CRUP-PP scores for enhancer")

  #Crup promoter scores for enhancer
  crupPPenh <- compute_crup_enhancer(regions,
                                       crupScores,
                                       promprob = TRUE)

  crupFeatures <- merge(crupFeatures,
                        crupPPenh,
                        by.x = "enhancer_id",
                        by.y = "enhancer_id",
                        all.x = TRUE)

  #Crup promoter scores for promoter
  crupPPprom <- compute_crup_promoter(regions,
                                      crupScores,
                                      promprob = TRUE)
  crupFeatures <- merge(crupFeatures,
                        crupPPprom,
                        by.x = "gene_id2",
                        by.y = "gene_id2",
                        all.x = TRUE)

  #Crup promoter scores for distance
  crupFeatures <- compute_crup_reg_distance_prom(crupFeatures,
                                                 crupScores,
                                                 betweenRanges)

  endPart()

  startPart("Getting the TPM values")
  features_table_all <- get_rnaseq(crupFeatures, tpmfile)

  features_table_all[is.na(features_table_all)] <- 0

  features_table_all <- features_table_all[, c("gene_id2",
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
                                               "TPM")]

  cat(paste0('time: ', format(Sys.time() - startTime), "\n"))
  endPart()
  return(features_table_all)

}
