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
#'                                                     featuresGeneric = generic_features)
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
                                    featuresGeneric) {
  start_time <- Sys.time()
  ## Pre-eliminary checks and computations


  ## Computing the crup scores
  startPart("Computing CRUP score features")

  ## Calling normalization step only on the chromosomes we have
  chr <- unique(featuresGeneric$chr)
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

  list_enh <- as.data.frame(unique(featuresGeneric$enhancer_id))
  colnames(list_enh) <- c("enhancer_id")

  list_prom <- as.data.frame(unique(featuresGeneric$gene_id2))
  colnames(list_prom) <- c("gene_id2")
  #Get Gencode and CCRes anntotations for the input genes and enhancers
  conn <- RSQLite::dbConnect(RSQLite::SQLite(),
                               system.file("extdata",
                              "Annotation.db",
                              package = "CENTRE"))
  #get chromosome tts new_start and new_end of input genes
  query <- paste("SELECT  gene_id1, chr, transcription_start, new_start, new_end FROM gencode WHERE gene_id1 in (",
		paste0(sprintf("'%s'", list_prom$gene_id2), collapse = ", "),")",sep="" )
  regions_prom <- RSQLite::dbGetQuery(conn, query)
  #get chr middle new_start new_end point of input enhancers
  query_enh <-  paste("SELECT  V5, V1, middle_point, new_start, new_end FROM ccres_enhancer WHERE V5 in (",
			paste0(sprintf("'%s'", list_enh$enhancer_id), collapse = ", "),")",sep="" )
  regions_enhancer <- RSQLite::dbGetQuery(conn, query_enh)
  RSQLite::dbDisconnect(conn)
  ##Making the ranges for the enhancers
  list_enh <- as.data.frame(unique(featuresGeneric$enhancer_id))
  colnames(list_enh) <- c("enhancer_id")

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
