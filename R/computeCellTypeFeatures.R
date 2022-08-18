#' Compute cell type specific features
#'
#' Computes the cell type specific features needed for the CENTRE classification
#' step.
#'
#' @param file Path to a file with gene and enhancer pairs.
#' @param metaData Dataframe indicating the paths to the ChIP-seq experiments.
#' More information on the format here `crupR::normalize()`
#' @param cores Number of cores to compute the CRUP score features
#' @param tpmpath Path to a file with the RNA-seq TPM values, with the names of
#' the genes given as ENSEMBLE ID's
#'
#' @return
#' A table containting the following computed features :
#'* CRUP enhancer score for enhancer region, promoter region and the region between the enhancer and the promoter
#'* CRUP promoter score for enhancer region, promoter region and the region between the enhancer and the promoter
#'* RNA-seq TPM values
#'
#'
#' @export
#'
#' @examples
#'files <- c(
#'"/project/CRUP_scores/CRUP_scores/ENCODE/Single_ended/Tissues/body_pancreas/H3K4me1/H3K4me1.bam",
#'"/project/CRUP_scores/CRUP_scores/ENCODE/Single_ended/Tissues/body_pancreas/H3K4me3/H3K4me3.bam",
#'"/project/CRUP_scores/CRUP_scores/ENCODE/Single_ended/Tissues/body_pancreas/H3K27ac/H3K27ac.bam" )
#'
#'features_generics <- "/project/CRUP_scores/CRUP_scores/ENCODE/Single_ended/Tissues/body_pancreas/Controls/features_generic.bam"
#'
#'metaData <- data.frame(HM = c("H3K4me1","H3K4me3","H3K27ac"),
#'                       condition = c(1,1,1), replicate = c(1,1,1),
#'                       bamFile = files, features_genericFile = rep(features_generics,3))
#'file <- "/project/CRUP_scores/EPI/example1.txt"
#'compute_features(file, metaData, cores)
#'
#'
computeCellTypeFeatures <- function(metaData, cores, sequencing, tpmfile,
                                    features_generic){
  ## Pre-eliminary checks and computations


  ## Computing the crup scores
  startPart("Computing CRUP score features")

  ## Calling normalization step only on the chromosomes we have
  chr <- unique(features_generic$chr)
  normalized <- crupR::normalize(metaData = metaData,
                                 condition = 1,
                                 replicate = 1,
                                 genome = "hg38",
                                 sequencing = "single",
                                 chroms = chr,
                                 cores = cores)
  #Get CRUP enhancer probabilities
  crup_scores_enh <- crupR::getEnhancers(data = normalized, cores = cores)
  crup_scores_enh <- crup_scores_enh$data_matrix


  ##Making the ranges for the enhancers
  list_enh <- as.data.frame(unique(features_generic$enhancer_id))
  colnames(list_enh) <- c("enhancer_id")

  regions_enhancer <-  merge(list_enh,
                    ccres_enhancer[,c('V1', 'V5', 'new_start', 'new_end')],
                    by.x='enhancer_id',
                    by.y='V5')


  ##Making the ranges for the genes

  list_prom <- as.data.frame(unique(features_generic$gene_id2))
  colnames(list_prom) <- c("gene_id2")
  regions_prom <- merge(list_prom,
                   gencode[,c('chr','gene_id1','new_start', 'new_end')],
                   by.x='gene_id2',
                   by.y='gene_id1')


  cat("Getting the CRUP-EP scores for enhancer, promoter and the regulatory
      distance")

  #Crup enhancer scores for enhancer
  crup_EP_enh <- compute_crup_enhancer(regions_enhancer, list_enh, crup_scores_enh)

  crup_features <-merge(features_generic,
                        crup_EP_enh,
                        by.x="enhancer_id",
                        by.y="cres_name",
                        all.x=TRUE)

  #CRUP enhancer scores for promoter
  crup_EP_prom <- compute_crup_promoter(regions_prom,list_prom, crup_scores_enh)

  crup_features <-merge(crup_features,
                        crup_EP_prom,
                        by.x="gene_id2",
                        by.y="gene_name",
                        all.x=TRUE)
  ##crup enhancer scores for distance
  crup_features <- compute_crup_reg_distance(crup_features, crup_scores_enh)


  ##Get CRUP promoter probabilities
  crup_scores_prom <- crupR::getEnhancers(data = normalized, cores = cores, promprob = T)
  crup_scores_prom <- crup_scores_prom$data_matrix

  cat("Getting the CRUP-PP scores for enhancer")

  #Crup promoter scores for enhancer

  crup_PP_enh <- compute_crup_enhancer(regions_enhancer, list_enh, crup_scores_prom, promprob = T)

  crup_features <-merge(crup_features,
                        crup_PP_enh,
                        by.x="enhancer_id",
                        by.y="cres_name",
                        all.x=TRUE)

  #Crup promoter scores for promoter

  crup_PP_prom <- compute_crup_promoter(regions_prom, list_prom, crup_scores_prom, promprob = T)
  crup_features <-merge(crup_features,
                        crup_PP_prom,
                        by.x="gene_id2",
                        by.y="gene_name",
                        all.x=TRUE)

  #Crup promoter scores for distance
  crup_features <- compute_crup_reg_distance(crup_features, crup_scores_prom, prom = T)


  endPart()

  startPart("Getting the TPM values")
  features_table_all <- get_rnaseq(crup_features, tpmfile)


  ###Some renaming and so on
  features_table_all[is.na(features_table_all)] <- 0
  print(head(features_table_all))
  features_table_all <- features_table_all[,c(1, 2,10,11, 12, 13, 14,15, 16, 17,
                                              18, 19,22, 23,24, 25, 26, 27, 28,
                                              29, 30, 31, 32, 33, 34, 35, 36)]

  print(head(features_table_all))
  colnames(features_table_all) <- c('gene_id2',
                                    'enhancer_id',
                                    'EP_prob_enh.1',
                                    'EP_prob_enh.2',
                                    'EP_prob_enh.3',
                                    'EP_prob_enh.4',
                                    'EP_prob_enh.5',
                                    'EP_prob_gene.1',
                                    'EP_prob_gene.2',
                                    'EP_prob_gene.3',
                                    'EP_prob_gene.4',
                                    'EP_prob_gene.5',
                                    'reg_dist_enh',
                                    'norm_reg_dist_enh',
                                    'PP_prob_enh.1',
                                    'PP_prob_enh.2',
                                    'PP_prob_enh.3',
                                    'PP_prob_enh.4',
                                    'PP_prob_enh.5',
                                    'PP_prob_gene.1',
                                    'PP_prob_gene.2',
                                    'PP_prob_gene.3',
                                    'PP_prob_gene.4',
                                    'PP_prob_gene.5',
                                    'reg_dist_prom',
                                    'norm_reg_dist_prom',
                                    'RNA_seq')

  print(head(features_table_all))
  endPart()
  return(features_table_all)

}
