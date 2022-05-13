#' Compute features
#'
#' Computes the features needed for the classification step.
#'
#' @param file Path to a file with either genes and enhancer pairs or with only
#' genes, separated by tabs. The genes need to be given as gencode IDs and
#' the enhancers as ENCODE CCREs IDs.
#' @param metaData Dataframe indicating the paths to the ChIP-seq experiments.
#' More information on the format here `crupR::normalize()`
#' @param cores Number of cores to compute the CRUP score features
#' @param tpmpath Path to a file with the RNA-seq TPM values, with the names of
#' the genes given as ENSEMBLE ID's
#'
#' @return
#' A table containting the following computed features :
#'* Distance between gene and enhancer
#'* CRUP enhancer score for enhancer region, promoter region and the region between the enhancer and the promoter
#'* CRUP promoter score for enhancer region, promoter region and the region between the enhancer and the promoter
#'* CAGE wilcoxon test score
#'* DNAase-seq wilcoxon test score
#'* CRUP correlation scores
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
#'inputs <- "/project/CRUP_scores/CRUP_scores/ENCODE/Single_ended/Tissues/body_pancreas/Controls/input.bam"
#'
#'metaData <- data.frame(HM = c("H3K4me1","H3K4me3","H3K27ac"),
#'                       condition = c(1,1,1), replicate = c(1,1,1),
#'                       bamFile = files, inputFile = rep(inputs,3))
#'file <- "/project/CRUP_scores/EPI/example1.txt"
#'compute_features(file, metaData, cores)
#'
#'
compute_features <- function(file, metaData, cores, tpmpath){

  startPart("Computing features")


  ##Checking if the file exists

  check_file(file)
  x <- read.table(file, sep = "\t", stringsAsFactors = F)
  check_file(tpmpath)
  tpmfile <- read.table(tpmpath, sep = "", stringsAsFactors = F, header = T)
  print(head(tpmfile))
  cat("Computing distance features")

  if (ncol(x) == 2) {
    cat("Input case with enhancer gene pairs")
    features_distances <- distances_gene_enhancer(x)


  }
  else{
    cat("Input case with genes")
    cat("calling map genes to chromosomes")
    returned_mat <- map_genes_to_chromosomes(x)
    features_distances <- distances_gene_enhancer(returned_mat)

  }


  cat("Computing CRUP score features")
  chr <- unique(features_distances$chr_gene)
  features_crup <- compute_features_crup_scores(x, metaData, cores, chr)

  cat("Computing Wilcoxon tests and CRUP correlations")
  x$gene_id2 <- gsub("\\..*","",x[,1])

  x$pair <- paste( x[,2],x$gene_id2, sep= "_")
  tpmvalue <- get_rnaseq(x, tpmfile)
  crup_cor <- crup_correlations(x)
  wilcoxon_features <- wilcoxon_test_features(x)
  print(wilcoxon_features)

  combined_tests<-scran::combinePValues(wilcoxon_features$cage_wilcoxon_test,
                                        wilcoxon_features$dhsexp_wilcoxon_test,
                                        wilcoxon_features$crupexp_wilcoxon_test,
                                        wilcoxon_features$dhsdhs_wilcoxon_test,
                                        method="fisher")
  combined_tests<-log(combined_tests)
  print(combined_tests)

  features_table <- cbind(
    features_distances[c("V1", "V2", "chr_gene", "distance")],
    features_crup[,3:ncol(features_crup)],
    combined_tests,
    crup_cor$crup_correlations,
    tpmvalue)

  features_table[is.na(features_table)] <- 0
  endPart()
  return(features_table)

}
