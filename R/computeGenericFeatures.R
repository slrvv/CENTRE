#' Compute generic features
#'
#' Computes the generic features needed by CENTRE to make predictions
#'
#' @param file Path to a file with either genes and enhancer pairs.
#' @param metaData Dataframe indicating the paths to the ChIP-seq experiments.
#' More information on the format here `crupR::normalize()`
#'
#'
#' @return
#' A table containting the following computed features :
#'* distance: Distance between gene and enhancer
#'* combined_tests: Combined value of the Wilcoxon tests (CAGE, DNase
#'expression, CRUP expression and DNase DNase)
#'* crup_COR: CRUP correlation scores
#'
#' @examples
#'
#'
#'
computeGenericFeatures <- function(x) {
  ## Pre-eliminary checks and computations

  x$gene_id1 <- gsub("\\..*", "", x[, 1])

  ## Computing the distance features

  startPart("Computing distance features")


  colnames(x) <- c("gene_id", "enhancer_id", "gene_id2")
  features_distances <- distances_gene_enhancer(x)

  cat("Removing pairs with distance over 500 Kb")
  features_distances <- features_distances[abs(features_distances$distance)
                                           <= 500000, ]
  endPart()

  ## Getting the values for the Wilcoxon tests and the CRUP correlations
  startPart("Get Wilcoxon tests and CRUP correlations")

  features_distances$pair <- paste(features_distances$enhancer_id,
                                   features_distances$gene_id2,
                                   sep = "_")

  wilcoxon_features <- wilcoxon_test_crup_cor(features_distances)

  ## Combining the values of the Wilcoxon tests
  wilcoxon_features$combined_tests <- scran::combinePValues(
    wilcoxon_features$cage_wilcoxon_test,
    wilcoxon_features$dhsexp_wilcoxon_test,
    wilcoxon_features$crupexp_wilcoxon_test,
    wilcoxon_features$dhsdhs_wilcoxon_test,
    method="fisher")


  wilcoxon_features$combined_tests <-log(wilcoxon_features$combined_tests)

  ## Return the table of features
  features_table <- wilcoxon_features[, c("gene_id2", "enhancer_id", "chr",
                                         "middle_point", "transcription_start",
                                         "distance", "cor_CRUP",
                                         "combined_tests")]

  features_table[is.na(features_table)] <- 0

  return(features_table)

}
