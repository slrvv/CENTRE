#' Compute generic features
#'
#' Computes the generic features needed by CENTRE to make predictions
#'
#' @param x Dataframe with the gene-enhancer pairs of interest
#'
#' @return
#' A table containting the following computed features :
#'* distance: Distance between gene and enhancer
#'* combined_tests: Combined value of the Wilcoxon tests (CAGE, DNase
#'expression, CRUP expression and DNase DNase)
#'* crup_COR: CRUP correlation scores
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
#' @export
#' @import utils
#' @importFrom metapod combineParallelPValues
#' @importFrom RSQLite dbConnect dbGetQuery dbDisconnect
computeGenericFeatures <- function(x) {
  start_time <- Sys.time()
  ## Pre-eliminary checks and computations

  x$gene_id1 <- gsub("\\..*", "", x[, 1])

  ## Computing the distance features

  startPart("Computing distance features")


  colnames(x) <- c("gene_id", "enhancer_id", "gene_id2")
  featuresDistances <- computeDistances(x)

  cat("Removing pairs with distance over 500 Kb")
  featuresDistances <- featuresDistances[abs(featuresDistances$distance)
                                           <= 500000, ]
  endPart()

  ## Getting the values for the Wilcoxon tests and the CRUP correlations
  startPart("Get Wilcoxon tests and CRUP correlations")

  featuresDistances$pair <- paste(featuresDistances$enhancer_id,
                                   featuresDistances$gene_id2,
                                   sep = "_")

  conn <- RSQLite::dbConnect(RSQLite::SQLite(),
                             system.file("extdata",
                                         "PrecomputedData.db",
                                         package = "CENTRE"))


  cageDf <- getPrecomputedValues("cage_test_data",
                                  "wilcoxtest_cage",
                                  featuresDistances,
                                  conn)
  dhsexpDf <- getPrecomputedValues("dhsexp_test_data",
                                    "wilcoxtest_dhs_exp",
                                    featuresDistances,
                                    conn)
  crupexpDf <- getPrecomputedValues("crupexp_test_data",
                                     "wilcoxtest_crup_exp",
                                     featuresDistances,
                                     conn)
  dhsdhsDf <- getPrecomputedValues("dhsdhs_test_data",
                                    "wilcoxtest_dhs_dhs",
                                    featuresDistances,
                                    conn)
  crupCorDf <- getPrecomputedValues("crup_cor",
                                      "cor_CRUP",
                                      featuresDistances,
                                      conn)

  RSQLite::dbDisconnect(conn)

  cageWilcoxTest <- cageDf[featuresDistances$pair, 1]
  dhsexpWilcoxTest <- dhsexpDf[featuresDistances$pair, 1]
  crupexpWilcoxTest <- crupexpDf[featuresDistances$pair, 1]
  dhsdhsWilcoxTest <- dhsdhsDf[featuresDistances$pair, 1]
  crupCor <- crupCorDf[featuresDistances$pair, 1]

  pValue <- list(cageWilcoxTest, dhsexpWilcoxTest,
                 crupexpWilcoxTest, dhsdhsWilcoxTest)
  pValue[is.na(pValue)] <- 0
  ## Combining the values of the Wilcoxon tests
  combined_tests <- metapod::combineParallelPValues(pValue, method="fisher")


  featuresDistances$combined_tests <- -log(combined_tests$p.value)
  featuresDistances$crup_cor <- crupCor
  ## Return the table of features
  featuresGeneric <- featuresDistances[, c("gene_id2", "enhancer_id", "chr",
                                         "middle_point", "transcription_start",
                                         "distance", "crup_cor",
                                         "combined_tests")]

  featuresGeneric[is.na(featuresGeneric)] <- 0
  cat(paste0('time: ', format(Sys.time() - start_time), "\n"))
  return(featuresGeneric)

}
