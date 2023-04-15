#' CENTRE predicton
#'
#' @param features_celltype The cell type specific features returned by
#' `CENTRE::computeCellTypeFeatures()`
#' @param features_generic The generic features returned by the function
#' `CENTRE::computeGenericFeatures()`
#' @param model Path to the model the predictions will be computed on,
#' The default is the CENTRE model.
#'
#' @return Dataframe containing the enhancer gene pairs and the probability of
#' them interacting based on CENTRE model
#'
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

#'# Finally compute the predictions
#'predictions <- centrePrediction(celltype_features,
#'generic_features)
#' @export
#' @importFrom stats predict
#' @import utils
#' @importFrom xgboost xgb.load xgb.DMatrix
#'
centrePrediction <- function(features_celltype,features_generic, model = NULL){
  #Merge the generic features and the cell type features

  start_time <- Sys.time()
  features_generic$distance <- abs(features_generic$distance) # make distance absolute distance
  #generate the pair id to merge both feature sets
  features_generic$pair <- paste(features_generic$enhancer_id,
                                 features_generic$gene_id2,
                                 sep = "_")
  features_celltype$pair <- paste(features_celltype$enhancer_id,
                                  features_celltype$gene_id2,
                                  sep = "_")

  #remove all non-feature columns except pair id
  features_generic <- features_generic[, c(6, 7, 8, 9)]
  features_celltype <- features_celltype[, - c(1, 2)]
  #mergeboth datasets
  features_all <- merge(features_celltype,
                        features_generic, by.x = "pair", by.y = "pair")

  ##Loading the xgboost model
  if (is.null(model)){
    model <- system.file("extdata", "centre2_final_model.txt", package = "CENTRE")
  } else {
    check_file(model)
  }
  xgb_model <- xgboost::xgb.load(model)


  ##Transforming data
  colnames(features_all) <- NULL
  pairs <- features_all[, 1]
  features_all <- features_all[, -1]


  feature_matrix <- data.matrix(features_all)

  test <- xgboost::xgb.DMatrix(data = feature_matrix)
  ##Predicting
  predictions <- predict(xgb_model, test)
  label <- as.numeric(predictions > 0.5)
  #Add the gene and enhancer id's
  predictions <- cbind(pairs, predictions, label)
  predictions <- as.data.frame(predictions)
  cat(paste0('time: ', format(Sys.time() - start_time), "\n"))
  return(predictions)
}
