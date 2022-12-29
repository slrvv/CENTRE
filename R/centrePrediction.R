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
#' #Start by providing genes with their ENSEMBL id
#' candidates <- read.table(system.file("extdata",
#'                         "exampleids.txt", package = "CENTRE"), header = T)
#'                         #Remember to give the columns the name "gene_id"
#'colnames(candidates) <- c("gene_id")
#'#Generate the candidate pairs
#'candidate_pairs <- createPairs(candidates)
#'#Compute the generic features for given names
#'generic_features <- computeGenericFeatures(candidate_pairs)
#'## Prepare the data needed for computing cell type features
#'files <- c(system.file("extdata","HeLa_H3K4me1.bam", package = "CENTRE"),
#'           system.file("extdata","HeLa_H3K4me3.bam", package = "CENTRE"),
#'           system.file("extdata","HeLa_H3K27ac.bam", package = "CENTRE"))
#'
#'inputs <- system.file("extdata", "HeLa_input.bam", package = "CENTRE")
#'
#'metaData <- data.frame(HM = c("H3K4me1", "H3K4me3", "H3K27ac"),
#'                      condition = c(1, 1, 1), replicate = c(1, 1, 1),
#'                      bamFile = files, inputFile = rep(inputs, 3))
#'
#'#More information on this step is found in the crupR documentation
#'
#'tpmfile <- read.table(system.file("extdata", "HeLa.tsv", package = "CENTRE"),
#'                     sep = "", stringsAsFactors = F, header = T)
#'
#'celltype_features <- computeCellTypeFeatures(metaData,
#'                                             cores = 1,
#'                                             "single",
#'                                             tpmfile,
#'                                             generic_features)
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
    model <- system.file("extdata", "centre2_final_model.txt")
  } else {
    check_file(model)
  }
  xgb_model <- xgboost::xgb.load(model)

  ##Transforming data
  feature_matrix <- as.matrix(features_all[, -1])
  test <- xgboost::xgb.DMatrix(data = feature_matrix)
  ##Predicting
  predictions <- predict(xgb_model, test)
  label <- as.numeric(predictions > 0.5)
  #Add the gene and enhancer id's
  predictions <- cbind(features_all[,1], predictions, label)
  return(predictions)
}
