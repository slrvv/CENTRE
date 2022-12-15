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
#' @export
#'
#' @examples
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
  if (is.NULL(model)){
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
