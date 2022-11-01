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
centrePrediction <- function(features_celltype,features_generic, model){
  ##Merge the generic features and the cell type features
  features_generic$distance <- abs(features_generic$distance)
  features_all <- cbind(features_celltype, features_generic[,c(6,7,8)])

  ##change it to a merge later

  ##Loading the xgboost model
  xgb_model <- xgboost::xgb.load(model)

  ##Transforming data
  feature_matrix <- as.matrix(features_all[, c(-1, -2)])
  test <- xgboost::xgb.DMatrix(data = feature_matrix)

  ##Predicting
  predictions <- predict(xgb_model, test)
  label <- as.numeric(predictions > 0.5)
  #Add the gene and enhancer id's
  predictions <- cbind(features_all[, c(1,2)], predictions, label)
  return(predictions)
}
