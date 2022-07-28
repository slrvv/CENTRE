centrePrediction <- function(features_generic, features_celltype, model){
  ##Merge the generic features and the cell type features
  features_generic$distance <- abs(features_generic$distance)
  features_all <- cbind(features_celltype, features_generic[,c(6,7,8)]) ##change it to a merge later

  ##Loading the xgboost model
  xgb_model <- xgboost::xgb.load(model)

  ##Transforming data
  feature_matrix <- as.matrix(features_all[, c(-1, -2)])
  test <- xgboost::xgb.DMatrix(data = feature_matrix)

  ##Predicting
  predictions <- predict(xgb_model, test)

  #Add the gene and enhancer id's
  predictions <- cbind(features_all[, c(1,2)], predictions)
  return(predictions)

}
