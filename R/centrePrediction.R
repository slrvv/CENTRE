#' CENTRE predicton
#'
#' Predicts if the enhancer gene pairs are interacting from the generic
#' features and cell type features.
#'
#' @param features_celltype The cell type specific features returned by
#' `CENTRE::computeCellTypeFeatures()`
#' @param features_generic The generic features returned by the function
#' `CENTRE::computeGenericFeatures()`
#' @param model Path to the model the predictions will be computed on,
#' The default is the CENTRE model.
#'
#' @return Data.frame containing the enhancer gene pairs and the probability of
#' them interacting based on CENTRE model
#'
#' @examples
#'\dontrun{
#' pairs <- readRDS(file = system.file("extdata",
#'        "input_cellType_pairs.rds",
#'        package = "CENTRE"
#'    ))
#'
#' generic_features <- computeGenericFeatures(pairs)
#'
#' celltype_features <- readRDS(file = system.file("extdata",
#'        "expected_cellType_HeLa_reduced.rds",
#'        package = "CENTRE"
#'    ))
#' # Finally compute the predictions
#' predictions <- centrePrediction(celltype_features, generic_features)
#'}
#' @export
#' @importFrom stats predict
#' @import utils
#' @importFrom xgboost xgb.load xgb.DMatrix
#' @importFrom dplyr select inner_join %>% .data
#' @importFrom R.utils gunzip
centrePrediction <- function(features_celltype,
    features_generic,
    model = NULL) {
    pair <- enhancer_id <- gene_id1 <- NULL
    distance <- crup_cor <- combined_tests <- NULL
    # Merge the generic features and the cell type features

    start_time <- Sys.time()
    message("Computing CENTRE predictions...")
    features_generic$distance <- abs(features_generic$distance)
    # generate the pair id to merge both feature sets
    features_generic$pair <- paste(features_generic$enhancer_id,
        features_generic$gene_id1,
        sep = "_"
    )
    features_celltype$pair <- paste(features_celltype$enhancer_id,
        features_celltype$gene_id1,
        sep = "_"
    )

    # remove all non-feature columns except pair id
    features_generic <- features_generic %>% dplyr::select(
        distance,
        crup_cor,
        combined_tests,
        pair
    )
    features_celltype <- features_celltype %>% dplyr::select(-c(
        gene_id1,
        enhancer_id
    ))
    # mergeboth datasets
    features_all <- dplyr::inner_join(features_celltype,
        features_generic,
        by = dplyr::join_by(pair)
    )

    ## Loading the xgboost model
    if (is.null(model)) {
        gz_file <- system.file("extdata",
            "centre2_final_model.json.gz",
            package = "CENTRE"
            )
        tmp_model_file <- tempfile(fileext = ".json")
        R.utils::gunzip(gz_file, destname=tmp_model_file, remove = FALSE)
        model <- tmp_model_file
    }
    xgb_model <- xgboost::xgb.load(model)


    ## Transforming data
    pairs <- features_all$pair
    features_all <- features_all %>% dplyr::select(-c(pair))
    colnames(features_all) <- NULL

    feature_matrix <- xgboost::xgb.DMatrix(data.matrix(features_all))
    ## Predicting
    score <- predict(xgb_model, feature_matrix)
    label <- as.numeric(score > 0.5)
    # Add the gene and enhancer id's
    predictions <- data.frame(
        pairs = pairs,
        score = score,
        label = label
    )
    message(paste0("time: ", format(Sys.time() - start_time), "\n"))
    return(predictions)
}
