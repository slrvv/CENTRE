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
#' pairs <- data.frame(
#'     gene_id1 = c("ENSG00000105281"),
#'     enhancer_id = c("EH38E1958626")
#' )
#'
#' generic_features <- CENTRE::computeGenericFeatures(pairs)
#'
#' # Compute Cell-type features
#' eh <- ExperimentHub::ExperimentHub()
#'
#' files <- c(
#'     system.file("extdata",
#'         "example/HeLa_H3K4me1.REF_chr19_reduced.bam",
#'         package = "CENTRE"
#'     ),
#'     system.file("extdata",
#'         "example/HeLa_H3K4me3.REF_chr19_reduced.bam",
#'         package = "CENTRE"
#'     ),
#'     system.file("extdata",
#'         "example/HeLa_H3K4me3.REF_chr19_reduced.bam",
#'         package = "CENTRE"
#'     )
#' )
#' inputs <- system.file("extdata",
#'     "example/HeLa_input.REF_chr19_reduced.bam",
#'     package = "CENTRE"
#' )
#' metaData <- data.frame(
#'     HM = c("H3K4me1", "H3K4me3", "H3K27ac"),
#'     condition = c(1, 1, 1), replicate = c(1, 1, 1),
#'     bamFile = files, inputFile = rep(inputs, 3)
#' )
#'
#' tpmpath <- unname(eh[["EH9545"]])
#' tpmfile <- read.table(tpmpath,
#'     sep = "",
#'     stringsAsFactors = FALSE, header = TRUE
#' )
#' tpmfile <- tpmfile[grep("E", tpmfile$gene_id), ]
#' celltype_features <- CENTRE::computeCellTypeFeatures(metaData,
#'     replicate = 1,
#'     input.free = FALSE,
#'     cores = 1,
#'     sequencing = "single",
#'     tpmData = tpmfile,
#'     pairs = pairs
#' )
#' # Finally compute the predictions
#' predictions <- centrePrediction(celltype_features, generic_features)
#'
#' @export
#' @importFrom stats predict
#' @import utils
#' @importFrom xgboost xgb.load xgb.DMatrix
#' @importFrom dplyr select inner_join %>%
centrePrediction <- function(features_celltype,
    features_generic,
    model = NULL) {
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
        xgb_model <- readRDS(system.file("extdata",
            "centre2_final_model.rds",
            package = "CENTRE"
        ))
    } else {
        xgb_model <- xgboost::xgb.load(model)
    }

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
