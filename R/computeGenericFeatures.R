#' Compute generic features
#'
#' Computes the generic features needed by CENTRE to make predictions.
#'
#' @param pairs Data.frame with the gene-enhancer pairs of interest
#'
#' @return
#' A data frame containing the following computed features :
#'* distance: Distance between gene and enhancer.
#'* combined_tests: Combined value of the Wilcoxon tests (CAGE, DNase
#' expression, CRUP expression and DNase DNase)
#'* crup_cor: CRUP correlation scores.
#'
#' @examples
#' # Create gene enhancer pairs
#' genes <- c(
#'     "ENSG00000130203.10",
#'     "ENSG00000280087.1"
#' )
#' pairs <- CENTRE::createPairs(genes)
#' generic_features <- CENTRE::computeGenericFeatures(pairs)
#' @export
#' @import utils
#' @importFrom AnnotationHub AnnotationHub
#' @importFrom ExperimentHub ExperimentHub
#' @importFrom CENTREannotation fetch_data
#' @importClassesFrom CENTREannotation CENTREannotDb
#' @importFrom CENTREprecomputed fetch_data_precomp
#' @importClassesFrom CENTREprecomputed CENTREprecompDb
#' @importFrom dplyr inner_join left_join join_by %>% rename select
computeGenericFeatures <- function(pairs) {
    startTime <- Sys.time()
    message("Computing CENTRE generic features\n\n")
    # Pre-eliminary checks and computations
    if (missing(pairs)) {
        stop("Need to provide a dataframe of enhancer and gene pairs.")
    }

    needed_names <- c("gene_id1", "enhancer_id")
    if (!all(needed_names %in% colnames(pairs))) {
        missing_cols <- setdiff(needed_names, colnames(pairs))
        stop(
            "The following expected columns are missing: ",
            paste(missing_cols, collapse = ", ")
        )
    }

    ## remove version identifier just in case user provided it.
    pairs$gene_id1 <- gsub("\\..*", "", pairs$gene_id1)
    ## Computing the distance features
    message("Computing distance features...\n")

    featuresDistances <- computeDistances(pairs)

    message("Removing pairs with distance over 500 Kb.\n")
    featuresDistances <- featuresDistances[featuresDistances$distance
    <= 500000, ]
    ## Getting the values for the Wilcoxon tests and the CRUP correlations
    message("Getting Wilcoxon tests and CRUP correlations...\n")

    featuresDistances$pair <- paste(featuresDistances$enhancer_id,
        featuresDistances$gene_id1,
        sep = "_"
    )
    # connect to the precomputed values database
    message("Load Wilcoxon tests from ExperimentHub.\n")
    combinedTestDf <- getPrecomputedValues(
        "combinedTestData",
        "combined_tests",
        featuresDistances
    )

    message("Load CRUP correlations from ExperimentHub.\n")
    crupCorDf <- getPrecomputedValues(
        "crup_cor",
        "cor_CRUP",
        featuresDistances
    )

    # join all the data into one dataframe
    featuresGeneric <- dplyr::left_join(
        x = featuresDistances,
        y = combinedTestDf,
        by = dplyr::join_by(pair)
    )
    featuresGeneric <- dplyr::left_join(
        x = featuresGeneric,
        y = crupCorDf,
        by = dplyr::join_by(pair)
    ) %>%
        dplyr::rename(crup_cor = cor_CRUP) %>%
        dplyr::select(
            gene_id1,
            enhancer_id,
            distance,
            crup_cor,
            combined_tests
        )

    featuresGeneric[is.na(featuresGeneric)] <- 0 ## NA values
    message(paste0("time: ", format(Sys.time() - startTime), "\n"))
    return(featuresGeneric)
}
