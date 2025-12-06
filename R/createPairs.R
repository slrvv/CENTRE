#' Create Pairs
#'
#' Creates all of the possible gene enhancer pairs at 500kb distance from the
#' given genes transcription start sites. The pairs can be also computed from
#' enhancers, in that mode we collect all gene enhancer pairs at 500kb distance
#' of the enhancer middle point.
#'
#' @param ids Vector with gene ENSEMBL id's or enhancer cCREs ids.
#' @param enhancerCentered Boolean value. If true the pairs are computed from
#' enhancers. In the default setting (false) pairs are computed from genes.
#' @return data.frame with two columns the ENSEMBL id's without their version
#' and the enhancer id's (ENCODE cCREs).
#'
#'
#' @examples
#' # Create gene enhancer pairs from genes
#' ids <- c(
#'     "ENSG00000130203.10",
#'     "ENSG00000171119.3"
#' )
#' pairs <- CENTRE::createPairs(ids)
#'
#' # Create gene enhancer pairs from enhancers
#' # ids_enh <- c("EH38E3750708", "EH38E2776554")
#' # pairs_enh <- CENTRE::createPairs(ids_enh, enhancerCentered = TRUE)
#'
#' @export
#' @import utils
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom regioneR extendRegions
#' @importFrom AnnotationHub AnnotationHub
#' @importFrom CENTREannotation fetch_data
#' @importClassesFrom CENTREannotation CENTREannotDb
createPairs <- function(ids, enhancerCentered = FALSE) {
    startTime <- Sys.time()

    if (missing(ids)) {
        stop("Need to provide a vector of SCREEN or ENSEMBL ID's.\n")
    }
    if (enhancerCentered == TRUE) {
        message("Computing pairs from the enhancer IDs...\n\n")
        ids <- data.frame(enhancer_id = ids)
        ccresOverlapping <- enhancerCenteredPairs(ids)
    } else {
        # check that the user named the column correctly
        message("Computing Pairs from gene IDs...\n\n")
        ids <- data.frame(gene_id = ids)
        ccresOverlapping <- geneCenteredPairs(ids)
    }

    message(paste0("time: ", format(Sys.time() - startTime), "\n"))
    return(ccresOverlapping)
}
