#' Create Pairs
#'
#' Creates all of the possible gene enhancer pairs at 500kb distance
#' from the given genes transcription start sites. The pairs can be also computed
#' from enhancers, in that mode we collect all gene enhancer pairs at 500kb distance
#' of the enhancer middle point.
#'
#' @param ids One column dataframe with gene ENSEMBL id's or enhancer cCREs ids
#' @param enhancerCentered Boolean value. If true the pairs are computed from 
#' enhancers. In the default setting (false) pairs are computed from genes.
#' @return dataframe with two columns the ENSEMBL id's without their version and
#' the enhancer id's (ENCODE cCREs)
#'
#'
#' @examples
#' #Create gene enhancer pairs
#' ids <- as.data.frame(c("ENSG00000130203.10",
#' "ENSG00000171119.3"))
#' colnames(ids) <- c("gene_id") #It is important to name the column gene_id or 
#' # enhancer_id
#' pairs <- CENTRE::createPairs(ids)
#' @export
#' @import utils
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom RSQLite dbConnect dbGetQuery dbDisconnect
#' @importFrom regioneR extendRegions
createPairs <- function(ids, enhancerCentered = FALSE) {
  startTime <- Sys.time()
  
  if (enhancerCentered == TRUE){
    stopifnot("The column needs to be named enhancer_id" =
                any(names(ids) == "enhancer_id"))
    cat("Pairs are computed from the enhancer IDs\n")
    ccresOverlapping <- enhancerCenteredPairs(ids)
   
  } else {
    #check that the user named the column correctly
    stopifnot("The column needs to be named gene_id" =
                any(names(ids) == "gene_id"))
    cat("Pairs are computed from the gene IDs\n")
    ccresOverlapping <- geneCenteredPairs(ids)
    
  }
  
  cat(paste0("time: ", format(Sys.time() - startTime), "\n"))
  return(ccresOverlapping)
}
