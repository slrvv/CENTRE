
#' Compute crup
#'
#' @param input matrix of enhancers and promoters
#' @param crupinput metaData matrix as defined in the crupR documentation
#' @param cores number of cores
#'
#' @return dataframe with the crup score computed for each region
#' @export
#' @import crupR
#'
#' @examples
compute_crup <- function(input, crupinput, cores){

  normalized <- crupR::normalize(metaData = crupinput, condition = 1, replicate = 1,
                                 genome = "hg38", sequencing = "single",
                                 chroms = c("chr1"), cores = cores) ##how do we deal with the type of sequencing?
  crup_scores <- crupR::getEnhancers(data = normalized, cores = cores)

  crup_scores <- as.data.frame(crup_scores$data_matrix)

  ccres_enhancer <- data.frame(ccres_enhancer)

  regions <- subset(ccres_enhancer, ccres_enhancer$V5 %in% input$V2, select = c(V5,new_start, new_end))

  score <- subset(crup_scores, start %in% regions$new_start, select = score )
  results <- cbind(regions, score)
  return(results)



}

