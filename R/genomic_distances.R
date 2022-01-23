
#' Compute genomic distances
#'
#' @param x A dataframe of either two columns (gene and enhancer pairs) or
#' one colum of genes
#'
#'
#' @return A dataframe with the gene and enhancer pairs and the genomic distance
#' between them
#'
#' @export
#'
#' @import GenomicRanges
#' @examples
#'
#'
distances <- function(x) {
  x <- unfactorize(x)
  gencode <- as.data.frame(gencode)
  ccres_enhancer <- as.data.frame(ccres_enhancer)
  if (ncol(x) == 2) {
    distance_dataframe <- distances_gene_enhancer(x)
    return(distance_dataframe)

  }
  else{

    returned_mat <- map_genes_to_chromosomes(x)
    distance_dataframe <- distances_gene_enhancer(returned_mat)
    return(distance_dataframe)



  }

}
