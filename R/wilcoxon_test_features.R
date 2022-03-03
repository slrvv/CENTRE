################################################################################
# function: gets the precomputed values of the wilcoxon test for our input gene
#           enhancer pairs
################################################################################

cage_wilcoxon_test <- function(x){
  x$gene_id2 <- gsub("\\..*","",x[,1])
  x$pair <- paste( x[,2],x$gene_id2, sep= "_")
  x$cage_wilcoxon_test <- all_genes[x$pair,3]
  return(x)

}
################################################################################
# function: gets the precomputed values of the CRUP correlation for our input
# gene enhancer pairs
################################################################################

crup_correlations <- function(x){
  x$gene_id2 <- gsub("\\..*","",x[,1])
  x$pair <- paste( x[,2],x$gene_id2, sep= "_")
  x$crup_correlations <- all_genes[x$pair,3]
  return(x)

}

################################################################################
# function: gets the precomputed values of the DNAase wilcoxon test for our input
# gene enhancer pairs
################################################################################

dhsexp_wilcoxon_test <- function(x){
  x$gene_id2 <- gsub("\\..*","",x[,1])
  x$pair <- paste( x[,2],x$gene_id2, sep= "_")
  x$dhsexp_wilcoxon_test  <- all_genes[x$pair,3]
  return(x)

}
