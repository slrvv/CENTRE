################################################################################
# function: gets the precomputed values of the wilcoxon test for our input gene
#           enhancer pairs
################################################################################

cage_wilcoxon_test <- function(x){
  x$cage_wilcoxon_test <- cage_test_data[x$pair,3]
  return(x)

}
################################################################################
# function: gets the precomputed values of the CRUP correlation for our input
# gene enhancer pairs
################################################################################

crup_correlations <- function(x){
  x$crup_correlations <- crup_cor[x$pair,3]
  return(x)

}

################################################################################
# function: gets the precomputed values of the DNAase wilcoxon test for our input
# gene enhancer pairs
################################################################################

dhsexp_wilcoxon_test <- function(x){
  x$dhsexp_wilcoxon_test  <-  dhsexp_test_data[x$pair,3]
  return(x)

}


