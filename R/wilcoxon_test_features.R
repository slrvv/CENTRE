
convert_na <- function(data){
  data[is.na(data)] <- 0
  return(data)

}



wilcoxon_test_features <- function(x){
  x$cage_wilcoxon_test <- cage_test_data[x$pair,3]
  x$dhsexp_wilcoxon_test  <-  dhsexp_test_data[x$pair,3]
  x$crupexp_wilcoxon_test  <-  crupexp_test_data[x$pair,3]
  x$dhsdhs_wilcoxon_test  <- dhsdhs_test_data[x$pair,3]
  return(x)
}




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

################################################################################
# function: gets the precomputed values of the CRUP expression wilcoxon test for our input
# gene enhancer pairs
################################################################################

crupexp_wilcoxon_test <- function(x){
  x$crupexp_wilcoxon_test  <-  crupexp_test_data[x$pair,3]
  return(x)

}

################################################################################
# function: gets the precomputed values of the DHS-DHS wilcoxon test for our input
# gene enhancer pairs
################################################################################

dhsdhs_wilcoxon_test <- function(x){
  x$dhsdhs_wilcoxon_test  <- dhsdhs_test_data[x$pair,3]
  return(x)

}


################################################################################
# function: get the RNA seq TPM values for our genes
################################################################################

get_rnaseq <- function(x, tpmfile){
  tpmfile$gene_id2 <- gsub("\\..*","",tpmfile[,1])
  rnaseq_tpmvalue <- subset(tpmfile, tpmfile$gene_id2 %in% x$gene_id2,
                    select = V3)
  return(rnaseq_tpmvalue)
}
