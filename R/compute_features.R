compute_features <- function(file, metaData, cores){

  ##Check if input and number of arguments is correct

  ##Checking if the file exists

  check_file(file)
  x <- read.table(file, header = T, sep = "\t", stringsAsFactors = F)

  ##Add the nice printing of stuff

  print("Computing features...")

  print("\n")
  print("Computing distance feature")
  features_distances <- distances(x)
  print("Computing CRUP score features")
  chr <- unique(features_distances$chr_gene)
  features_crup <- compute_features_crup_scores(x, metaData, cores, chr)
  print("Computing Wilcoxon tests and CRUP correlations")
  x$gene_id2 <- gsub("\\..*","",x[,1])
  x$pair <- paste( x[,2],x$gene_id2, sep= "_")
  wilcoxon_cage <- cage_wilcoxon_test(x)
  crup_cor <- crup_correlations(x)
  wilcoxon_dnase <- dhsexp_wilcoxon_test(x)


  features_table <- cbind(
    features_distances[c("V1", "V2", "chr_gene", "distance")],
    features_crup[,3:ncol(features_crup)],
    wilcoxon_cage$cage_wilcoxon_test,
    crup_cor$crup_correlations,
    wilcoxon_dnase$dhsexp_wilcoxon_test)

  features_table[is.na(features_table)] <- 0

  return(features_table)
}
