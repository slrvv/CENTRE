###############################################################################
# function: convert factor type columns into character type
###############################################################################
unfactorize <- function(df) {
  for (i in which(sapply(df, class) == "factor")) {
    df[[i]] <- as.character(df[[i]])
  }

  return(df)
}

###############################################################################
# function: get the distance from gene to enhancer
###############################################################################
distances_gene_enhancer <- function(x) {
  x <- unfactorize(x)
  gencode <- unfactorize(gencode)
  x$chr_gene <- gencode[gencode$gene_id %in% x[, 1], 1]
  x$chr_enh <- ccres_enhancer[ccres_enhancer$V5 %in% x[,2],1]
  x$gene_tts <- gencode[gencode$gene_id %in% x[, 1], 10]
  x$enh_middlepoint <- ccres_enhancer[ccres_enhancer$V5 %in% x[, 2], 11]

  for(i in 1:nrow(x)){
    if(x[i,]$chr_gene != x[i,]$chr_enh){
      warning("Gene and enhancer are not in the same chromosome,
            distance will not be computed")
      x <- x[- i,] #what happens if all are in diff chromosomes
    }
  }

  x$distance <- abs(x$enh_middlepoint - x$gene_tts)

  return(x)
}

###############################################################################
# function: extend genes by 50000bp and map the gene region
# to overlapping enhancers
###############################################################################


map_genes_to_chromosomes <- function(x){
  x <- unfactorize(x)
  #get chromosome and tts of our genes
  x$chr_gene <- gencode[gencode$gene_id %in% x[, 1], 1]
  x$gene_tts <- gencode[gencode$gene_id %in% x[, 1], 10]

  #extend start coorditantes of tts
  x$start_tts <- x$gene_tts
  x$start_tts[x$gene_tts <= 50000] <- 1
  x$start_tts[x$gene_tts > 50000] <- x$start_tts[ x$gene_tts > 50000] - 50000

  #extend end coordinates of tts
  x$end_tts <- x$gene_tts
  x$end_tts[x$gene_tts + 50000 > chromosomes[chromosomes[, 1] %in% x$chr_gene, 2]] <- chromosomes[chromosomes[, 1] %in% x$chr_gene, 2]
  x$end_tts[x$gene_tts + 50000 <= chromosomes[chromosomes[, 1] %in% x$chr_gene, 2]] <-
    x$end_tts[x$gene_tts + 50000 <= chromosomes[chromosomes[, 1] %in% x$chr_gene, 2]] + 50000


  genes_range <- with(x, GenomicRanges::GRanges(chr_gene, IRanges::IRanges(start_tts, end_tts)))

  enhancer_range<-  with(ccres_enhancer, GenomicRanges::GRanges(V1, IRanges::IRanges(start=new_start,end=new_end)))

  overlaps <- GenomicRanges::findOverlaps(genes_range, enhancer_range,
                                          ignore.strand = T)

  cres_overlaping <-data.frame(gene=overlaps@from,enhancer=overlaps@to)
  cres_overlaping$gene_id <- x[cres_overlaping$gene,1]
  cres_overlaping$enhancer_id <- ccres_enhancer$V5[cres_overlaping$enhancer]
  cres_overlaping <- cres_overlaping[,3:4]
  return(cres_overlaping)


}

###############################################################################
# function: compute the distances between genes and enhancers
###############################################################################


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
    print(returned_mat)
    distance_dataframe <- distances_gene_enhancer(returned_mat)
    return(distance_dataframe)



  }

}
