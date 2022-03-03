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
  distances <- c()
  for (i in 1:nrow(x)) {
    if(ccres[ccres$V5 == x[i, 2], 1] != gencode[gencode$gene_id == x[i, 1], 1] ){
      warning("Gene and enhancer are not in the same chromosome, distance will not be computed")
    }
    else {
      eldist <- abs(ccres[ccres$V5 == x[i, 2], 11]
                    - gencode[gencode$gene_id == x[i, 1], 10])

      distances <- c(distances, eldist)
    }

  }
  final <- as.data.frame(cbind(x[, 1], x[, 2], distances))
  final <- unfactorize(final)
  return(final)
}

###############################################################################
# function: extend genes by 50000bp and map the gene region
# to overlapping enhancers
###############################################################################

### this can be done in a much nicer way
map_genes_to_chromosomes <- function(x){
  start <- c()
  end <- c()
  chrs <- c()
  for (i in 1:nrow(x)) {
    tts <- gencode[gencode$gene_id == x[i, 1], 10]
    chr <- gencode[gencode$gene_id == x[i, 1], 1]
    chrs <- c(chrs, chr)
    if (tts <= 50000) {
      start <- c(start, 1)
    } else {
      start <- c(start, tts - 50000)
    }

    if (tts + 50000 > chromosomes[chromosomes[, 1] == chr, 2]) {
      end <- c(end, chromosomes[chromosomes[, 1] == chr, 2])
    } else {
      end <- c(end, tts + 50000)
    }
  }

  genesdata <- as.data.frame(cbind(chr, start, end))
  colnames(genesdata) <- c("seqnames", "start", "end")
  genes_range <- GenomicRanges::makeGRangesFromDataFrame(genesdata)
  GenomicRanges::elementMetadata(genes_range) <- x[, 1]


  enhancerdata <- ccres_enhancer[, c("V1", "V2", "V3")]
  colnames(enhancerdata) <- c("seqnames", "start", "end")
  enhancer_range <- GenomicRanges::makeGRangesFromDataFrame(enhancerdata)


  GenomicRanges::elementMetadata(enhancer_range) <- ccres_enhancer[, "V5"]
  overlaps <- GenomicRanges::findOverlaps(genes_range, enhancer_range,
                                          ignore.strand = T)

  returned_mat <- as.data.frame(cbind(x[overlaps@from, ],
                                      ccres_enhancer[overlaps@to, "V5"]))
  return(returned_mat)


}
