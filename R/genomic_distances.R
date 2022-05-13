##################################################################
#definition: output messages
##################################################################

done <- function() {
  cat(".. done.\n")
}
skip <- function() {
  cat("\t ..")
}
startPart <- function(m) {
  cat(paste0("\n--- ", m, " ---\n\n"))
}
endPart <- function() {
  cat("\n\t>>> All done!\n")
}


###############################################################################
# lookup table: end position of the chromosomes in hg38 of human genome
###############################################################################

chromosomes <- data.frame(
  chr= c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
          "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
          "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22",
          "chrX", "chrY"),
  position = c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
               159345973, 145138636, 138394717,133797422, 135086622, 133275309,
               114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
               58617616, 64444167, 46709983, 50818468, 156040895, 57227415)
)


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
  x$end_tts[x$gene_tts + 50000  >= chromosomes[chromosomes[, 1] %in% x$chr_gene, 2]] <- chromosomes[chromosomes[, 1] %in% x$chr_gene, 2]
  x$end_tts[x$gene_tts + 50000 < chromosomes[chromosomes[, 1] %in% x$chr_gene, 2]] <- x$end_tts[x$gene_tts + 50000 < chromosomes[chromosomes[, 1] %in% x$chr_gene, 2]] + 50000
  genes_range <- with(x, GenomicRanges::GRanges(chr_gene, IRanges::IRanges(start = start_tts, end = end_tts)))

  enhancer_range<-  with(ccres_enhancer, GenomicRanges::GRanges(V1, IRanges::IRanges(start=new_start,end=new_end)))

  overlaps <- GenomicRanges::findOverlaps(genes_range, enhancer_range,
                                          ignore.strand = T)

  cres_overlaping <-data.frame(gene=overlaps@from,enhancer=overlaps@to)
  print(cres_overlaping)
  cres_overlaping$gene_id <- x[cres_overlaping$gene,1]
  cres_overlaping$enhancer_id <- ccres_enhancer$V5[cres_overlaping$enhancer]
  cres_overlaping <- cres_overlaping[,3:4]
  return(cres_overlaping)


}



###############################################################################
# function: check if file exists
###############################################################################

check_file <- function(f) {
  if (!(file.exists(f))) {
    message <- paste0("File ", f, " does not exist.\n")
    stop(message)
  }
}
