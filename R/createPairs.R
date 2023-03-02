#' Create Pairs
#'
#' Creates all of the possible gene enhancer pairs at 500kb distance
#' from the given genes
#' transcription start sites.
#'
#' @param gene One column dataframe with gene ENSEMBL id's
#' @return dataframe with two columns the ENSEMBL id's without their version and
#' the enhancer id's
#'
#'
#' @examples
#' #Start by providing genes with their ENSEMBL id
#' candidates <- read.table(system.file("extdata",
#' "exampleids.txt", package = "CENTRE"), header = T)
#' #Remember to give the columns the name "gene_id"
#' colnames(candidates) <- c("gene_id")
#' #Generate the candidate pairs
#' candidate_pairs <- createPairs(candidates)
#' @export
#' @import utils
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges

createPairs <- function(gene) {
  start_time <- Sys.time()

  colnames(gene) <- c("gene_id")

  gene$gene_id1 <- gsub("\\..*", "", gene$gene_id)



  #get chromosome and tts of our genes
  ##This should be changed to annotation hub v40 but it is missing
  gene <- merge(gene,
                  gencode[, c("chr","gene_id1","transcription_start")],
                  by.x = "gene_id1",
                  by.y = "gene_id1")

  gene$startTts <- integer(nrow(gene))
  gene$endTts <- integer(nrow(gene))

  for (i in seq_along(gene)) {
    ##extend 500 kb to the left of tts

    if (gene$transcription_start[i] <= 500000) {
      ##We check that the extension doesnt fall outside of the chromosome
      gene$startTts[i] <- 1
    } else{
      gene$startTts[i] <- gene$transcription_start[i] - 500000
    }

    ##extend 500kb to the right of tts

    chrSize <- chromosomes[chromosomes[, 1] %in% gene$chr[i], 2]

    if (gene$transcription_start[i] + 500000 >= chrSize) {
      ##We check that the extension doesnt fall outside of the chromosome
      gene$endTts[i] <- chrSize
    } else {

      gene$endTts[i] <- gene$transcription_start[i] + 500000
    }

  }

  genesRange <- with(gene,
                      GenomicRanges::GRanges(chr,
                                             IRanges::IRanges(start = startTts,
                                                              end = endTts)))
  #Implement it with annotation Hub
  enhancerRange<-  with(ccres_enhancer,
                         GenomicRanges::GRanges(V1,
                                                IRanges::IRanges(start = new_start,
                                                                 end = new_end)))

  overlaps <- GenomicRanges::findOverlaps(genesRange, enhancerRange,
                                          ignore.strand = TRUE)

  ccresOverlapping <-data.frame(gene = overlaps@from, enhancer = overlaps@to)
  ccresOverlapping$gene_id1 <- gene[ccresOverlapping$gene, 1]
  ccresOverlapping$enhancer_id <- ccres_enhancer$V5[ccresOverlapping$enhancer]
  ccresOverlapping <- ccresOverlapping[, 3:4]
  cat(paste0('time: ', format(Sys.time() - start_time), "\n"))
  return(ccresOverlapping)


}
