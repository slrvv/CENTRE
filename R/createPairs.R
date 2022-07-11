#' Create Pairs
#'
#' Creates all of the possible gene enhancer pairs at 500kb distance from the given genes
#' transcription start sites.
#'
#' @param gene One column dataframe with gene ENSEMBL id's
#' @return dataframe with two columns the ENSEMBL id's without their version and
#' the enhancer id's
#' @export
#'
#' @examples
createPairs <- function(gene){

  colnames(gene) <- c("gene_id")

  gene$gene_id1 <- gsub("\\..*","",gene$gene_id)



  #get chromosome and tts of our genes
  gene<- merge(gene,
                  gencode[,c('chr','gene_id1','transcription_start')],
                  by.x='gene_id1',
                  by.y= 'gene_id1')

  start_tts <- integer(nrow(gene))
  end_tts <- integer(nrow(gene))
  for (i in 1:nrow(gene)){
    ##extend 500 kb to the left of tts

    if (gene$transcription_start[i] <= 500000){
      ##We check that the extension doesnt fall outside of the chromosome
      start_tts[i] <- 1
    } else{
      start_tts[i] <- gene$transcription_start[i] - 500000
    }

    ##extend 500kb to the right of tts

    chr_size <- chromosomes[chromosomes[,1] %in% gene$chr[i], 2]

    if (gene$transcription_start[i] + 500000 >= chr_size) {
      ##We check that the extension doesnt fall outside of the chromosome
      end_tts[i] <- chr_size
    } else {

      end_tts[i] <- gene$transcription_start[i] + 500000
    }

  }

  gene<- cbind(gene, start_tts, end_tts)

  genes_range <- with(gene,
                      GenomicRanges::GRanges(chr,
                                             IRanges::IRanges(start = start_tts,
                                                              end = end_tts)))

  enhancer_range<-  with(ccres_enhancer,
                         GenomicRanges::GRanges(V1,
                                                IRanges::IRanges(start=new_start,
                                                                 end=new_end)))

  overlaps <- GenomicRanges::findOverlaps(genes_range, enhancer_range,
                                          ignore.strand = T)

  cres_overlaping <-data.frame(gene=overlaps@from,enhancer=overlaps@to)
  cres_overlaping$gene_id1 <- gene[cres_overlaping$gene,1]
  cres_overlaping$enhancer_id <- ccres_enhancer$V5[cres_overlaping$enhancer]
  cres_overlaping <- cres_overlaping[,3:4]
  return(cres_overlaping)


}




