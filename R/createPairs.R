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
#' #Create gene enhancer pairs
#' genes <- as.data.frame(c("ENSG00000130203.10",
#' "ENSG00000171119.3"))
#' colnames(genes) <- c("gene_id") #It is important to name the column gene_id
#' pairs <- CENTRE::createPairs(genes)
#' @export
#' @import utils
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom RSQLite dbConnect dbGetQuery dbDisconnect
#' @importFrom regioneR extendRegions
createPairs <- function(gene) {
  startTime <- Sys.time()

  #check that the user named the column correctly
  stopifnot("The column needs to be named gene_id" =
              any(names(gene) == "gene_id"))

  ## remove the "." version id of ENSEMBL ids
  gene$gene_id1 <- gsub("\\..*", "", gene$gene_id)


  ## connect to our GENCODE v40 database to get tts of the genes
  conn <- RSQLite::dbConnect(RSQLite::SQLite(),
                             system.file("extdata",
                            "Annotation.db",
			    package = "CENTRE"))
  #get chromosome and tts of our genes

  query <- paste("SELECT  gene_id1, chr, transcription_start FROM gencode WHERE gene_id1 in (",
  paste0(sprintf("'%s'", gene$gene_id1), collapse = ", "),")",sep="" )
  gene <- RSQLite::dbGetQuery(conn, query)

  #Select all of the annotation for ccres v3
  ccresEnhancer <- RSQLite::dbGetQuery(conn, "SELECT * FROM ccres_enhancer")
  RSQLite::dbDisconnect(conn)

  genesRange <- with(gene,
                     GenomicRanges::GRanges(chr,
                                            IRanges::IRanges(start = transcription_start,
                                                             end = transcription_start)))

  #extend the gene region 500Kb to the left of TTS and to the right
  genesRange <- regioneR::extendRegions(genesRange,
                                        extend.start=500000,
                                        extend.end=500000)

  enhancerRange<-  with(ccresEnhancer,
                         GenomicRanges::GRanges(V1,
                                                IRanges::IRanges(start = new_start,
                                                                 end = new_end)))


  # find the enhancers that overlap the extended gene region
  overlaps <- GenomicRanges::findOverlaps(genesRange, enhancerRange,
                                          ignore.strand = TRUE)

  ccresOverlapping <-data.frame(gene = overlaps@from, enhancer = overlaps@to)
  ccresOverlapping$gene_id1 <- gene$gene_id1[ccresOverlapping$gene]
  ccresOverlapping$enhancer_id <- ccresEnhancer$V5[ccresOverlapping$enhancer]

  #select the gene and enhancer ids
  ccresOverlapping <- ccresOverlapping[, c("gene_id1", "enhancer_id")]

  ### add a function to exclude any pairs that are not in the same chromosome

  cat(paste0('time: ', format(Sys.time() - startTime), "\n"))
  return(ccresOverlapping)
}
