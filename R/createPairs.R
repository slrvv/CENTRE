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
<<<<<<< Updated upstream
createPairs <- function(gene) {
  start_time <- Sys.time()

  colnames(gene) <- c("gene_id")

  gene$gene_id1 <- gsub("\\..*", "", gene$gene_id)


=======
#' @importFrom regioneR extendRegions
createPairs <- function(gene) {
  startTime <- Sys.time()

  #check that the user named the column correctly
  stopifnot("The column needs to be named gene_id" =
              any(names(gene) == "gene_id"))

  ## remove the "." version id of ENSEMBL ids
  gene$gene_id1 <- gsub("\\..*", "", gene$gene_id)


  ## connect to our GENCODE v40 database to get tts of the genes
>>>>>>> Stashed changes
  conn <- RSQLite::dbConnect(RSQLite::SQLite(),
                             system.file("extdata",
                            "Annotation.db",
			    package = "CENTRE"))
<<<<<<< Updated upstream
  #get chromosome and tts of our genes
  query <- paste("SELECT  gene_id1, chr, transcription_start FROM gencode WHERE gene_id1 in (",
  paste0(sprintf("'%s'", gene$gene_id1), collapse = ", "),")",sep="" )
  gene <- RSQLite::dbGetQuery(conn, query)
  #Select all of the annotation for ccres
  ccres_enhancer <- RSQLite::dbGetQuery(conn, "SELECT * FROM ccres_enhancer")
  RSQLite::dbDisconnect(conn)
  gene$startTts <- integer(nrow(gene))
  gene$endTts <- integer(nrow(gene))
  for (i in seq_len(nrow(gene)))
 {
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
=======

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
>>>>>>> Stashed changes
                         GenomicRanges::GRanges(V1,
                                                IRanges::IRanges(start = new_start,
                                                                 end = new_end)))

<<<<<<< Updated upstream
=======
  # find the enhancers that overlap the extended gene region
>>>>>>> Stashed changes
  overlaps <- GenomicRanges::findOverlaps(genesRange, enhancerRange,
                                          ignore.strand = TRUE)

  ccresOverlapping <-data.frame(gene = overlaps@from, enhancer = overlaps@to)
<<<<<<< Updated upstream
  ccresOverlapping$gene_id1 <- gene[ccresOverlapping$gene, 1]
  ccresOverlapping$enhancer_id <- ccres_enhancer$V5[ccresOverlapping$enhancer]
  ccresOverlapping <- ccresOverlapping[, 3:4]
  cat(paste0('time: ', format(Sys.time() - start_time), "\n"))
  return(ccresOverlapping)


=======

  ccresOverlapping$gene_id1 <- gene$gene_id1[ccresOverlapping$gene]
  ccresOverlapping$enhancer_id <- ccresEnhancer$V5[ccresOverlapping$enhancer]

  #select the gene and enhancer ids
  ccresOverlapping <- ccresOverlapping[, c("gene_id1", "enhancer_id")]

  ### add a function to exclude any pairs that are not in the same chromosome

  cat(paste0('time: ', format(Sys.time() - startTime), "\n"))
  return(ccresOverlapping)

>>>>>>> Stashed changes
}
