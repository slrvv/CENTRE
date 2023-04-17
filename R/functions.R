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
  chr = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
         "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
         "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22",
         "chrX", "chrY"),
  position = c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
               159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
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
# function: extend 500Kbp
###############################################################################
extend <- function(gene){
  ##extend 500 kb to the left of tts
  if (gene$transcription_start[i] <= 500000) {
    ##We check that the extension doesnt fall outside of the chromosome
    start_tts[i] <- 1
  } else{
    start_tts[i] <- gene$transcription_start[i] - 500000
  }

  ##extend 500kb to the right of tts

  chr_size <- chromosomes[chromosomes[, 1] %in% gene$chr[i], 2]

  if (gene$transcription_start[i] + 500000 >= chr_size) {
    ##We check that the extension doesnt fall outside of the chromosome
    end_tts[i] <- chr_size
  } else {

    end_tts[i] <- gene$transcription_start[i] + 500000
  }


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

###############################################################################
# function: get the distance from gene to enhancer
###############################################################################
computeDistances <- function(x) {
  # connect to annotation dataBase
  conn <- RSQLite::dbConnect(RSQLite::SQLite(),
                              system.file("extdata",
                             "Annotation.db",
                             package = "CENTRE"))
  #get chromosome and tts of our genes
  query <- paste("SELECT  gene_id1, chr, transcription_start FROM gencode WHERE gene_id1 in (",
  paste0(sprintf("'%s'", x$gene_id2), collapse = ", "),")",sep="" )
  gencode <- RSQLite::dbGetQuery(conn, query)
  
  #get chr and middle point of enhancers
  query_enh <-  paste("SELECT  V5, V1, middle_point FROM ccres_enhancer WHERE V5 in (",
  paste0(sprintf("'%s'", x$enhancer_id), collapse = ", "),")",sep="" )
  ccres_enhancer <- RSQLite::dbGetQuery(conn, query_enh)
  RSQLite::dbDisconnect(conn)

  #Get the chr gene_id and transcription_start from gencode annotation
  #Getting the chrosomes and the middle points for the provided enhancers
  result <- merge(x,
                  ccres_enhancer[, c("V1", "V5", "middle_point")],
                   by.x = "enhancer_id",
                   by.y = "V5") #change V1 and V5 to more meaningful names
  #Getting the chrosomes and transcription start sites for the provided genes
  result <- merge(result,
                  gencode[, c("chr", "gene_id1", "transcription_start")],
                  by.x = "gene_id2",
                  by.y = "gene_id1")

  cat("Removing all gene enhancer pairs that are not in the same chromosome.\n")

  result <- result[(result$V1 == result$chr), ]

  result$distance <- result$middle_point - result$transcription_start

  return(result)
}





###############################################################################
# function: get scores for enhancers
###############################################################################
compute_crup_enhancer <- function(regions_enhancer,
                                  list_enh,
                                  crup_scores,
                                  promprob = F) {

  #Overlapping the  enhancer ranges with the crup scores
  enhancer_ranges <- with(regions_enhancer,
                          GenomicRanges::GRanges(V1,
                                                 IRanges::IRanges(start = new_start,
                                                                  end = new_end)))

  hits_crup <- GenomicRanges::findOverlaps(enhancer_ranges, crup_scores)
  cres_EP <- data.frame(cres = hits_crup@from, EP = hits_crup@to)

  cres_EP$cres_name <- list_enh$enhancer_id[cres_EP$cres]
  if (promprob == T) {
    cres_EP$PP_enhancer <- GenomicRanges::elementMetadata(crup_scores)$probP[cres_EP$EP]
    cres_EP$cres_name <- factor(cres_EP$cres_name)
    cres_EP$bin <- rep(1:5, times = nrow(regions_enhancer))

    trial <- stats::reshape(cres_EP[,3:5],
                   idvar = "cres_name",
                   timevar = "bin",
                   direction = "wide",
                   v.names = "PP_enhancer")
  } else {
    cres_EP$EP_enhancer <- GenomicRanges::elementMetadata(crup_scores)$prob[cres_EP$EP]
    cres_EP$cres_name <- factor(cres_EP$cres_name)
    cres_EP$bin <- rep(1:5, times = nrow(regions_enhancer))


    trial <-stats::reshape(cres_EP[,3:5],
                   idvar = "cres_name",
                   timevar = "bin",
                   direction = "wide",
                   v.names = "EP_enhancer")
  }



  return(trial)




}

###############################################################################
# function: get  scores for promoters
###############################################################################
compute_crup_promoter <- function(regions_prom,
                                  list_prom,
                                  crup_scores,
                                  promprob = F) {
  #Overlapping with CRUP scores
  genes_ranges <- with(regions_prom,
                       GenomicRanges::GRanges(chr,
                                              IRanges::IRanges(start = new_start,
                                                               end = new_end)))

  hits_crup <- GenomicRanges::findOverlaps(genes_ranges, crup_scores)
  cres_EP <- data.frame(promoter = hits_crup@from, EP = hits_crup@to)
  cres_EP$gene_name <- list_prom$gene_id[cres_EP$promoter]
  if (promprob == T) {
    cres_EP$PP_promoter <- GenomicRanges::elementMetadata(crup_scores)$probP[cres_EP$EP]
    cres_EP$gene_name <- factor(cres_EP$gene_name)
    #Returning the probabilities in bins
    cres_EP$bin <- rep(1:5, nrow(regions_prom))
    trial <- stats::reshape(cres_EP[,3:5],
                   idvar = "gene_name",
                   timevar = "bin",
                   direction = "wide",
                   v.names = "PP_promoter")
  } else {
    cres_EP$EP_promoter <- GenomicRanges::elementMetadata(crup_scores)$prob[cres_EP$EP]
    cres_EP$gene_name <- factor(cres_EP$gene_name)
    #Returning the probabilities in bins
    cres_EP$bin <- rep(1:5, nrow(regions_prom))
    trial <- stats::reshape(cres_EP[,3:5],
                   idvar = "gene_name",
                   timevar = "bin",
                   direction = "wide",
                   v.names = "EP_promoter")

  }

  return(trial)
}

###############################################################################
# function: get scores for distance between enhancer and promoter
###############################################################################

compute_crup_reg_distance_enh <- function(input, prediction) {
  ##Check if the distances are negative and flip the start and end around
  input$bstart <- input$middle_point
  input$bstart[input$distance > 0] <- input$transcription_start[input$distance > 0]
  input$bend <- input$transcription_start
  input$bend[input$distance > 0] <- input$middle_point[input$distance > 0]

  #Make the gene enhancer pairs into ranges
  input$pair <- paste0(input$gene_id, "_", input$enhancer_id)
  input <- input[!(duplicated(input$pair)), ]
  between_ranges <- with(input,
                         GenomicRanges::GRanges(chr,
                                                IRanges::IRanges(start = bstart,
                                                                 end = bend)))

  hits_enh <- GenomicRanges::findOverlaps(between_ranges, prediction)
  cres_EP <- data.frame(between = hits_enh@from,
                        EP_reg_distance = GenomicRanges::elementMetadata(prediction)$prob[hits_enh@to])

  bins <- as.data.frame(table(cres_EP$between))

  cres_EP1 <- cres_EP[cres_EP$EP_reg_distance > 0.5, ]
 
  if (nrow(cres_EP1) == 0) {

    bins_pos <- as.data.frame(matrix(c(seq(1:nrow(input)),
                                         rep(0, times = nrow(input))),
                                       nrow = nrow(input),
                                       ncol =  2))
    colnames(bins_pos) <- c("Var1", "Freq")
  } else {
    bins_pos <- as.data.frame(table(cres_EP1$between))
  }

  all_bins <- merge(bins, bins_pos, by.x = "Var1", by.y = "Var1", all.x = TRUE)
  all_bins[is.na(all_bins)] <- 0
  colnames(all_bins) <- c("pair", "bins", "bins_pos")

  input$reg_dist_enh <- all_bins$bins_pos
  input$norm_reg_dist_enh <- all_bins$bins_pos / all_bins$bins



  return(input)



}

###############################################################################
# function: get scores for distance between enhancer and promoter
###############################################################################

compute_crup_reg_distance_prom <- function(input, prediction) {
  ##Check if the distances are negative and flip the start and end around
  input$bstart <- input$middle_point
  input$bstart[input$distance > 0] <- input$transcription_start[input$distance > 0]
  input$bend <- input$transcription_start
  input$bend[input$distance > 0] <- input$middle_point[input$distance > 0]

  #Make the gene enhancer pairs into ranges
  input$pair <- paste0(input$gene_id, "_", input$enhancer_id)
  input <- input[!(duplicated(input$pair)), ]
  between_ranges <- with(input,
                         GenomicRanges::GRanges(chr,
                                                IRanges::IRanges(start = bstart,
                                                                 end = bend)))

  hits_enh <- GenomicRanges::findOverlaps(between_ranges, prediction)
  cres_EP <- data.frame(between = hits_enh@from,
                        EP_reg_distance = GenomicRanges::elementMetadata(prediction)$probP[hits_enh@to])

  bins <- as.data.frame(table(cres_EP$between))

  cres_EP1 <- cres_EP[cres_EP$EP_reg_distance > 0.5, ]

  if (nrow(cres_EP1) == 0) {

    bins_pos <- as.data.frame(matrix(c(seq(1:nrow(input)),
                                         rep(0, times = nrow(input))),
                                       nrow = nrow(input),
                                       ncol =  2))
    colnames(bins_pos) <- c("Var1", "Freq")
  } else {
    bins_pos <- as.data.frame(table(cres_EP1$between))
  }


  all_bins <- merge(bins, bins_pos, by.x = "Var1", by.y = "Var1", all.x = TRUE)
  all_bins[is.na(all_bins)] <- 0
  colnames(all_bins) <- c("pair", "bins", "bins_pos")

  input$reg_dist_prom <- all_bins$bins_pos
  input$norm_reg_dist_prom <- all_bins$bins_pos / all_bins$bins


  return(input)



}

################################################################################
# function: Makes the query for Precomputed.db
################################################################################

queryMaker <- function(table, feature, x) {
  query <- paste("SELECT ", feature, ", pair FROM ", table, " WHERE pair in (",
                paste0(sprintf("'%s'", x$pair), collapse = ", "), ")", sep = "")
  return(query)
}

################################################################################
# function: gets the precomputed values of the Wilcoxon tests from
# PrecomputedData.db
################################################################################

getPrecomputedValues <- function(table, feature, x, conn) {

  query <- queryMaker(table, feature, x)


  df_return <- RSQLite::dbGetQuery(conn, query)
  rownames(df_return) <- df_return$pair
  df_return$pair <- NULL
  return(df_return)
}




################################################################################
# function: get the RNA seq TPM values for our genes
################################################################################

get_rnaseq <- function(x, tpmfile) {
  tpmfile$gene_id2 <- gsub("\\..*", "", tpmfile[, 1])

  x <- merge(x, tpmfile[,c(3,4)],by.x = "gene_id2", by.y = "gene_id2" )



  return(x)
}
