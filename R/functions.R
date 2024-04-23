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
  for (i in which(vapply(df, class) == "factor")) {
    df[[i]] <- as.character(df[[i]])
  }

  return(df)
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
  paste0(sprintf("'%s'", x$gene_id2), collapse = ", "), ")", sep = "")
  gencode <- RSQLite::dbGetQuery(conn, query)
  #get chr and middle point of enhancers
  query_enh <-  paste("SELECT  V5, V1, middle_point FROM ccres_enhancer WHERE V5 in (",
  paste0(sprintf("'%s'", x$enhancer_id), collapse = ", "), ")", sep = "")
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
  result$distance <- abs(result$middle_point - result$transcription_start)
  return(result)
}


###############################################################################
# function: get scores for enhancers
###############################################################################
compute_crup_enhancer <- function(regions_enhancer,
                                  crup_scores,
                                  promprob = FALSE) {

  #Overlapping the  enhancer ranges with the crup scores
  enhancer_ranges <- with(regions_enhancer,
                          GenomicRanges::GRanges(V1,
                                                 IRanges::IRanges(start = new_start.x,
                                                                  end = new_end.x),
                                                 enhancer_id = enhancer_id))

  enhancer_ranges <- unique(enhancer_ranges)
  hits_crup <- GenomicRanges::findOverlaps(enhancer_ranges, crup_scores)
  cres_EP <- data.frame(cres = hits_crup@from, EP = hits_crup@to)
  cres_EP$enhancer_id <- GenomicRanges::elementMetadata(enhancer_ranges)$enhancer_id[cres_EP$cres]

  if (promprob == TRUE) {
    # get the PP crup scores
    cres_EP$PP_prob_enh <- GenomicRanges::elementMetadata(crup_scores)$probP[cres_EP$EP]
    cres_EP$enhancer_id <- factor(cres_EP$enhancer_id)
    cres_EP$bin <- rep(1:5, times = length(enhancer_ranges))

    trial <- stats::reshape(cres_EP[, 3:5],
                   idvar = "enhancer_id",
                   timevar = "bin",
                   direction = "wide",
                   v.names = "PP_prob_enh")
  } else {
    #get the EP crup scores
    cres_EP$EP_prob_enh <- GenomicRanges::elementMetadata(crup_scores)$prob[cres_EP$EP]
    cres_EP$enhancer_id <- factor(cres_EP$enhancer_id)
    cres_EP$bin <- rep(1:5, times = length(enhancer_ranges))

    trial <- stats::reshape(cres_EP[, 3:5],
                   idvar = "enhancer_id",
                   timevar = "bin",
                   direction = "wide",
                   v.names = "EP_prob_enh")
  }

  return(trial)
}

###############################################################################
# function: get  scores for promoters
###############################################################################
compute_crup_promoter <- function(regions_prom,
                                  crup_scores,
                                  promprob = FALSE) {
  #Overlapping with CRUP scores
  gene_ranges <- with(regions_prom,
                       GenomicRanges::GRanges(chr,
                                              IRanges::IRanges(start = new_start.y,
                                                               end = new_end.y),
                                              gene_id2 = gene_id2))

  gene_ranges <- unique(gene_ranges)
  hits_crup <- GenomicRanges::findOverlaps(gene_ranges, crup_scores)
  cres_EP <- data.frame(promoter = hits_crup@from, EP = hits_crup@to)
  cres_EP$gene_id2 <- GenomicRanges::elementMetadata(gene_ranges)$gene_id2[cres_EP$promoter]

  if (promprob == TRUE) {
    cres_EP$PP_prob_gene <- GenomicRanges::elementMetadata(crup_scores)$probP[cres_EP$EP]
    cres_EP$gene_id2 <- factor(cres_EP$gene_id2)
    cres_EP$bin <- rep(1:5, length(gene_ranges))

    trial <- stats::reshape(cres_EP[, 3:5],
                   idvar = "gene_id2",
                   timevar = "bin",
                   direction = "wide",
                   v.names = "PP_prob_gene")
  } else {
    cres_EP$EP_prob_gene <- GenomicRanges::elementMetadata(crup_scores)$prob[cres_EP$EP]
    cres_EP$gene_id2 <- factor(cres_EP$gene_id2)
    cres_EP$bin <- rep(1:5, length(gene_ranges))

    trial <- stats::reshape(cres_EP[, 3:5],
                   idvar = "gene_id2",
                   timevar = "bin",
                   direction = "wide",
                   v.names = "EP_prob_gene")
  }

  return(trial)
}

###############################################################################
# function: make between_ranges GRanges object for reg distance calculations
###############################################################################
createBetweenRanges <- function(regions) {
  ##Check if the distances are negative and flip the start and end around
  ##compute distance as middle point - tss so if distance is negative it means
  ## tss > middle point.
  ## in the positive case bstart should be the tts and in the negative it should be
  ## the middle point
  regions$bstart <- regions$middle_point ##assign every value to the middle point
  regions$bstart[regions$distance > 0] <- regions$transcription_start[regions$distance > 0]
  ## assign bstart to the tss in cases where distance is positive
  regions$bend <- regions$transcription_start ## the same for bend
  regions$bend[regions$distance > 0] <- regions$middle_point[regions$distance > 0]
  regions$chr <- regions$V1
  #Make the gene enhancer pairs into ranges
  regions$pair <- paste(regions$enhancer_id, regions$gene_id2, sep = "_")

  between_ranges <- with(regions,
                         GenomicRanges::GRanges(chr,
                                                IRanges::IRanges(start = bstart,
                                                                 end = bend),
                                                pair = pair))
  between_ranges <- unique(between_ranges)
  return(between_ranges)
}

###############################################################################
# function: get scores for region between enhancer and promoter
###############################################################################

compute_crup_reg_distance_enh <- function(input, prediction, between_ranges) {
  ##overlap the ranges objects with predictions
  hits_enh <- GenomicRanges::findOverlaps(between_ranges, prediction)
  cres_EP <- data.frame(between = GenomicRanges::elementMetadata(between_ranges)$pair[hits_enh@from],
                        EP_reg_distance = GenomicRanges::elementMetadata(prediction)$prob[hits_enh@to])

  bins <- as.data.frame(table(cres_EP$between))
  cres_EP1 <- cres_EP[cres_EP$EP_reg_distance > 0.5, ]

  if (nrow(cres_EP1) != 0) {
    bins_pos <- as.data.frame(table(cres_EP1$between))
    all_bins <- merge(bins, bins_pos, by.x = "Var1", by.y = "Var1", all.x = TRUE)
    ## cases in which bins_pos is 0 will have an NA value which will be 0 in the
    ## next step
    all_bins[is.na(all_bins)] <- 0
    colnames(all_bins) <- c("pair", "bins", "bins_pos")
    input <- merge(input, all_bins, by.x = "pair", by.y = "pair")
    input$reg_dist_enh <- input$bins_pos
    input$norm_reg_dist_enh <- input$bins_pos / input$bins

  } else {
    #avoid uncommon edge case in which all ET pairs have EP_reg_distance above 0.5 is
    # 0
    input$reg_dist_enh <- 0
    input$norm_reg_dist_enh <- 0
  }
  return(input)
}

###############################################################################
# function: get PP scores for region between enhancer and promoter
###############################################################################

compute_crup_reg_distance_prom <- function(input, prediction, between_ranges) {

  hits_prom <- GenomicRanges::findOverlaps(between_ranges, prediction)
  cres_PP <- data.frame(between = GenomicRanges::elementMetadata(between_ranges)$pair[hits_prom@from],
                        PP_reg_distance = GenomicRanges::elementMetadata(prediction)$probP[hits_prom@to])

  bins <- as.data.frame(table(cres_PP$between))
  cres_PP1 <- cres_PP[cres_PP$PP_reg_distance > 0.5, ]

  if (nrow(cres_PP1) != 0) {
    bins_pos <- as.data.frame(table(cres_PP1$between))
    all_bins <- merge(bins, bins_pos, by.x = "Var1", by.y = "Var1", all.x = TRUE)
    ## cases in which bins_pos is 0 will have an NA value which will be 0 in the
    ## next step
    all_bins[is.na(all_bins)] <- 0
    colnames(all_bins) <- c("pair", "bins", "bins_pos")
    input <- merge(input, all_bins, by.x = "pair", by.y = "pair")
    input$reg_dist_prom <- input$bins_pos.y
    input$norm_reg_dist_prom <- input$bins_pos.y / input$bins.y

  } else {
    #avoid uncommon edge case in which all ET pairs have EP_reg_distance above 0.5 is
    # 0
    input$reg_dist_prom <- 0
    input$norm_reg_dist_prom <- 0
  }
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
  x <- merge(x, tpmfile[, c(3, 4)], by.x = "gene_id2", by.y = "gene_id2")
  return(x)
}

################################################################################
# function: create a region annotation dataframe for the compute cell type
# features function
################################################################################

createRegionsDf <- function(listProm, listEnh, pairs) {

  conn <- RSQLite::dbConnect(RSQLite::SQLite(),
                             system.file("extdata",
                                         "Annotation.db",
                                         package = "CENTRE"))
  #get chromosome tts new_start and new_end of input genes
  query <- paste("SELECT  gene_id1, chr, transcription_start, new_start, new_end FROM gencode WHERE gene_id1 in (",
                 paste0(sprintf("'%s'", listProm$gene_id2), collapse = ", "),
                 ")", sep = "")
  regionsProm <- RSQLite::dbGetQuery(conn, query)
  #get chr middle new_start new_end point of input enhancers
  queryEnh <-  paste("SELECT  V5, V1, middle_point, new_start, new_end FROM ccres_enhancer WHERE V5 in (",
                     paste0(sprintf("'%s'", listEnh$enhancer_id),
                            collapse = ", "),
                     ")", sep = "")
  regionsEnhancer <- RSQLite::dbGetQuery(conn, queryEnh)
  RSQLite::dbDisconnect(conn)

  ##create a dataframe with the middle point newstart and newend for each of the
  ##pairs
  regions <- merge(pairs,
                   regionsEnhancer,
                   by.x = "enhancer_id",
                   by.y = "V5",
                   all.x = TRUE)

  regions <- merge(regions,
                   regionsProm,
                   by.x = "gene_id2",
                   by.y = "gene_id1",
                   all.x = TRUE)

  #add distance value to sort start and end for regulatory distance calculations.
  regions$distance <- regions$middle_point - regions$transcription_start
  return(regions)
}

################################################################################
# function : helper function to download the files needed after install
################################################################################

downloader <- function(file, method) {
  url <- paste0("http://owww.molgen.mpg.de/~CENTRE_data/", file)
  cat(paste0("Downloading ", file, "\n"))
  exit <- download.file(url,
                        destfile = paste(system.file("extdata",
                                                     package = "CENTRE")
                                       , file,
                                       sep = "/") ,
                        method = method)
  if (exit != 0) {
    stop(paste0("Download of ",
                file, " failed. Non-zero exit status."))
  }

  f <- system.file("extdata", file, package = "CENTRE")
  if (!file.exists(f)) {
    stop(paste0("Download of ", file,
                " failed or file was saved in the wrong directory."))
  }
}

