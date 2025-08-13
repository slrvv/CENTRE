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
# function: get the distance from gene to enhancer
###############################################################################
computeDistances <- function(x) {
  # connect to annotation dataBase
  ah <- AnnotationHub::AnnotationHub()
  CENTREannotgeneDb <- ah[["AH116730"]]
  CENTREannotenhDb <- ah[["AH116731"]]
  #get chromosome and tts of our genes
  gencode <- CENTREannotation::fetch_data(CENTREannotgeneDb,
                                       columns = c("gene_id1", 
                                                   "chr", 
                                                   "transcription_start"), 
                                       entries = x$gene_id1, 
                                       column_filter = "gene_id1")
  #get chr and middle point of enhancers
  ccres_enhancer <- CENTREannotation::fetch_data(CENTREannotenhDb, 
                                                columns = c("enhancer_id", 
                                                   "chr", 
                                                   "middle_point"), 
                                                entries = x$enhancer_id, 
                                                column_filter = "enhancer_id")
  
  #Get the chr gene_id and transcription_start from gencode annotation
  #Getting the chrosomes and the middle points for the provided enhancers
  result <- dplyr::inner_join(x = x,
                  y = ccres_enhancer[, c("chr", "enhancer_id", "middle_point")],
                  by = dplyr::join_by("enhancer_id" == "enhancer_id")) 

  #Getting the chrosomes and transcription start sites for the provided genes
  result <- dplyr::inner_join(x = result,
                  y = gencode[, c("chr", "gene_id1", "transcription_start")],
                  by = dplyr::join_by("gene_id1" == "gene_id1"),
                  suffix = c(".enh", ".gene")) 
  

  message("Removing all gene enhancer pairs that are not in the same chromosome.\n")
  result <- result[(result$chr.enh == result$chr.gene), ]
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
# function: gets the precomputed values of the Wilcoxon tests from
# PrecomputedData.db
################################################################################

getPrecomputedValues <- function(table, feature, x) {

  eh <- ExperimentHub::ExperimentHub()
  precompDb <- eh[["EH9540"]]
  df_return <- CENTREprecomputed::fetch_data(precompDb, 
                                            table = table, 
                                            columns = c("pair", feature), 
                                            entries = x$pair, 
                                            column_filter = "pair")
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

  ah <- AnnotationHub::AnnotationHub
  CENTREannotgeneDb <- ah[["AH116730"]]
  CENTREannotenhDb <- ah[["AH116731"]]

  regionsProm <- CENTREannotation::fetch_data(CENTREannotgeneDb,
                                       columns = c("gene_id1", 
                                                   "chr", 
                                                   "transcription_start", 
                                                   "new_start", 
                                                   "new_end"), 
                                       entries = listProm$gene_id2, 
                                       column_filter = "gene_id1")
  #get chr middle new_start new_end point of input enhancers
  queryEnh <-  paste("SELECT  V5, V1, middle_point, new_start, new_end FROM ccres_enhancer WHERE V5 in (",
                     paste0(sprintf("'%s'", listEnh$enhancer_id),
                            collapse = ", "),
                     ")", sep = "")
  regionsEnhancer <- CENTREannotation::fetch_data(CENTREannotenhDb,
                                       columns = c("enhancer_id", 
                                                   "chr", 
                                                   "middle_point", 
                                                   "new_start", 
                                                   "new_end"), 
                                       entries = listEnh$enhancer_id, 
                                       column_filter = "enhancer_id")

  ##create a dataframe with the middle point newstart and newend for each of the
  ##pairs
  regions <- dplyr::left_join(x = pairs,
                   y = regionsEnhancer,
                   by = c("enhancer_id" == "enhancer_id"))

  regions <- dplyr::left_join(x = regions,
                   y = regionsProm,
                   by = c("gene_id2" == "gene_id1"),
                   suffix = c(".enh", ".gene"))

  #add distance value to sort start and end for regulatory distance calculations.
  regions$distance <- regions$middle_point - regions$transcription_start
  return(regions)
}



################################################################################
# function : gene centered pairs
################################################################################

geneCenteredPairs <- function(gene){
  
  ## remove the "." version id of ENSEMBL ids
  gene$gene_id1 <- gsub("\\..*", "", gene$gene_id)
  
  
  ## connect to our GENCODE v40 database to get tts of the genes
  
  ah <- AnnotationHub::AnnotationHub()
  
  CENTREannotgeneDb <- ah[["AH116730"]]
  
  #get chromosome and tts of our genes
  gene <- CENTREannotation::fetch_data(CENTREannotgeneDb,
                                       columns = c("gene_id1", 
                                                   "chr", 
                                                   "transcription_start"), 
                                       entries = gene$gene_id1, 
                                       column_filter = "gene_id1")
  
  genesRange <- with(gene,
                     GenomicRanges::GRanges(chr,
                                            IRanges::IRanges(start = transcription_start,
                                                             end = transcription_start),
                                            gene_id1 = gene_id1))
  
  #extend the gene region 500Kb to the left of TTS and to the right
  genesRange <- regioneR::extendRegions(genesRange,
                                        extend.start = 500000,
                                        extend.end = 500000)
  
  #Select all of the annotation for ccres v3
  
  CENTREannotenhDb <- ah[["AH116731"]]
  
  #to make the ranges that is overlapped smaller retrieve only the genes in 
  #the same list of chromosomes as we have in the input gene dataframe
  
  chrList <- unique(gene$chr)
  ccresEnhancer <- CENTREannotation::fetch_data(CENTREannotenhDb,
                                       columns = c("enhancer_id", 
                                                   "chr", 
                                                   "new_start", 
                                                   "new_end"), 
                                       entries = chrList, 
                                       column_filter = "chr")
 
  
  enhancerRange <-  with(ccresEnhancer,
                         GenomicRanges::GRanges(chr,
                                                IRanges::IRanges(start = new_start,
                                                                 end = new_end),
                                                enhancer_id = enhancer_id))
  
  
  # find the enhancers that overlap the extended gene region
  overlaps <- GenomicRanges::findOverlaps(genesRange, enhancerRange,
                                          ignore.strand = TRUE)
  
  ccresOverlapping <- data.frame(gene_id1 = GenomicRanges::elementMetadata(genesRange)$gene_id1[overlaps@from],
                                 enhancer_id = GenomicRanges::elementMetadata(enhancerRange)$enhancer_id[overlaps@to])
  
  return(ccresOverlapping)
  
}



################################################################################
# function : enhancer centered pairs
################################################################################

enhancerCenteredPairs <- function(enhancer){
  
  #get chromosome and middle point of our enhancers
  
  ah <- AnnotationHub::AnnotationHub()

  CENTREannotenhDb <- ah[["AH116731"]]
  
  enhancer <- CENTREannotation::fetch_data(CENTREannotenhDb,
                                           columns = c("enhancer_id", "chr", "middle_point"),
                                           entries =  enhancer$enhancer_id,
                                           column_filter = "enhancer_id")

  enhancerRange <- with(enhancer,
                        GenomicRanges::GRanges(chr,
                                               IRanges::IRanges(start = middle_point,
                                                                end = middle_point),
                                               enhancer_id = enhancer_id))
  
  #extend the enhancer region 500Kb to the left of middle point and to the right
  enhancerRange <- regioneR::extendRegions(enhancerRange,
                                           extend.start = 500000,
                                           extend.end = 500000)
  
  #Select all of the annotation from gencode
  CENTREannotgeneDb <- ah[["AH116730"]]
  
  #to make the ranges that is overlapped smaller retrieve only the genes in 
  #the same list of chromosomes as we have in the input enhancer dataframe
  
  chrList <- unique(enhancer$chr)
  gene <- CENTREannotation::fetch_data(CENTREannotgeneDb,
                                       columns = c("gene_id", 
                                                   "chr", 
                                                   "new_start", 
                                                   "new_end"), 
                                       entries = chrList, 
                                       column_filter = "chr")
  
  gene$gene_id1 <- gsub("\\..*", "", gene$gene_id)
  
  

  genesRange <-  with(gene,
                         GenomicRanges::GRanges(chr,
                                                IRanges::IRanges(start = new_start,
                                                                 end = new_end),
                                                gene_id1 = gene_id1))
  
  
  # find the enhancers that overlap the extended gene region
  overlaps <- GenomicRanges::findOverlaps(enhancerRange,
                                          genesRange,
                                          ignore.strand = TRUE)
  
  ccresOverlapping <- data.frame(gene_id1 = GenomicRanges::elementMetadata(genesRange)$gene_id1[overlaps@to],
                        enhancer_id = GenomicRanges::elementMetadata(enhancerRange)$enhancer_id[overlaps@from])
  return(ccresOverlapping)
}
