#helper functions of the computeCellTypeFeatures function


#'@description Create a data.frame containing the positional information of 
#'enhancers and genes in the pairs data.frame outputted by createPairs()
#'
#'@param listProm list of unique promoters in pairs dataframe 
#'@param listEnh list of unique enhancers in pairs dataframe
#'@param pairs output of createPairs
#'@noRd
createRegionsDf <- function(listProm, listEnh, pairs) {
  
  ah <- AnnotationHub::AnnotationHub()
  CENTREannotgeneDb <- ah[["AH116730"]]
  CENTREannotenhDb <- ah[["AH116731"]]
  
  regionsProm <- CENTREannotation::fetch_data(CENTREannotgeneDb,
                                              columns = c("gene_id1", 
                                                          "chr", 
                                                          "transcription_start", 
                                                          "new_start", 
                                                          "new_end"), 
                                              entries = listProm, 
                                              column_filter = "gene_id1")
  
  #get chr middle new_start new_end point of input enhancers
  regionsEnhancer <- CENTREannotation::fetch_data(CENTREannotenhDb,
                                                  columns = c("enhancer_id", 
                                                              "chr", 
                                                              "middle_point", 
                                                              "new_start", 
                                                              "new_end"), 
                                                  entries = listEnh, 
                                                  column_filter = "enhancer_id")
  
  ##create a dataframe with the middle point newstart and newend for each of the
  ##pairs
  regions <- dplyr::left_join(x = pairs,
                              y = regionsEnhancer,
                              by = dplyr::join_by(enhancer_id))
  
  regions <- dplyr::left_join(x = regions,
                              y = regionsProm,
                              by = dplyr::join_by(gene_id1),
                              suffix = c(".enh", ".gene"))
  
  #add distance value to sort start and end for regulatory distance calculations.
  regions$distance <- regions$middle_point - regions$transcription_start
  return(regions)
}

#'@description Get the CRUP scores (Enhancer Probability or Promoter probability)
#'for each enhancer.
#'
#'@param regionsEnhancer data.frame output of createRegionsDf
#'@param crupScores GRanges containing the CRUP scores. 
#'Output of crupR::getEnhancers()
#'@param promprob TRUE if the CRUP scores are Promoter probabilities (CRUP-PP),
#'@noRd
computeCrupEnhancer <- function(regionsEnhancer,
                                crupScores,
                                promprob = FALSE) {
  
  #Overlapping the  enhancer ranges with the crup scores

  enhancerRanges <- with(regionsEnhancer,
                         GenomicRanges::GRanges(chr.enh,
                                                IRanges::IRanges(start = new_start.enh,
                                                                 end = new_end.enh),
                                                enhancer_id = enhancer_id))
  enhancerRanges <- unique(enhancerRanges)
  hitsCrup <- GenomicRanges::findOverlaps(enhancerRanges, crupScores)
  cresEP <- data.frame(cres = hitsCrup@from, EP = hitsCrup@to)
  cresEP$enhancer_id <- GenomicRanges::elementMetadata(enhancerRanges)$enhancer_id[cresEP$cres]
  
  if (promprob == TRUE) {
    # get the PP crup scores
    cresEP$PP_prob_enh <- GenomicRanges::elementMetadata(crupScores)$probP[cresEP$EP]
    cresEP$enhancer_id <- factor(cresEP$enhancer_id)
    cresEP$bin <- rep(1:5, times = length(enhancerRanges))
    
    trial <- stats::reshape(cresEP[, 3:5],
                            idvar = "enhancer_id",
                            timevar = "bin",
                            direction = "wide",
                            v.names = "PP_prob_enh")
  } else {
    #get the EP crup scores
    cresEP$EP_prob_enh <- GenomicRanges::elementMetadata(crupScores)$prob[cresEP$EP]
    cresEP$enhancer_id <- factor(cresEP$enhancer_id)
    cresEP$bin <- rep(1:5, times = length(enhancerRanges))
    
    trial <- stats::reshape(cresEP[, 3:5],
                            idvar = "enhancer_id",
                            timevar = "bin",
                            direction = "wide",
                            v.names = "EP_prob_enh")
  }
  
  return(trial)
}

#'@description Get the CRUP scores (Enhancer Probability or Promoter probability)
#'for each promoter.
#'
#'@param regionsEnhancer data.frame output of createRegionsDf.
#'@param crupScores GRanges containing the CRUP scores. 
#'Output of crupR::getEnhancers()
#'@param promprob TRUE if the CRUP scores are Promoter probabilities (CRUP-PP),
#'@noRd
computeCrupPromoter <- function(regionsProm,
                                crupScores,
                                promprob = FALSE) {
  #Overlapping with CRUP scores
  geneRanges <- with(regionsProm,
                     GenomicRanges::GRanges(chr.gene,
                                            IRanges::IRanges(start = new_start.gene,
                                                             end = new_end.gene),
                                            gene_id1 = gene_id1))
  
  geneRanges <- unique(geneRanges)
  hitsCrup <- GenomicRanges::findOverlaps(geneRanges, crupScores)
  cresEP <- data.frame(promoter = hitsCrup@from, EP = hitsCrup@to)
  cresEP$gene_id1 <- GenomicRanges::elementMetadata(geneRanges)$gene_id1[cresEP$promoter]
  
  if (promprob == TRUE) {
    cresEP$PP_prob_gene <- GenomicRanges::elementMetadata(crupScores)$probP[cresEP$EP]
    cresEP$gene_id1 <- factor(cresEP$gene_id1)
    cresEP$bin <- rep(1:5, length(geneRanges))
    
    trial <- stats::reshape(cresEP[, 3:5],
                            idvar = "gene_id1",
                            timevar = "bin",
                            direction = "wide",
                            v.names = "PP_prob_gene")
  } else {
    cresEP$EP_prob_gene <- GenomicRanges::elementMetadata(crupScores)$prob[cresEP$EP]
    cresEP$gene_id1 <- factor(cresEP$gene_id1)
    cresEP$bin <- rep(1:5, length(geneRanges))
    
    trial <- stats::reshape(cresEP[, 3:5],
                            idvar = "gene_id1",
                            timevar = "bin",
                            direction = "wide",
                            v.names = "EP_prob_gene")
  }
  
  return(trial)
}

#'@description create GRanges for the regions between enhancer and gene in 
#'each of the pairs.
#'
#'@param regionsEnhancer data.frame output of createRegionsDf.
#'@noRd
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
  regions$chr <- regions$chr.enh
  #Make the gene enhancer pairs into ranges
  regions$pair <- paste(regions$enhancer_id, regions$gene_id1, sep = "_")
  betweenRanges <- with(regions,
                        GenomicRanges::GRanges(chr,
                                               IRanges::IRanges(start = bstart,
                                                                end = bend),
                                               pair = pair))
  betweenRanges <- unique(betweenRanges)

  return(betweenRanges)
}


#'@description compute CRUP regulatory distance based on enhancer probabilities.
#'
#'@param input data.frame with all the annotation on the gene-enhancer pairs
#'@param prediction GRanges of CRUP scores
#'@param betweenRanges GRanges of the regions between each enhancer-gene pair
#'@noRd
computeCrupRegDistanceEnh <- function(input, prediction, betweenRanges) {
  ##overlap the ranges objects with predictions
  hitsEnh <- GenomicRanges::findOverlaps(betweenRanges, prediction)
  cresEP <- data.frame(between = GenomicRanges::elementMetadata(betweenRanges)$pair[hitsEnh@from],
                        EP_reg_distance = GenomicRanges::elementMetadata(prediction)$prob[hitsEnh@to])
  
  bins <- as.data.frame(table(cresEP$between))
  cresEP1 <- cresEP[cresEP$EP_reg_distance > 0.5, ]
  
  if (nrow(cresEP1) != 0) {
    binsPos <- as.data.frame(table(cresEP1$between))
    allBins <- dplyr::left_join(bins,
                                binsPos,
                                dplyr::join_by(Var1)) %>% replace(is.na(.), 0)

    ## cases in which bins_pos is 0 will have an NA value which will be 0 in the
    ## next step
    colnames(allBins) <- c("pair", "bins", "bins_pos")
    input <- dplyr::inner_join(input, 
                               allBins,
                               dplyr::join_by(pair))
    input$reg_dist_enh <- input$bins_pos
    input$norm_reg_dist_enh <- input$bins_pos / input$bins
    
  } else {
    #avoid uncommon edge case in which all ET pairs have EP_reg_distance above 0.5 is
    # 0
    input$reg_dist_enh <- 0
    input$bins <- 0
    input$norm_reg_dist_enh <- 0
    input$bins_pos <- 0
  }
  return(input)
}

#'@description compute CRUP regulatory distance based on promoter probabilities.
#'
#'@param input data.frame with all the annotation on the gene-enhancer pairs
#'@param prediction GRanges of CRUP scores with the Promoter Probabilities
#'@param betweenRanges GRanges of the regions between each enhancer-gene pair
#'@noRd
computeCrupRegDistanceProm <- function(input, prediction, betweenRanges) {
  
  hitsProm <- GenomicRanges::findOverlaps(betweenRanges, prediction)
  cresPP <- data.frame(between = GenomicRanges::elementMetadata(betweenRanges)$pair[hitsProm@from],
                        PP_reg_distance = GenomicRanges::elementMetadata(prediction)$probP[hitsProm@to])
  
  bins <- as.data.frame(table(cresPP$between))
  cresPP1 <- cresPP[cresPP$PP_reg_distance > 0.5, ]
  if (nrow(cresPP1) != 0) {
    binsPos <- as.data.frame(table(cresPP1$between))
    allBins <- dplyr::left_join(bins,
                                binsPos,
                                dplyr::join_by(Var1)) %>% replace(is.na(.), 0)
    ## cases in which bins_pos is 0 will have an NA value which will be 0 in the
    ## next step
    colnames(allBins) <- c("pair", "bins", "bins_pos")
    input <- dplyr::inner_join(input, 
                               allBins,
                               dplyr::join_by(pair))
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

#'@description collect RNA-seq TPM from the user provided data and add it to
#'the data.frame with all of the features
#'
#'@param x data.frame collecting the computed features
#'@param tpmfile data.frame of TPM RNA-seq values given by the user
#'@noRd
getRNAseq <- function(x, tpmfile) {
  tpmfile$gene_id1 <- gsub("\\..*", "", tpmfile[, 1])
  x <- dplyr::left_join(x,
                        tpmfile[, c(3, 4)], 
                        dplyr::join_by(gene_id1))
  return(x)
}

#'@description collect and join to main data.frame all the features based on 
#'Enhancer probabilities
#'
#'@param regions createRegionsDf
#'@param crupScores GRanges of CRUP scores
#'@noRd
getEPFeatures <- function(regions, crupScores, pairs){
  #Crup enhancer scores for enhancer
  crupEPenh <- computeCrupEnhancer(regions,
                                   crupScores)
  crupEPFeaturesEnh <- dplyr::left_join(pairs,
                                        crupEPenh,
                                        by = dplyr::join_by(enhancer_id))
  
  #CRUP enhancer scores for promoter
  crupEPprom <- computeCrupPromoter(regions,
                                    crupScores)
  crupEPFeatures <- dplyr::left_join(crupEPFeaturesEnh,
                                     crupEPprom,
                                     by = dplyr::join_by(gene_id1))
  
  #create the betweenRanges objects that is used for the distance calculations
  betweenRanges <- createBetweenRanges(regions)
  crupEPFeatures <- computeCrupRegDistanceEnh(crupEPFeatures,
                                                  crupScores,
                                                  betweenRanges)
  return(crupEPFeatures)
  
}


#'@description collect and join to main data.frame all the features based on 
#'Promoter probabilities
#'
#'@param regions createRegionsDf
#'@param crupScores GRanges of CRUP scores with Promoter Probabilities
#'@noRd
getPPFeatures <- function(regions, crupScores, crupEPFeatures){
  #Crup enhancer scores for enhancer
  crupPPenh <- computeCrupEnhancer(regions,
                                   crupScores, 
                                   promprob = TRUE)
  
  crupFeatures <- dplyr::left_join(crupEPFeatures,
                                   crupPPenh,
                                   by = dplyr::join_by(enhancer_id))
  
  #CRUP enhancer scores for promoter
  crupPPprom <- computeCrupPromoter(regions,
                                    crupScores, 
                                    promprob = TRUE)
  
  crupFeatures <- dplyr::left_join(crupFeatures,
                                   crupPPprom,
                                   by = dplyr::join_by(gene_id1))
  
  #create the betweenRanges objects that is used for the distance calculations
  betweenRanges <- createBetweenRanges(regions)
  crupFeatures <- computeCrupRegDistanceProm(crupFeatures,
                                                crupScores,
                                                betweenRanges)
  return(crupFeatures)
  
}


#'@description reformat table with all cell type features
#'
#'@param featuresDf Table with all cell type features
#'@noRd
reformatDf <- function(featuresDf){
  featuresDfReformat <- featuresDf %>% select(gene_id1,
                                                       enhancer_id,
                                                       EP_prob_enh.1,
                                                       EP_prob_enh.2,
                                                       EP_prob_enh.3,
                                                       EP_prob_enh.4,
                                                       EP_prob_enh.5,
                                                       EP_prob_gene.1,
                                                       EP_prob_gene.2,
                                                       EP_prob_gene.3,
                                                       EP_prob_gene.4,
                                                       EP_prob_gene.5,
                                                       reg_dist_enh,
                                                       norm_reg_dist_enh,
                                                       PP_prob_enh.1,
                                                       PP_prob_enh.2,
                                                       PP_prob_enh.3,
                                                       PP_prob_enh.4,
                                                       PP_prob_enh.5,
                                                       PP_prob_gene.1,
                                                       PP_prob_gene.2,
                                                       PP_prob_gene.3,
                                                       PP_prob_gene.4,
                                                       PP_prob_gene.5,
                                                       reg_dist_prom,
                                                       norm_reg_dist_prom,
                                                       TPM) %>% replace(is.na(.), 0)

  return(featuresDfReformat)
}
