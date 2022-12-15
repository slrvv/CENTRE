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
distances_gene_enhancer <- function(x) {


  #Getting the chrosomes and the middle points for the provided enhancers
  result <- merge(x,
                  ccres_enhancer[,c('V1','V5','middle_point')],
                  by.x='enhancer_id',
                  by.y='V5') #change V1 and V5 to more meaningful names
  #Getting the chrosomes and transcription start sites for the provided genes
  result <- merge(result,
                  gencode[,c('chr','gene_id1','transcription_start')],
                  by.x='gene_id1',
                  by.y= 'gene_id1')

  cat("Removing all gene enhancer pairs that are not in the same chromosome.\n")

  result <- result[(result$V1 == result$chr),]

  result$distance <- result$middle_point - result$transcription_start


  return(result)
}

###############################################################################
# function: extend genes by 50000bp and map the gene region
# to overlapping enhancers
###############################################################################




###############################################################################
# function: get scores for enhancers
###############################################################################
compute_crup_enhancer <- function(input, crup_scores) {

  #Getting the extended regions for the provided enhancers and making it into a
  #GRanges object
  regions <-  merge(input,
                    ccres_enhancer[,c('V1', 'V5', 'new_start', 'new_end')],
                    by.x='enhancer_id',
                    by.y='V5')


  enhancer_ranges <- with(regions,
                       GenomicRanges::GRanges(V1,
                                              IRanges::IRanges(start = new_start,
                                                               end = new_end)))
  #Overlapping the  enhancer ranges with the crup scores
  hits_crup <- GenomicRanges::findOverlaps(enhancer_ranges, crup_scores)

  cres_EP <- data.frame(cres = hits_crup@from, EP = hits_crup@to)
  cres_EP$cres_name <- input$enhancer_id[cres_EP$cres]
  cres_EP$EP_enhancer <- GenomicRanges::elementMetadata(crup_scores)$score[cres_EP$EP]##change to score

  #Returning the probabilities in bins
  cres_EP$bin<-rep(1:5,nrow(regions))
  print(cres_EP)
  trial<-reshape(cres_EP[,3:5],
                 direction="wide",
                 timevar = "bin",
                 idvar="cres_name",
                 v.names="EP_enhancer")
  return(trial)

}

###############################################################################
# function: get  scores for promoters
###############################################################################
compute_crup_promoter <- function(input, crup_scores) {
  #Getting the extended regions for the provided promoters and making it into a
  #GRanges object
  regions <- merge(input,
                  gencode[,c('chr','gene_id1','new_start', 'new_end')],
                  by.x='gene_id1',
                  by.y='gene_id1')

  #Overlapping the  promoter ranges with the crup scores
  genes_ranges <- with(regions,
                       GenomicRanges::GRanges(chr,
                                              IRanges::IRanges(start = new_start,
                                                               end = new_end)))

  hits_crup <- GenomicRanges::findOverlaps(genes_ranges, crup_scores)
  cres_EP <- data.frame(promoter = hits_crup@from, EP = hits_crup@to)
  cres_EP$gene_name <- input$gene_id[cres_EP$promoter]
  cres_EP$EP_promoter <- GenomicRanges::elementMetadata(crup_scores)$score[cres_EP$EP]

  #Returning the probabilities in bins
  cres_EP$bin<-rep(1:5,nrow(regions))
  trial<-reshape(cres_EP[,3:5],
                 idvar="gene_name",
                 timevar = "bin",
                 direction="wide",
                 v.names="EP_promoter")
  return(trial)
}

###############################################################################
# function: get scores for distance between enhancer and promoter
###############################################################################

compute_crup_reg_distance <- function(input, prediction) {

  ##Check if the distances are negative and flip the start and end around
  print(input)
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
                        EP_reg_distance = GenomicRanges::elementMetadata(prediction)$score[hits_enh@to])

  bins<-as.data.frame(table(cres_EP$between))


  cres_EP1<-cres_EP[cres_EP$EP_prob>0.5,]

  if (nrow(cres_EP1) == 0){

    bins_pos <- as.data.frame(matrix(c(c(seq(1:nrow(input)), rep(0, times = nrow(input))),
                                       nrow=nrow(input),
                                       ncol= 2)))
    colnames(bins_pos) <- c("Var1", "Freq")
  } else {
    bins_pos<-as.data.frame(table(cres_EP1$between))
  }


  all_bins<-merge(bins,bins_pos,by.x="Var1",by.y="Var1",all.x=TRUE)
  all_bins[is.na(all_bins)] <- 0

  colnames(all_bins)<-c("pair","bins","bins_pos")

  reg_dist_enh<- all_bins$bins_pos
  norm_reg_dist_enh<-all_bins$bins_pos/all_bins$bins

  reg_distance <- cbind(reg_dist_enh, norm_reg_dist_enh)


  return(reg_distance)



}



###############################################################################
# function: compute CRUP features
###############################################################################

compute_features_crup_scores <- function(input, crupinput, cores, chromosome){
  print(chromosome)
  normalized <- crupR::normalize(metaData = crupinput, condition = 1, replicate = 1,
                                 genome = "hg38", sequencing = "single",
                                 chroms = chromosome , cores = cores) ##how do we deal with the type of sequencing?

  crup_scores_enh <- crupR::getEnhancers(data = normalized, cores = cores)

  crup_scores_enh <- crup_scores_enh$data_matrix

  crup_scores_prom <- crupR::getEnhancers(data = normalized, cores = cores, promprob = T)
  crup_scores_prom <- crup_scores_prom$data_matrix

  crup_EP_enh <- compute_crup_enhancer(input, crup_scores_enh)
  crup_EP_prom <- compute_crup_promoter(input, crup_scores_enh)
  crup_EP_reg <- compute_crup_reg_distance(input, crup_scores_enh)
  crup_PP_enh <- compute_crup_enhancer(input, crup_scores_prom)
  crup_PP_prom <- compute_crup_promoter(input, crup_scores_prom)
  crup_PP_reg <- compute_crup_reg_distance(input, crup_scores_prom)

  input <- cbind(input, crup_EP_enh[2:6])
  input <- cbind(input, crup_EP_prom[2:6])
  input <- cbind(input, crup_EP_reg)
  input <- cbind(input, crup_PP_enh[2:6])
  input <- cbind(input, crup_PP_prom[2:6])
  input <- cbind(input, crup_PP_reg)

  return(input)

}

################################################################################
# function: gets the precomputed values of the Wilcoxon tests
################################################################################

wilcoxon_test_features <- function(x){

  x$cage_wilcoxon_test <- cage_test_data[x$pair,3]
  x$dhsexp_wilcoxon_test  <-  dhsexp_test_data[x$pair,3]
  x$crupexp_wilcoxon_test  <-  crupexp_test_data[x$pair,3]
  x$dhsdhs_wilcoxon_test  <- dhsdhs_test_data[x$pair,3]
  return(x)
}



################################################################################
# function: gets the precomputed values of the CRUP correlation for our input
# gene enhancer pairs
################################################################################

crup_correlations <- function(x){
  x$cor_CRUP <- crup_cor[x$pair,3]
  return(x)

}


################################################################################
# function: get the RNA seq TPM values for our genes
################################################################################

get_rnaseq <- function(x, tpmfile){
  tpmfile$gene_id2 <- gsub("\\..*","",tpmfile[,1])
  rnaseq_tpmvalue <- subset(tpmfile, tpmfile$gene_id2 %in% x$gene_id2,
                            select = V3)
  return(rnaseq_tpmvalue)
}
