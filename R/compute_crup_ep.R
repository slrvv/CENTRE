###############################################################################
# function: get enhancer scores for enhancers
###############################################################################
compute_crup_EP_enhancer <- function(input, crup_scores) {
  print("Computing CRUP EP scores for enhancer...")
  print(input)
  ccres_enhancer <- as.data.frame(ccres_enhancer)

  # This should go outside so you dont compute the regions twice
  regions <- subset(ccres_enhancer,
                    ccres_enhancer$V5 %in% input$enhancer_id,
                    select = c(V1, V5, new_start, new_end))

  genes_ranges <- with(regions,
                       GenomicRanges::GRanges(V1,
                                              IRanges::IRanges(start = new_start,
                                                               end = new_end)))
  #make this into one function for crup EP and PP

  hits_crup <- GenomicRanges::findOverlaps(genes_ranges, crup_scores)
  cres_EP <- data.frame(cres = hits_crup@from, EP = hits_crup@to)
  cres_EP$cres_name <- input$enhancer_id[cres_EP$cres]
  cres_EP$EP_enhancer <- GenomicRanges::elementMetadata(crup_scores)$score[cres_EP$EP]
  cres_EP$bin<-rep(1:5,nrow(regions))
  trial<-reshape(cres_EP[,3:5],
                 direction="wide",
                 timevar = "bin",
                 idvar="cres_name",
                 v.names="EP_enhancer")
  return(trial)



}

###############################################################################
# function: get enhancer scores for promoters
###############################################################################
compute_crup_EP_promoter <- function(input, crup_scores) {
  print("Computing CRUP EP scores for promoter...")
  gencode <- as.data.frame(gencode)

  regions <- subset(gencode,
                    gencode$gene_id %in% input[, 1],
                    select = c(chr, gene_id, new_start, new_end))

  genes_ranges <- with(regions,
                       GenomicRanges::GRanges(chr,
                                              IRanges::IRanges(start = new_start,
                                                               end = new_end)))

  hits_crup <- GenomicRanges::findOverlaps(genes_ranges, crup_scores)
  cres_EP <- data.frame(promoter = hits_crup@from, EP = hits_crup@to)
  cres_EP$gene_name <- input$gene_id[cres_EP$promoter]
  cres_EP$EP_promoter <- GenomicRanges::elementMetadata(crup_scores)$score[cres_EP$EP]

  cres_EP$bin<-rep(1:5,nrow(regions))
  trial<-reshape(cres_EP[,3:5],
                 idvar="gene_name",
                 timevar = "bin",
                 direction="wide",
                 v.names="EP_promoter")
  return(trial)


}

###############################################################################
# function: get enhancer scores for distance between enhancer and promoter
###############################################################################

compute_crup_EP_reg_distance <- function(input, prediction) {
  print("Computing CRUP EP scores for regulatory distance...")
  ccres_enhancer <- unfactorize(ccres_enhancer)
  gencode <- unfactorize(gencode)
  middle_point <- ccres_enhancer[ccres_enhancer$V5 %in% input$enhancer_id, ]$middle_point
  tts <- gencode[gencode$gene_id %in% input$gene_id, ]$transcription_start

  input$chr <- gencode[gencode$gene_id %in% input$gene_id, ]$chr

  input$distance <- middle_point - tts
  input$bstart <- middle_point
  input$bstart[input$distance > 0] <- tts[input$distance > 0]
  input$bend <- tts
  input$bend[input$distance > 0] <- middle_point[input$distance > 0]


  input$pair <- paste0(input$gene_id, "_", input$enhancer_id)#173509
  input <- input[!(duplicated(input$pair)), ]#169049
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
# function: get promoter scores for enhancer
###############################################################################
compute_crup_PP_enhancer <- function(input, crup_scores) {
  print("Computing CRUP PP scores for enhancer...")
  ccres_enhancer <- as.data.frame(ccres_enhancer)

  regions <- subset(ccres_enhancer, ccres_enhancer$V5 %in% input$enhancer_id,
                    select = c(V1, V5, new_start, new_end))

  genes_ranges <- with(regions, GenomicRanges::GRanges(V1, IRanges::IRanges(start = new_start, end = new_end)))

  hits_crup <- GenomicRanges::findOverlaps(genes_ranges, crup_scores)
  cres_EP <- data.frame(cres = hits_crup@from, EP = hits_crup@to)
  cres_EP$cres_name <- input$enhancer_id[cres_EP$cres]
  cres_EP$PP_enhancer <- GenomicRanges::elementMetadata(crup_scores)$score[cres_EP$EP]
  cres_EP$bin<-rep(1:5,nrow(regions))
  trial<-reshape(cres_EP[,3:5],idvar="cres_name",timevar = "bin",direction="wide",v.names="PP_enhancer")
  return(trial)


}

###############################################################################
# function: get promoter scores for promoter
###############################################################################

compute_crup_PP_promoter <- function(input, crup_scores) {
  print("Computing CRUP PP scores for promoter...")
  gencode <- as.data.frame(gencode)

  regions <- subset(gencode, gencode$gene_id %in% input[, 1],
                    select = c(chr, gene_id, new_start, new_end))

  genes_ranges <- with(regions, GenomicRanges::GRanges(chr, IRanges::IRanges(start = new_start, end = new_end)))

  hits_crup <- GenomicRanges::findOverlaps(genes_ranges, crup_scores)
  cres_EP <- data.frame(promoter = hits_crup@from, EP = hits_crup@to)
  cres_EP$gene_name <- input$gene_id[cres_EP$promoter]
  cres_EP$PP_promoter<- GenomicRanges::elementMetadata(crup_scores)$score[cres_EP$EP]
  cres_EP$bin<-rep(1:5,nrow(regions))
  trial<-reshape(cres_EP[,3:5],idvar="gene_name",timevar = "bin",direction="wide",v.names="PP_promoter")
  return(trial)

}

###############################################################################
# function: get promoter scores for regulatory distance
###############################################################################
compute_crup_PP_reg_distance <- function(input, prediction) {
  print("Computing CRUP PP scores for regulatory distance...")
  ccres_enhancer <- unfactorize(ccres_enhancer)
  gencode <- unfactorize(gencode)
  middle_point <- ccres_enhancer[ccres_enhancer$V5 %in% input$enhancer_id, ]$middle_point
  tts <- gencode[gencode$gene_id %in% input$gene_id, ]$transcription_start

  input$chr <- gencode[gencode$gene_id %in% input$gene_id, ]$chr

  input$distance <- middle_point - tts
  input$bstart <- middle_point
  input$bstart[input$distance > 0] <- tts[input$distance > 0]
  input$bend <- tts
  input$bend[input$distance > 0] <- middle_point[input$distance > 0]


  input$pair <- paste0(input$gene_id, "_", input$enhancer_id)#173509
  input <- input[!(duplicated(input$pair)), ]#169049

  between_ranges <- with(input, GenomicRanges::GRanges(chr,IRanges::IRanges(start=bstart, end=bend)))
  hits_enh <- GenomicRanges::findOverlaps(between_ranges, prediction)

  cres_EP <- data.frame(between = hits_enh@from,
              PP_reg_distance = GenomicRanges::elementMetadata(prediction)$score[hits_enh@to])

  bins<-as.data.frame(table(cres_EP$between))



  cres_EP1<-cres_EP[cres_EP$EP_prob>0.5,]

  if (nrow(cres_EP1) == 0){

    bins_pos <- as.data.frame(matrix(c(seq(1:nrow(input)), rep(0, times = nrow(input))),
                                     nrow=nrow(input),
                                     ncol= 2))
    colnames(bins_pos) <- c("Var1", "Freq")
    print(bins_pos)
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
# function: compute features
###############################################################################

compute_features_crup_scores <- function(input, crupinput, cores, chromosome){
  normalized <- crupR::normalize(metaData = crupinput, condition = 1, replicate = 1,
                                 genome = "hg38", sequencing = "single",
                                 chroms = chromosome , cores = cores) ##how do we deal with the type of sequencing?

  crup_scores_enh <- crupR::getEnhancers(data = normalized, cores = cores)

  crup_scores_enh <- crup_scores_enh$data_matrix

  crup_scores_prom <- crupR::getEnhancers(data = normalized, cores = cores, promprob = T)
  crup_scores_prom <- crup_scores_prom$data_matrix

  crup_EP_enh <- compute_crup_EP_enhancer(input, crup_scores_enh)
  crup_EP_prom <- compute_crup_EP_promoter(input, crup_scores_enh)
  crup_EP_reg <- compute_crup_EP_reg_distance(input, crup_scores_enh)
  crup_PP_enh <- compute_crup_PP_enhancer(input, crup_scores_prom)
  crup_PP_prom <- compute_crup_PP_promoter(input, crup_scores_prom)
  crup_PP_reg <- compute_crup_PP_reg_distance(input, crup_scores_prom)

  input <- cbind(input, crup_EP_enh[2:6])
  input <- cbind(input, crup_EP_prom[2:6])
  input <- cbind(input, crup_EP_reg)
  input <- cbind(input, crup_PP_enh[2:6])
  input <- cbind(input, crup_PP_prom[2:6])
  input <- cbind(input, crup_PP_reg)

  return(input)

}
