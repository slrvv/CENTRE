###############################################################################
# function: get enhancer scores for enhancers
###############################################################################
compute_crup_EP_enhancer <- function(input, crup_scores) {

  ccres_enhancer <- as.data.frame(ccres_enhancer)

  regions <- subset(ccres_enhancer, ccres_enhancer$V5 %in% input$V2,
                    select = c(V1, V5, new_start, new_end))

  genes_ranges <- with(regions, GenomicRanges::GRanges(
    V1, IRanges::IRanges(start = new_start, end = new_end)))

  hits_crup <- GenomicRanges::findOverlaps(genes_ranges, crup_scores)
  cres_EP <- data.frame(cres = hits_crup@from, EP = hits_crup@to)
  cres_EP$cres_name <- ccres_enhancer$V5[cres_EP$cres]
  cres_EP$EP_prob <- GenomicRanges::elementMetadata(crup_scores)$score[cres_EP$EP]
  normilized_score <- aggregate(EP_prob ~ cres_name, cres_EP, sum)

  return(normilized_score)



}

###############################################################################
# function: get enhancer scores for promoters
###############################################################################
compute_crup_EP_promoter <- function(input, crup_scores) {

  gencode <- as.data.frame(gencode)

  regions <- subset(gencode, gencode$gene_id %in% input[, 1],
                    select = c(chr, gene_id, new_start, new_end))

  genes_ranges <- with(regions, GenomicRanges::GRanges(
    chr, IRanges::IRanges(start = new_start, end = new_end)))

  hits_crup <- GenomicRanges::findOverlaps(genes_ranges, crup_scores)
  cres_EP <- data.frame(promoter = hits_crup@from, EP = hits_crup@to)
  cres_EP$gene_name <- gencode$gene_id[cres_EP$promoter]
  cres_EP$EP_prob <- GenomicRanges::elementMetadata(crup_scores)$score[cres_EP$EP]
  normilized_score  <- aggregate(EP_prob ~ gene_name, cres_EP, sum)

  return(normilized_score)


}

###############################################################################
# function: get enhancer scores for distance between enhancer and promoter
###############################################################################

compute_crup_EP_reg_distance <- function(input, prediction) {
  ccres_enhancer <- unfactorize(ccres_enhancer)
  gencode <- unfactorize(gencode)
  middle_point <- ccres_enhancer[ccres_enhancer$V5 %in% input$V2, ]$middle_point
  tts <- gencode[gencode$gene_id %in% input$V1, ]$transcription_start

  input$chr <- gencode[gencode$gene_id %in% input$V1, ]$chr

  input$distance <- middle_point - tts
  input$bstart <- middle_point
  input$bstart[input$distance > 0] <- tts[input$distance > 0]
  input$bend <- tts
  input$bend[input$distance > 0] <- middle_point[input$distance > 0]


  input$pair <- paste0(input$V1, "_", input$V2)#173509
  input <- input[!(duplicated(input$pair)), ]#169049

  between_ranges <- with(input, GenomicRanges::GRanges(
    chr, IRanges::IRanges(start = bstart, end = bend)))
  hits_enh <- GenomicRanges::findOverlaps(between_ranges, prediction)
  cres_EP <- data.frame(between = hits_enh@from,
                      EP_prob = GenomicRanges::elementMetadata(prediction)$score[hits_enh@to])
  normilized_score <- aggregate(EP_prob ~ between, cres_EP, mean)
  return(normilized_score)



}
###############################################################################
# function: get promoter scores for enhancer
###############################################################################
compute_crup_PP_enhancer <- function(input, crup_scores) {

  ccres_enhancer <- as.data.frame(ccres_enhancer)

  regions <- subset(ccres_enhancer, ccres_enhancer$V5 %in% input$V2,
                    select = c(V1, V5, new_start, new_end))

  genes_ranges <- with(regions, GenomicRanges::GRanges(
    V1, IRanges::IRanges(start = new_start, end = new_end)))

  hits_crup <- GenomicRanges::findOverlaps(genes_ranges, crup_scores)
  cres_EP <- data.frame(cres = hits_crup@from, EP = hits_crup@to)
  cres_EP$cres_name <- ccres_enhancer$V5[cres_EP$cres]
  cres_EP$EP_prob <- GenomicRanges::elementMetadata(crup_scores)$score[cres_EP$EP]
  normilized_score <- aggregate(EP_prob ~ cres_name, cres_EP, sum)

  return(normilized_score)



}

###############################################################################
# function: get promoter scores for promoter
###############################################################################

compute_crup_PP_promoter <- function(input, crup_scores) {
  gencode <- as.data.frame(gencode)

  regions <- subset(gencode, gencode$gene_id %in% input[, 1],
                    select = c(chr, gene_id, new_start, new_end))

  genes_ranges <- with(regions, GenomicRanges::GRanges(
    chr, IRanges::IRanges(start = new_start, end = new_end)))

  hits_crup <- GenomicRanges::findOverlaps(genes_ranges, crup_scores)
  cres_EP <- data.frame(promoter = hits_crup@from, EP = hits_crup@to)
  cres_EP$gene_name <- gencode$gene_id[cres_EP$promoter]
  cres_EP$EP_prob <- GenomicRanges::elementMetadata(crup_scores)$score[cres_EP$EP]
  normilized_score <- aggregate(EP_prob ~ gene_name, cres_EP, sum)

  return(normilized_score)


}

###############################################################################
# function: get promoter scores for regulatory distance
###############################################################################
compute_crup_PP_reg_distance <- function(input, prediction) {
  ccres_enhancer <- unfactorize(ccres_enhancer)
  gencode <- unfactorize(gencode)
  middle_point <- ccres_enhancer[ccres_enhancer$V5 %in% input$V2, ]$middle_point
  tts <- gencode[gencode$gene_id %in% input$V1, ]$transcription_start

  input$chr <- gencode[gencode$gene_id %in% input$V1, ]$chr

  input$distance <- middle_point - tts
  input$bstart <- middle_point
  input$bstart[input$distance > 0] <- tts[input$distance > 0]
  input$bend <- tts
  input$bend[input$distance > 0] <- middle_point[input$distance > 0]


  input$pair <- paste0(input$V1, "_", input$V2)#173509
  input <- input[!(duplicated(input$pair)), ]#169049

  between_ranges <- with(input, GenomicRanges::GRanges(
    chr,IRanges::IRanges(start=bstart, end=bend)))
  hits_enh <- GenomicRanges::findOverlaps(between_ranges, prediction)

  cres_EP <- data.frame(between = hits_enh@from,
              EP_prob = GenomicRanges::elementMetadata(prediction)$score[hits_enh@to])

  normilized_score <- aggregate(EP_prob ~ between, cres_EP, mean)
  return(normilized_score)

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

  crup_EP_enh <- compute_crup_EP_enhancer(input, crupinput, cores, crup_scores_enh)
  crup_EP_prom <- compute_crup_EP_promoter(input, crupinput, cores, crup_scores_enh)
  crup_EP_reg <- compute_crup_EP_reg_distance(input, crupinput, cores, crup_scores_enh)
  crup_PP_enh <- compute_crup_EP_enhancer(input, crupinput, cores, crup_scores_prom)
  crup_PP_prom <- compute_crup_EP_promoter(input, crupinput, cores, crup_scores_prom)
  crup_PP_reg <- compute_crup_EP_reg_distance(input, crupinput, cores, crup_scores_prom)

  input$EP_enh <- crup_EP_enh$EP_prob
  input$EP_prom <- crup_EP_prom$EP_prob
  input$EP_distance <- crup_EP_reg$EP_prob
  input$PP_enh <- crup_PP_enh$EP_prob
  input$PP_prom <- crup_PP_prom$EP_prob
  input$PP_distance <- crup_PP_reg$EP_prob

  return(input)

}
