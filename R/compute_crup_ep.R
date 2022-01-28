
compute_crup_EP_enhancer <- function(input, crupinput, cores){

  normalized <- crupR::normalize(metaData = crupinput, condition = 1, replicate = 1,
                                 genome = "hg38", sequencing = "single",
                                 chroms = c("chr1"), cores = cores) ##how do we deal with the type of sequencing?
  crup_scores <- crupR::getEnhancers(data = normalized, cores = cores)

  crup_scores <- crup_scores$data_matrix

  ccres_enhancer <- as.data.frame(ccres_enhancer)

  regions <- subset(ccres_enhancer, ccres_enhancer$V5 %in% input$V2, select = c(V1,V5,new_start, new_end))

  genes_ranges <- with(regions, GenomicRanges::GRanges(V1, IRanges::IRanges(start=new_start,end=new_end)))

  hits_crup<- GenomicRanges::findOverlaps(genes_ranges,crup_scores)
  cres_EP<-data.frame(cres=hits_crup@from,EP=hits_crup@to)
  cres_EP$cres_name <-ccres_enhancer$V5[cres_EP$cres]
  cres_EP$EP_prob<-GenomicRanges::elementMetadata(crup_scores)$score[cres_EP$EP]
  normilized_score<-aggregate(EP_prob ~ cres_name, cres_EP, sum)

  return(normilized_score)



}

