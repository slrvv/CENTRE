get_data <- function() {
  gencode <-  read.table(
    "/project/CRUP_scores/bed_files_V38/gencode_v38.one.transcript.bed",
    header = T, sep = "\t")
  
  ccres <- read.table(
    "/project/CRUP_scores/bed_files_V38/GRCh38-cCREsV3_all500.bed",
    header = T, sep = "\t")
  
  ccres_enhancer <- ccres[ccres$V6
                          %in% c("dELS", "pELS", "pELS,CTCF-bound",
                                 "dELS,CTCF-bound"), ]
  chromosomes <- read.table("/project/CRUP_scores/Scripts/EPI/chromosomes.txt",
                            sep = "\t")
  usethis::use_data(gencode, ccres, ccres_enhancer, chromosomes, internal= T, overwrite = T)
}
