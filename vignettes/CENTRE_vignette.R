## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE,
                      comment = "#>")

## ---- eval=FALSE--------------------------------------------------------------
#  CENTRE::downloadPrecomputedData(method = "wget", all = FALSE)
#  # Make sure whatever method you use to download is available on your system

## -----------------------------------------------------------------------------
genes <- as.data.frame(c("ENSG00000130203.10",
                         "ENSG00000171119.3"))
colnames(genes) <- c("gene_id") #It is important to name the column gene_id
pairs <- CENTRE::createPairs(genes)

## -----------------------------------------------------------------------------
colnames(pairs) <- c("gene_id", "enhancer_id")
generic_features <- CENTRE::computeGenericFeatures(pairs)

