# *CENTRE - Short description*
CENTRE is a tool for identifying and characterizing gene enhancer interactions.


### Contact

lopez_s@molgen.mpg.de

### Citation


## Requirements
- R (tested 4.0.0)
- crupR
- GenomicRanges and IRanges
- scran

## User Provided Data

CENTRE computes features based on Histone ChIP-seq (H3K27ac, H3K4me3 and H3K4me1
) data and RNA-seq data that the user provides. As well as, a file with the 
genes of interest or the genes and enhancer pairs of interest. An example for 
all of these files in provided inside the package.

### User data : 

- Histone ChIP-sep in BAM format for H3K27ac, H3K4me3 and H3K4me1
- RNA-seq ENCODE4 gene quantifications in TSV format. This file will have two 
columns one 
  with the ENSEMBL ID's and one with the TPM values. 
- Plain text file with the genes or gene and enhancer pairs of interest with 
  columns separated by tabs. 


## CENTRE Provided Data

CENTRE uses precomputed datasets that the user needs to download from the 
following link (insert url) and add to the /data folder. 


## How to run CENTRE

### Installing CENTRE
Install this R package. It is not yet published so we recommend blabla

### General workflow

Two use cases : 
`
- User has only genes: `createPairs()` -> `computeGenericFeatures()` -> 
`computeCellTypeFeatures()` -> `centreClassification()`
- User has gene enhancer pairs : `computeGenericFeatures()` -> 
`computeCellTypeFeatures()` -> `centreClassification()`

### Create Pairs
The function `createPairs()` finds all enhancers in 500KB from the transcription
start site of the genes provided and creates all possible enhancer gene pairs.
For more information on how this function runs check the manuals in the /man folder.


### Compute generic features
The function `computeGenericFeatures()` computes the features that are not cell 
type specific. For more information on how this function runs check the manuals
in the /man folder.

### Compute cell type specific features
The function `computeCellTypeFeatures()` computes the features that are not cell 
type specific. For more information on how this function runs check the manuals
in the /man folder.


### Example run


```



```

