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
) data and RNA-seq data that the user provides. As well as, a file with the genes
of interest or the genes and enhancer pairs of interest. An example for all of these
files in provided inside the package.

###User data : 

- Histone ChIP-sep in BAM format for H3K27ac, H3K4me3 and H3K4me1
- RNA-seq ENCODE4 gene quantifications in TSV format. This file will have two columns one 
  with the ENSEMBL ID's and one with the TPM values. 
- Plain text file with the genes or gene and enhancer pairs of interest with 
  columns separated by tabs. 


## CENTRE Provided Data

CENTRE uses precomputed datasets that the user needs to download from the following
link (insert url) and add to the /data folder. 


##How to run CENTRE

### General workflow

1. Install this R package. It is not yet published so we recommend blabla
2. Add the data into the /data folder
3. Provide the experiments that 
3. Compute the features
4. Make the predictions

### Compute features
The function `compute_features()` computes the features needed for CENTRE to 
make predictions. For more information on how this function runs check the manuals
in the /man folder.

### Make predictions


### Example run
