#' Compute cell type specific features
#'
#' Computes the cell type specific features needed for the CENTRE classification
#' step.
#'
#' @param metaData Data.frame indicating the paths to the ChIP-seq experiments.
#' More information on the format here `crupR::normalize
#' @param replicate The number of replicates of the ChIP-seq experiments
#' that need to be normalized.
#' @param input.free Boolean value indicating whether a Control/Input ChIP-seq
#' experiment is provided to go with the Histone Modification ChIP-seq
#' experiments.
#' If the parameter is set to FALSE the normalization of ChIP-seq experiments
#' will be run in input.free mode.
#' @param cores Number of cores to compute the CRUP score features
#' @param sequencing Type of sequencing of the ChIP-seq experiments "paired" or
#' "single". The parameter takes single as default
#' @param tpmData Dataframe of two columns one with the RNA-seq TPM values,
#' one with the names of the genes given as ENSEMBLE ID's
#' @param chr NULL or a vector of chromosomes. Use only if the crupR
#' normalization should be done for certain chromosomes. If this parameter is
#' not used crupR normalization is done for all chromosomes.
#' Using only certain chromosomes for normalization might change results
#' and is not the intented used of crupR or CENTRE.
#' @param pairs The output of `CENTRE::createPairs()`.
#'
#'
#' @return
#' A table containing the following computed features :
#'* CRUP enhancer score for enhancer region, promoter region and the region
#' between the enhancer and the promoter.
#'* CRUP promoter score for enhancer region, promoter region and the region
#' between the enhancer and the promoter.
#'* TPM values from the RNA-seq experiment given.
#'
#'
#' @examples
#' #for the sake of runtime we load pairs that
#' #were precomputed for the example. The user should run createPairs()
#' #to get them.
#' \dontrun{
#' pairs <- readRDS(file = system.file("extdata",
#'        "input_cellType_pairs.rds",
#'        package = "CENTRE"
#'    ))
#'
#' # Compute Cell-type features
#' eh <- ExperimentHub::ExperimentHub()
#'
#' files <- c(
#'     system.file("extdata",
#'         "example/HeLa_H3K4me1.REF_chr19_reduced.bam",
#'         package = "CENTRE"
#'     ),
#'     system.file("extdata",
#'         "example/HeLa_H3K4me3.REF_chr19_reduced.bam",
#'         package = "CENTRE"
#'     ),
#'     system.file("extdata",
#'         "example/HeLa_H3K4me3.REF_chr19_reduced.bam",
#'         package = "CENTRE"
#'     )
#' )
#' inputs <- system.file("extdata",
#'     "example/HeLa_input.REF_chr19_reduced.bam",
#'     package = "CENTRE"
#' )
#' metaData <- data.frame(
#'     HM = c("H3K4me1", "H3K4me3", "H3K27ac"),
#'     condition = c(1, 1, 1), replicate = c(1, 1, 1),
#'     bamFile = files, inputFile = rep(inputs, 3)
#' )
#'
#' tpmpath <- unname(eh[["EH9545"]])
#' tpmfile <- read.table(tpmpath,
#'     sep = "", stringsAsFactors = FALSE,
#'     header = TRUE
#' )
#' tpmfile <- tpmfile[grep("E", tpmfile$gene_id), ]
#' celltype_features <- CENTRE::computeCellTypeFeatures(metaData,
#'     replicate = 1,
#'     input.free = FALSE,
#'     cores = 1,
#'     chr = "chr19",
#'     sequencing = "single",
#'     tpmData = tpmfile,
#'     pairs = pairs
#' )
#'}
#' @export
#' @importFrom crupR normalize getEnhancers
#' @import utils
#' @importFrom GenomicRanges GRanges findOverlaps elementMetadata
#' @importFrom IRanges IRanges
#' @importFrom stats reshape
#' @importFrom AnnotationHub AnnotationHub
#' @importFrom CENTREannotation fetch_data
#' @importClassesFrom CENTREannotation CENTREannotDb
#' @importFrom dplyr left_join join_by inner_join %>% select .data
#' @importFrom BiocParallel MulticoreParam
computeCellTypeFeatures <- function(metaData,
    replicate,
    input.free = FALSE,
    cores,
    sequencing = "single",
    tpmData,
    chr = NULL,
    pairs) {
    startTime <- Sys.time()
    ## Computing the crup scores
    message("Computing CRUP score features...\n")

    if (missing(pairs)) {
        stop("Need to provide a dataframe of enhancer and gene pairs.")
    }

    needed_names <- c("gene_id1", "enhancer_id")
    if (!all(needed_names %in% colnames(pairs))) {
        missing_cols <- setdiff(needed_names, colnames(pairs))
        stop(
            "Error: The following expected columns are missing: ",
            paste(missing_cols, collapse = ", ")
        )
    }

    if (missing(tpmData)) {
        stop("Need to provide a dataframe of RNA-seq TPM values.")
    }
    ## Calling normalization step only on the chromosomes we have
    normalized <- crupR::normalize(
        metaData = metaData,
        condition = 1,
        replicate = replicate,
        mapq = 10,
        input.free = input.free,
        genome = "hg38",
        sequencing = sequencing,
        chroms = chr,
        BPPARAM = BiocParallel::MulticoreParam(workers = cores)
    )
    # Get CRUP enhancer probabilities
    crupScores <- crupR::getEnhancers(data = normalized, all = TRUE)
    listEnh <- unique(pairs$enhancer_id)
    listProm <- unique(pairs$gene_id1)
    # Get Gencode and CCRes anntotations for the input genes and enhancers
    regions <- createRegionsDf(listProm, listEnh, pairs)
    pairs$pair <- paste(pairs$enhancer_id, pairs$gene_id1, sep = "_")

    message(paste0(
        "Getting the CRUP-EP scores for enhancer, promoter and the",
        "\n",
        "regulatory distance...\n"
    ))

    crupEPFeatures <- getEPFeatures(regions, crupScores, pairs)


    message(paste0(
        "Getting the CRUP-PP scores for enhancer, promoter and the",
        "\n",
        "regulatory distance...\n"
    ))

    # Crup promoter scores for distance
    # Compute the promoter probability from probA and probE
    # In CRUP probA is the probability of a region being an active reg. element
    # probE is the probability of a region being an active enhancer
    crupScores$probP <- crupScores$probA * (1 - crupScores$probE)
    crupFeatures <- getPPFeatures(regions, crupScores, crupEPFeatures)

    message("Getting the TPM values.\n")
    featuresTableAll <- getRNAseq(crupFeatures, tpmData)

    featuresTableAll <- reformatDf(featuresTableAll)

    message(paste0("time: ", format(Sys.time() - startTime), "\n"))
    return(featuresTableAll)
}
