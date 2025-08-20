# helper functions of main package functions except for computeCellTypeFeatures

#' @description Create the candidate enhancer gene pairs at 500KB from the input
#' genes. Helper function of createPairs()
#'
#' @param gene vector of ENSEMBL ID's given by the user
#' @noRd
geneCenteredPairs <- function(gene) {
    ## remove the "." version id of ENSEMBL ids
    gene$gene_id1 <- gsub("\\..*", "", gene$gene_id)


    ## connect to our GENCODE v40 database to get tts of the genes

    message("Getting gene annotation.\n")
    ah <- AnnotationHub::AnnotationHub()

    CENTREannotgeneDb <- ah[["AH116730"]]

    # get chromosome and tts of our genes
    gene <- CENTREannotation::fetch_data(CENTREannotgeneDb,
        columns = c(
            "gene_id1",
            "chr",
            "transcription_start"
        ),
        entries = gene$gene_id1,
        column_filter = "gene_id1"
    )

    genesRange <- with(
        gene,
        GenomicRanges::GRanges(chr,
            IRanges::IRanges(
                start = transcription_start,
                end = transcription_start
            ),
            gene_id1 = gene_id1
        )
    )

    # extend the gene region 500Kb to the left of TTS and to the right
    genesRange <- regioneR::extendRegions(genesRange,
        extend.start = 500000,
        extend.end = 500000
    )

    # Select all of the annotation for ccres v3

    message("Getting enhancer annotation.\n")
    CENTREannotenhDb <- ah[["AH116731"]]

    # to make the ranges that is overlapped smaller retrieve only the genes in
    # the same list of chromosomes as we have in the input gene dataframe

    chrList <- unique(gene$chr)
    ccresEnhancer <- CENTREannotation::fetch_data(CENTREannotenhDb,
        columns = c(
            "enhancer_id",
            "chr",
            "new_start",
            "new_end"
        ),
        entries = chrList,
        column_filter = "chr"
    )


    enhancerRange <- with(
        ccresEnhancer,
        GenomicRanges::GRanges(chr,
            IRanges::IRanges(
                start = new_start,
                end = new_end
            ),
            enhancer_id = enhancer_id
        )
    )


    # find the enhancers that overlap the extended gene region
    overlaps <- GenomicRanges::findOverlaps(genesRange, enhancerRange,
        ignore.strand = TRUE
    )

    geneMeta <- GenomicRanges::elementMetadata(genesRange)
    enhMeta <- GenomicRanges::elementMetadata(enhancerRange)
    ccresOverlapping <- data.frame(
        gene_id1 = geneMeta$gene_id1[overlaps@from],
        enhancer_id = enhMeta$enhancer_id[overlaps@to]
    )

    return(ccresOverlapping)
}



#' @description Create the candidate enhancer gene pairs at 500KB from the input
#' enhancers. Helper function of createPairs()
#'
#' @param enhancer vector of ENCODE cCREs enhancer ID's gievn bu the user
#' @noRd
enhancerCenteredPairs <- function(enhancer) {
    # get chromosome and middle point of our enhancers
    message("Getting enhancer annotation.\n")
    ah <- AnnotationHub::AnnotationHub()

    CENTREannotenhDb <- ah[["AH116731"]]

    enhancer <- CENTREannotation::fetch_data(CENTREannotenhDb,
        columns = c(
            "enhancer_id",
            "chr",
            "middle_point"
        ),
        entries = enhancer$enhancer_id,
        column_filter = "enhancer_id"
    )

    enhancerRange <- with(
        enhancer,
        GenomicRanges::GRanges(chr,
            IRanges::IRanges(
                start = middle_point,
                end = middle_point
            ),
            enhancer_id = enhancer_id
        )
    )

    # extend the enhancer region 500Kb to the left of middle point and to the right
    enhancerRange <- regioneR::extendRegions(enhancerRange,
        extend.start = 500000,
        extend.end = 500000
    )

    # Select all of the annotation from gencode
    message("Getting gene annotation.\n")
    CENTREannotgeneDb <- ah[["AH116730"]]

    # to make the ranges that is overlapped smaller retrieve only the genes in
    # the same list of chromosomes as we have in the input enhancer dataframe

    chrList <- unique(enhancer$chr)
    gene <- CENTREannotation::fetch_data(CENTREannotgeneDb,
        columns = c(
            "gene_id",
            "chr",
            "new_start",
            "new_end"
        ),
        entries = chrList,
        column_filter = "chr"
    )

    gene$gene_id1 <- gsub("\\..*", "", gene$gene_id)



    genesRange <- with(
        gene,
        GenomicRanges::GRanges(chr,
            IRanges::IRanges(
                start = new_start,
                end = new_end
            ),
            gene_id1 = gene_id1
        )
    )


    # find the enhancers that overlap the extended gene region
    overlaps <- GenomicRanges::findOverlaps(enhancerRange,
        genesRange,
        ignore.strand = TRUE
    )
    geneMeta <- GenomicRanges::elementMetadata(genesRange)
    enhMeta <- GenomicRanges::elementMetadata(enhancerRange)
    ccresOverlapping <- data.frame(
        gene_id1 = geneMeta$gene_id1[overlaps@to],
        enhancer_id = enhMeta$enhancer_id[overlaps@from]
    )
    return(ccresOverlapping)
}



#' @description Get the precomputed values from the ExperimentHub database
#'
#' @param table Which of the tables in the precomputedTests database to fetch
#' from.
#' @param feature Column in the table to fetch from.
#' @param x data.frame containing column of pair identifiers.
#' @noRd
getPrecomputedValues <- function(table, feature, x) {
    eh <- ExperimentHub::ExperimentHub()
    # connect to precomputed data through experimentHub
    suppressMessages(precompDb <- eh[["EH9540"]])
    # fetch the needed data from the database using CENTREprecomputed package
    dfReturn <- CENTREprecomputed::fetch_data(precompDb,
        table = table,
        columns = c("pair", feature),
        entries = x$pair,
        column_filter = "pair"
    )
    return(dfReturn)
}


#' @description Compute distance between enhancer middle point and TTS,
#' Helper function of computeGenericFeatures()
#'
#' @param x data.frame containing the enhancer gene pairs, with columns gene_id1
#' and enhancer_id
#' @noRd
computeDistances <- function(x) {
    # connect to annotation dataBase
    ah <- AnnotationHub::AnnotationHub()
    CENTREannotgeneDb <- ah[["AH116730"]]
    CENTREannotenhDb <- ah[["AH116731"]]
    # get chromosome and tts of our genes
    gencode <- CENTREannotation::fetch_data(CENTREannotgeneDb,
        columns = c(
            "gene_id1",
            "chr",
            "transcription_start"
        ),
        entries = x$gene_id1,
        column_filter = "gene_id1"
    )
    # get chr and middle point of enhancers
    ccres_enhancer <- CENTREannotation::fetch_data(CENTREannotenhDb,
        columns = c(
            "enhancer_id",
            "chr",
            "middle_point"
        ),
        entries = x$enhancer_id,
        column_filter = "enhancer_id"
    )

    # Get the chr gene_id and transcription_start from gencode annotation
    # Getting the chrosomes and the middle points for the provided enhancers
    result <- dplyr::inner_join(
        x = x,
        y = ccres_enhancer[, c(
            "chr",
            "enhancer_id",
            "middle_point"
        )],
        by = dplyr::join_by(enhancer_id)
    )

    # Getting the chrosomes and transcription start sites for the provided genes
    result <- dplyr::inner_join(
        x = result,
        y = gencode[, c(
            "chr",
            "gene_id1",
            "transcription_start"
        )],
        by = dplyr::join_by(gene_id1),
        suffix = c(".enh", ".gene")
    )


    message("Removing all pairs that are not in the same chromosome.\n")
    result <- result[(result$chr.enh == result$chr.gene), ]
    result$distance <- abs(result$middle_point - result$transcription_start)
    return(result)
}

#' @description Get the precomputed values from the ExperimentHub database
#'
#' @param table Which of the tables in the precomputedTests database to fetch
#' from.
#' @param feature Column in the table to fetch from.
#' @param x data.frame containing column of pair identifiers.
#' @noRd
getPrecomputedValues <- function(table, feature, x) {
    eh <- ExperimentHub::ExperimentHub()
    # connect to precomputed data through experimentHub
    suppressMessages(precompDb <- eh[["EH9540"]])
    # fetch the needed data from the database using CENTREprecomputed package
    dfReturn <- CENTREprecomputed::fetch_data(precompDb,
        table = table,
        columns = c("pair", feature),
        entries = x$pair,
        column_filter = "pair"
    )
    return(dfReturn)
}


#' @description Compute distance between enhancer middle point and TTS,
#' Helper function of computeGenericFeatures()
#'
#' @param x data.frame containing the enhancer gene pairs, with columns gene_id1
#' and enhancer_id
#' @noRd
computeDistances <- function(x) {
    # connect to annotation dataBase
    ah <- AnnotationHub::AnnotationHub()
    CENTREannotgeneDb <- ah[["AH116730"]]
    CENTREannotenhDb <- ah[["AH116731"]]
    # get chromosome and tts of our genes
    gencode <- CENTREannotation::fetch_data(CENTREannotgeneDb,
        columns = c(
            "gene_id1",
            "chr",
            "transcription_start"
        ),
        entries = x$gene_id1,
        column_filter = "gene_id1"
    )
    # get chr and middle point of enhancers
    ccres_enhancer <- CENTREannotation::fetch_data(CENTREannotenhDb,
        columns = c(
            "enhancer_id",
            "chr",
            "middle_point"
        ),
        entries = x$enhancer_id,
        column_filter = "enhancer_id"
    )

    # Get the chr gene_id and transcription_start from gencode annotation
    # Getting the chrosomes and the middle points for the provided enhancers
    result <- dplyr::inner_join(
        x = x,
        y = ccres_enhancer[, c(
            "chr",
            "enhancer_id",
            "middle_point"
        )],
        by = dplyr::join_by(enhancer_id)
    )

    # Getting the chrosomes and transcription start sites for the provided genes
    result <- dplyr::inner_join(
        x = result,
        y = gencode[, c(
            "chr",
            "gene_id1",
            "transcription_start"
        )],
        by = dplyr::join_by(gene_id1),
        suffix = c(".enh", ".gene")
    )


    message("Removing all pairs that are not in the same chromosome.\n")
    result <- result[(result$chr.enh == result$chr.gene), ]
    result$distance <- abs(result$middle_point - result$transcription_start)
    return(result)
}
