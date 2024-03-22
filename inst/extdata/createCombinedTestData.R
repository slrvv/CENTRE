# ------------------------------------------------------------------------------
# Script Name: Creating combined test data base
# Purpose: Creates the combinedTest.db database, previously
# the combined pvalues were computed during the running time
# of CENTRE
# Author: Sara Lopez Ruiz de Vargas
# Date: 2023-12-05
# Version: 1.0
# ------------------------------------------------------------------------------

# necessary functions

################################################################################
# function: Makes the query for Precomputed.db
################################################################################

queryMaker <- function(table, feature, x) {
  query <- paste("SELECT ", feature, ", pair FROM ", table, " WHERE pair in (",
                 paste0(sprintf("'%s'", x$pair), collapse = ", "), ")", sep = "")
  return(query)
}

################################################################################
# function: gets the precomputed values of the Wilcoxon tests from
# PrecomputedData.db
################################################################################

getPrecomputedValues <- function(table, feature, x, conn) {

  query <- queryMaker(table, feature, x)
  df_return <- RSQLite::dbGetQuery(conn, query)
  rownames(df_return) <- df_return$pair
  df_return$pair <- NULL
  print(head(df_return))
  return(df_return)
}


# First you need to create all posible pairs based on the GENCODE database

devtools::load_all() ## use createPairs function from here

conn <- RSQLite::dbConnect(RSQLite::SQLite(),
                           system.file("extdata",
                                       "Annotation.db",
                                       package = "CENTRE"))

query <- "SELECT gene_id FROM gencode"
gene <- RSQLite::dbGetQuery(conn, query)
benchmark <- readRDS(file= system.file("extdata",
                                       "input_generic_features.rds",
                                       package = "CENTRE"))

colnames(gene) <- c("gene_id")

RSQLite::dbDisconnect(conn)

pairs <- createPairs(gene)
#dim(pairs)
#Checking that we dont have pair duplicates
cat(paste0("All pairs are unique ", nrow(pairs) == nrow(unique(pairs)) ))

## Get the values of each test for all pairs
pairs$pair <- paste(pairs$enhancer_id,
                                pairs$gene_id,
                                sep = "_")

conn <- RSQLite::dbConnect(RSQLite::SQLite(),
                           system.file("extdata",
                                       "PrecomputedData.db",
                                       package = "CENTRE"))

cageDf <- getPrecomputedValues("cage_test_data",
                               "wilcoxtest_cage",
                               pairs,
                               conn)

dhsexpDf <- getPrecomputedValues("dhsexp_test_data",
                                 "wilcoxtest_dhs_exp",
                                 pairs,
                                 conn)
crupexpDf <- getPrecomputedValues("crupexp_test_data",
                                  "wilcoxtest_crup_exp",
                                  pairs,
                                  conn)
dhsdhsDf <- getPrecomputedValues("dhsdhs_test_data",
                                 "wilcoxtest_dhs_dhs",
                                 pairs,
                                 conn)

RSQLite::dbDisconnect(conn)


### remove the rownames and make it a variable

cageDf$pair <- rownames(cageDf)
rownames(cageDf) <- NULL

pValuedata <- merge(pairs, cageDf, by.x = "pair", by.y="pair", all.x = T )

dhsexpDf$pair <- rownames(dhsexpDf)
rownames(dhsexpDf) <- NULL


pValuedata <- merge(pValuedata, dhsexpDf, by.x = "pair", by.y= "pair", all.x = T )

crupexpDf$pair <- rownames(crupexpDf)
rownames(crupexpDf) <- NULL

pValuedata <- merge(pValuedata, crupexpDf, by.x = "pair", by.y= "pair", all.x = T )


dhsdhsDf$pair <- rownames(dhsdhsDf)
rownames(dhsdhsDf) <- NULL


pValuedata <- merge(pValuedata, dhsdhsDf, by.x = "pair", by.y= "pair", all.x = T )

pValuedata[pValuedata$pair=="EH38E3350767_ENSG00000059728",]
#check pValues are the correct ones for this particular pair

pValue <- list(pValuedata$wilcoxtest_cage, pValuedata$wilcoxtest_dhs_exp,
               pValuedata$wilcoxtest_crup_exp, pValuedata$wilcoxtest_dhs_dhs)

## Combining the values of the Wilcoxon tests
combinedtests <- metapod::combineParallelPValues(pValue, method="fisher")



pValuedata$combined_tests <- -log(combinedtests$p.value) ## at this point there is a mixing up of the ids



rownames(pValuedata) <- pValuedata$pair

pValuedata <- pValuedata[, - 1]

saveRDS(pValuedata,
        file = "/project/CRUP_scores/CENTRE/inst/extdata/combinedTestData.rds" )
