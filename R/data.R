#' CAGE test data
#'
#' Result of the Wilcoxon tests on the CAGE data of the enhancer gene pairs.
#'
#' @format A data frame with three variables: \code{gene_id} Gene ID,
#' \code{symbol38} Enhancer ID, \code{wilcoxtest_cage} result of the test ).
"cage_test_data"

#' CRUP correlations
#'
#' Correlation between the CRUP scores of enhancer and genes  in the pairs
#'
#' @format A data frame with three variables: \code{gene_id} Gene ID,
#' \code{symbol38} Enhancer ID, \code{cor_CRUP} value of the correlation ).
"crup_cor"

#' CRUP expression test
#'
#' Results of the Wilcoxon test on the CRUP scores and RNA expression of the pairs
#'
#' @format A data frame with three variables: \code{gene_id} Gene ID,
#' \code{symbol38} Enhancer ID, \code{wilcoxtest_crup_exp} result of the test ).
"crupexp_test_data"

#' DNase-DNase test
#'
#' Results of the Wilcoxon tests on the DNase data of the pairs
#'@format A data frame with three variables: \code{gene_id} Gene ID,
#' \code{symbol38} Enhancer ID, \code{wilcoxtest_dhs_dhs} result of the test ).
"dhsdhs_test_data"

#' DHS expression test
#'
#' Results of the Wilcoxon test on the DNase data and RNA expression of the pairs
#'
#' @format A data frame with three variables: \code{gene_id} Gene ID,
#' \code{symbol38} Enhancer ID, \code{wilcoxtest_dhs_exp} result of the test ).
"dhsexp_test_data"
