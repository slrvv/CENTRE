#' Download PrecomputedDataLight.db
#'
#' Downloads the PrecomputedDataLight.db database that is needed for the
#' computeGenericFeatures() function. This databased contains the values
#' of the Wilcoxon Rank tests and the CRUP correlations.
#'
#' @param method Method to be used for downloading files. Check download.file()
#' manuals to see the current available methods
#' @param all Boolean value. When set to true, the PrecomputedDataAllTests.db
#' database will be downloaded too. This database holds all of the p-values
#' which are combined in "combined_tests" column.
#' @return Integer code 0 if download was succesful and the file was saved in
#' the correct directory. In the case of failure errors will be raised.
#' @examples
#' # methods accepts "internal", "wininet" (Windows only) ,"wget", "curl" etc.
#' # Assuming the user has wget available
#' downloadPrecomputedData(method="wget")
#' @export
#' @import utils

downloadPrecomputedData <- function(method, all = FALSE) {
  start_time <- Sys.time()
  #Download PrecomputedDataLight.db
  downloader("PrecomputedDataLight.db", method)
  #Download Annotation.db
  downloader("Annotation.db", method)
  if (all == TRUE) {
    downloader("PrecomputedAllTests.db", method)
  }
  cat(paste0("time: ", format(Sys.time() - start_time), "\n"))
  return(0)
}
