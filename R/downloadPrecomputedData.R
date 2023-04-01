#' Download PrecomputedData.db
#' 
#' Downloads the PrecomputedData.db database that is needed for the 
#' computeGenericFeatures() function. This databased contains the values
#' of the Wilcoxon Rank tests and the CRUP correlations.
#' 
#' @param method Method to be used for downloading files. Check download.file()
#' manuals to see the current available methods 
#' @return Integer code 0 if download was succesful and the file was saved in the
#' correct directory. In the case of failure errors will be raised.
#' @examples
#' # methods accepts "internal", "wininet" (Windows only) ,"wget", "curl" etc.
#' # Assuming the user has wget available
#' downloadPrecomputedData(method="wget")
#' @export
#' @import utils

downloadPrecomputedData <- function(method) {
	start_time <- Sys.time()
	url <- "http://owww.molgen.mpg.de/~CENTRE_data/PrecomputedData.db"
	cat("Downloading PrecomputedData.db\n")
	exit <- download.file(url, 
				destfile=system.file("extdata",
						"PrecomputedData.db",
						package = "CENTRE"),
						method = method)         
	if (exit != 0 ) {
		stop("Download of PrecomputedData.db failed. Non-zero exit status.")
        }

	f <- system.file("extdata","PrecomputedData.db", package = "CENTRE") 
	if (!file.exists(f)){
		stop("Download of of PrecomputedData.db failed or file was saved in the wrong directory.")
	}
	##Download sysdata.rda
	url_sys <- "http://owww.molgen.mpg.de/~CENTRE_data/Annotation.db"
	cat("Downloading sysdata.rda\n")
	exit_sys <- download.file(url_sys,
				destfile=system.file("extdata",
				"Annotation.db",
				package="CENTRE"),
				method = method)
	if (exit_sys != 0 ) {
                 stop("Download Annotation.db failed. Non-zero exit status.")
         }

         fs <- system.file("extdata","PrecomputedData.db", package = "CENTRE")
         if (!file.exists(fs)){
                 stop("Download Annotation.db failed or file was saved in the wrong directory.")
         }

	cat(paste0('time: ', format(Sys.time() - start_time), "\n"))
	return(0)
}
