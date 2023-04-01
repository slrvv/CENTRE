# Installation of dependencies required by CENTRE


pkgInstall <- function(pkg) { 

  if(! isTRUE(pkg %in% .packages(all.available=T))) { 
    
    x <- tryCatch({eval(parse(text = sprintf("install.packages(\"%s\", dependencies = T)", pkg)))})

    if(is.null(x)) {

	
  	if (!'BiocManager' %in% installed.packages()){
		if (!requireNamespace("BiocManager", quietly = TRUE))
    		install.packages("BiocManager")
		BiocManager::install()
	}
	x <- tryCatch({eval(parse(text = sprintf("BiocManager::install(\"%s\")", pkg)))})
    }

    if(is.null(x)) {
	cat(paste0("Unable to install package: ",pkg,".\n"));
  	q();
    }else{
	    eval(parse(text = sprintf("require(\"%s\")", pkg)))
    }
  }
}


pkgTest <- function(pkg){
  if (!pkg %in% installed.packages()) {
	cat(paste0("package ",pkg, " not found. Install now."))
	return(F)
  }else{
	return(T)	
  }
} 


pkgLoad <- function(pkg) {
  if(pkgTest(pkg) == F){
 	pkgInstall(pkg)
  }else{
	suppressMessages(library(pkg, character.only = TRUE))

  }
}

pkg <- c("GenomicRanges", "IRanges", "RSQLite", "metapod", "stats", "xgboost", "devtools")
for (i in seq_along(pkg)){
	pkgLoad(pkg[i])
}

devtools::install_github("akbariomgba/crupR")
