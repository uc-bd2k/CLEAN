`getFunctionalCategories` <-
function(functionalCategories, species=c("Hs", "Rn", "Mm")) {
	if(missing(species)|is.null(species)) {
		species <- "Hs"
		warning("Argument `species` should be one of `Hs`, `Mm`, `Rn`). Using the default `Hs`.")
	}
	species <- tryCatch(match.arg(species), 
		error = function(e) {warning("Argument `species` should be one of `Hs`, `Mm`, `Rn`). Using the default `Hs`.");"Hs"})
	require(AnnotationDbi)
    fCategories <- list()
    require(paste("CLEAN", species, sep = "."), character.only = TRUE)
    for(i in 1:length(functionalCategories)) {
           e1 <- new.env()
           data(list = paste(functionalCategories[i], species, sep = "."), envir = e1)
           if(length(ls(e1)) > 0) {
               fCategories <- c(fCategories, list(get(paste(functionalCategories[i], species, sep = "."), envir = e1)))
               names(fCategories)[length(fCategories)] <- functionalCategories[i]
        }
    }
    fCategories
}

