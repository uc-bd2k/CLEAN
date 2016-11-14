`RS` <-
function(sigvals, geneids, 
    database="GO", functionalCategories=NULL, species="Hs", 
	min.g=10, minGenesInCategory=NULL, 
	max.g=NA, maxGenesInCategory=NULL, 
	sig.cutoff=0.1, sigFDR=NULL,
	averageMultipleProbes=TRUE,
	allGenesInCategoriesAsBackground=FALSE,
	two.sided=FALSE, na.rm=TRUE,
	verbose=TRUE) {

	## sanity checks etc.

	if(!missing(database)) {
		warning("Argument `database` deprecated.  Please use `functionalCategories` instead.")
		if(missing(functionalCategories)) functionalCategories <- database
	}
	if (is.null(functionalCategories)) functionalCategories <- "GO"

	if(!(missing(min.g) & missing(max.g))) {
		warning("Argument `min.g`, `max.g` deprecated.  Please use `minGenesInCategory`, `maxGenesInCategory` instead.")
		if(missing(minGenesInCategory)) minGenesInCategory <- min.g
		if(missing(maxGenesInCategory)) maxGenesInCategory <- max.g
	}
	
	if(!missing(sig.cutoff)) {
		warning("Argument `sig.cutoff` deprecated.  Please use `sigFDR` instead.")
		if(missing(sigFDR)) sigFDR <- sig.cutoff
	}

	RandomSet(sigvals=sigvals, 
		geneids=geneids,
		functionalCategories=functionalCategories, 
		species=species,
		minGenesInCategory=minGenesInCategory,
		maxGenesInCategory=maxGenesInCategory,
		sigFDR=sigFDR,
		averageMultipleProbes=averageMultipleProbes,
		allGenesInCategoriesAsBackground=allGenesInCategoriesAsBackground,
		two.sided=two.sided,
		na.rm=na.rm,
		verbose=verbose)
}
