`geneSigInList` <- 
function(geneList,allGenes,CategoryID2GeneID,refGenes = NULL, 
	minGenesInCategory=10,maxGenesInCategory=1000,verbose=TRUE, inBkg=FALSE,
	sigFDR = 0.1) {
     
	if(is.null(refGenes)) refGenes <- unique(allGenes)
	refCategories<-initializeCategories(CategoryID2GeneID, AllGeneIDs=refGenes,verbose=verbose, 
		minGenesInCategory = minGenesInCategory, maxGenesInCategory = maxGenesInCategory, inBkg = inBkg)
	sigCategories<-significantCategories(unique(geneList), refCategories=refCategories, verbose = verbose)
	if(verbose) print(sigCategories$categories[1:10,])
	list(sigCategories$categories[sigCategories$categories[,"FisherFDR"]<sigFDR,])
}


