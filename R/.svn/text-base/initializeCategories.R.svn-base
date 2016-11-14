`initializeCategories` <-
function(CategoryID2GeneID, AllGeneIDs=NULL, minGenesInCategory=10, maxGenesInCategory=1000, inBkg=TRUE, verbose=TRUE, filterEmptyCategories=TRUE){
	
	if(is.null(AllGeneIDs)) {
		AllGeneIDs <- unique(unlist(CategoryID2GeneID))
	}
	if(verbose) cat("Total Number of categories \n",length(CategoryID2GeneID),"\n")
	
	l <- sapply(CategoryID2GeneID, length)
	CategoryID2GeneID <- CategoryID2GeneID[l >= minGenesInCategory]
	
	#Getting frequencies of genes on the microarray in each category
	if(filterEmptyCategories) {
		GenesInCategories <- sapply(CategoryID2GeneID, function(x) length(intersect(AllGeneIDs, x)))
		#GenesInCategories <- sapply(CategoryID2GeneID, function(x) sum(!is.na(match(AllGeneIDs, x))))
	} else {
		GenesInCategories <- sapply(CategoryID2GeneID, length)
	}
	nonEmptyCategories <- GenesInCategories > 0
	GenesInCategories <- GenesInCategories[nonEmptyCategories]
	allGenesInCategories <- intersect(unique(unlist(CategoryID2GeneID[nonEmptyCategories])), AllGeneIDs)
	nGenesInCategories <- length(allGenesInCategories)
	if(verbose) cat("Number of genes on microarray represented in functional categories\n",nGenesInCategories,"\n")
	
	CategoriesWithGenes <- CategoryID2GeneID[names(GenesInCategories)]

	if(verbose) cat("Number of categories with genes from the list\n",length(GenesInCategories),"\n")

	minMaxCategories <- names(GenesInCategories[GenesInCategories >= minGenesInCategory & GenesInCategories <= maxGenesInCategory])
	nMinMaxCategories <- length(minMaxCategories)
	if(verbose) cat("Number of categories with at least", minGenesInCategory, "and not more than", maxGenesInCategory, "genes\n", nMinMaxCategories, "\n")
	#Getting all genes in such categories
	CategoriesWithGenes <- CategoriesWithGenes[minMaxCategories]
	allGenes <- intersect(unique(unlist(CategoriesWithGenes)), allGenesInCategories)

	if(inBkg) {
		AllGeneIDs <- allGenes
		if(verbose) cat("Only genes occuring in considered categories are used in the analysis\n")
	} 
	nAllGenes<-length(AllGeneIDs)
	
	if(verbose) cat("Number of genes in the background list of categories\n", nAllGenes,"\n")

	CategoriesWithGenes <- lapply(CategoriesWithGenes, function(l) intersect(AllGeneIDs, l))
	#CategoriesWithGenes <- lapply(CategoriesWithGenes, function(l) l[match(AllGeneIDs, l, 0)])

	list(GenesInCategories=GenesInCategories[names(CategoriesWithGenes)],
		CategoriesWithGenes=CategoriesWithGenes,
		allGenesInCategories=allGenesInCategories,
		nGenesInCategories=nGenesInCategories, 
		AllGeneIDs=AllGeneIDs)
}

