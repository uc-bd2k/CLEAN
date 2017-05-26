`significantCategories` <-
function(genes, refCategories, lean=TRUE, verbose = TRUE){

	#Selecting categories with at least minGenesInCategory and no more than maxGenesInCategory genes
#	minMaxCategories<-names(refCategories$GenesInCategories[refCategories$GenesInCategories>=minGenesInCategory & 
#														refCategories$GenesInCategories<=maxGenesInCategory])
	minMaxCategories <- names(refCategories$GenesInCategories)

	nMinMaxCategories<-length(minMaxCategories)
#	if(verbose) cat("Number of categories with more than ",minGenesInCategory," genes\n",nMinMaxCategories,"\n")

	#Getting all genes in such categories
#	allGenes<-refCategories$CategoriesWithGenes[minMaxCategories]
#	CategoriesWithGenesMinMax<-allGenes
#	allGenes<-intersect(unique(unlist(CategoriesWithGenesMinMax)),refCategories$allGenesInCategories)

	allGenes <- refCategories$AllGeneIDs

	genesInUsedCategories<-intersect(genes,allGenes)

#	nGenes<-length(genes)
#	if(inBkg) nGenes<-length(genesInUsedCategories)
	nGenes<-length(genesInUsedCategories)
#	if(inBkg) nAllGenes<-length(allGenes)
#	else nAllGenes<-length(refCategories$AllGeneIDs)
	nAllGenes<-length(allGenes)
	nComplementGenes<-nAllGenes-nGenes
	#if(inBkg & verbose) cat("Only genes in occuring in considered categories are used in the analysis\n")
#	if(verbose) cat("Number of genes in the background list of categories\n",nAllGenes,"\n")
	if(verbose) cat("Number of genes\n",nGenes,"\n")

	#Testing categories with at least one gene
	nGenesInCategoryAll<-unlist(lapply(refCategories$CategoriesWithGenes,function(x) length(intersect(x,genes))))
	nonEmptyCategories<-nGenesInCategoryAll>0
	if(verbose) cat("Number of categories with significant genes\n",sum(nonEmptyCategories),"\n")
	nGenesInCategoryAll<-nGenesInCategoryAll[nonEmptyCategories]
	CategoriesWithGenes<-names(nGenesInCategoryAll)
	matchWithMinMax<-match(CategoriesWithGenes,minMaxCategories)
	CategoriesWithGenes<-CategoriesWithGenes[!is.na(matchWithMinMax)]
	nGenesInCategoryAll<-nGenesInCategoryAll[!is.na(matchWithMinMax)]
	if(verbose) cat("Number of categories with significant genes tested \n",sum(!is.na(matchWithMinMax)),"\n")


	#nGenesInCategories<-NULL
	#nGenesInCategoriesAll<-NULL
	#FisherPValues<-NULL
	#logOR<-NULL
	
#	if(inBkg) m <- unlist(lapply(refCategories$CategoriesWithGenes[CategoriesWithGenes], 
#		function(l) length(intersect(allGenes, l))))
#	else m <- unlist(lapply(refCategories$CategoriesWithGenes[CategoriesWithGenes], 
#		function(l) length(intersect(refCategories$AllGeneIDs, l))))
	if(sum(!is.na(matchWithMinMax)) > 0) {
		m <- unlist(lapply(refCategories$CategoriesWithGenes[CategoriesWithGenes], length))
		n <- nAllGenes - m
		k <- nGenes
		#a <- unlist(lapply(refCategories$CategoriesWithGenes[CategoriesWithGenes], function(l) length(intersect(genes, l))))
		a <- nGenesInCategoryAll
		FisherPValues <- phyper(a - 1, m, n, k, lower.tail = FALSE)
		nGenesInCategories <- a
		nGenesInCategoriesAll <- m
		logOR <- log((a+1)*(n - k + a + 1)/((k - a + 1)*(m - a + 1)))
		#no fudge factor if sum(twoBytwo)==4 --> ever an issue?
	} else {
		FisherPValues <- NULL
		nGenesInCategories <- NULL
		nGenesInCategoriesAll <- NULL
		logOR <- NULL
	}
#	for(x in CategoriesWithGenes){
#		currentCategoryGenes<-unique(unlist(refCategories$CategoriesWithGenes[x]))
#		nGenesInCategory<-length(intersect(genes,currentCategoryGenes))
#		nGenesNOTinCategory<-nGenes-nGenesInCategory
#		if(inBkg) nComplementGenesInCategory<-length(intersect(setdiff(allGenes,genes),currentCategoryGenes))
#		else nComplementGenesInCategory<-length(intersect(setdiff(refCategories$AllGeneIDs,genes),currentCategoryGenes))
#		nComplementGenesNOTinCategory<-nComplementGenes-nComplementGenesInCategory
#		twoBytwo<-matrix(c(nGenesInCategory,nGenesNOTinCategory,nComplementGenesInCategory,nComplementGenesNOTinCategory), byrow=TRUE,ncol=2)
#		FisherPValues<-c(FisherPValues,fisher.test(twoBytwo, alternative = "greater")$p.value)
#		if(sum(twoBytwo)==4) logOR<-c(logOR,log((twoBytwo[1,1]*twoBytwo[2,2])/(twoBytwo[1,2]*twoBytwo[2,1])))
#		else logOR<-c(logOR,log(((twoBytwo[1,1]+1)*(twoBytwo[2,2]+1))/((twoBytwo[1,2]+1)*(twoBytwo[2,1]+1))))
#		nGenesInCategories<-c(nGenesInCategories,nGenesInCategory)
#		nGenesInCategoriesAll<-c(nGenesInCategoriesAll,nComplementGenesInCategory+nGenesInCategory)
#	}
#browser()
	nPs<-length(FisherPValues)
	if(nPs>0){
		pOrder<-order(FisherPValues)
		pOriginal<-order(pOrder)
		FisherFDR<-FisherPValues[pOrder]*nMinMaxCategories/(1:nPs)
		for(i in 1:nPs) FisherFDR[i]<-min(FisherFDR[i:nPs],1)
		FisherFDR<-FisherFDR[pOriginal]
		FishBonf=FisherPValues*nMinMaxCategories
		FisherBonf=ifelse(FishBonf<1,FishBonf,1)
		oPvalues<-order(FisherPValues)
		CategoriesWithAllGenes<-refCategories$CategoriesWithGenes[CategoriesWithGenes]
		if(lean) list(categories=data.frame(list(ID=CategoriesWithGenes[oPvalues],list(FisherPValue=unlist(FisherPValues[oPvalues]),FisherFDR=FisherFDR[oPvalues],FisherBonf=FisherBonf[oPvalues],
			nGenesInCategory=nGenesInCategories[oPvalues],nAllGenesInCategory=nGenesInCategoriesAll[oPvalues],logOR=logOR[oPvalues])),stringsAsFactors=F))
		else list(categories=data.frame(list(ID=CategoriesWithGenes[oPvalues],list(FisherPValue=unlist(FisherPValues[oPvalues]),FisherFDR=FisherFDR[oPvalues],FisherBonf=FisherBonf[oPvalues],
			nGenesInCategory=nGenesInCategories[oPvalues],nAllGenesInCategory=nGenesInCategoriesAll[oPvalues],logOR=logOR[oPvalues])),stringsAsFactors=F),
			geneNumbers=list(nGenes=nGenes,nAllGenes=nAllGenes),genes=genes,allGenes=allGenes,CategoriesWithAllGenes=CategoriesWithAllGenes)
	}
	else NULL
}

