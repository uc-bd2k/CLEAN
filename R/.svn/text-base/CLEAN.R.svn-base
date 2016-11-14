`CLEAN` <-
function(tree, allGenes, CategoryID2GeneID, CategoryID2Desc, refGenes, maxSize=1000, minSize=10, maxNumOfClusters=NULL, sigFDR = 0.1,
		minGenesInCategory=10, maxGenesInCategory=1000, inBkg=TRUE, estimateNullDistribution = FALSE, na.rm=TRUE,
		saveDataObjects = FALSE, verbose = TRUE) {
	
	if(is.matrix(CategoryID2GeneID) | is.data.frame(CategoryID2GeneID)) {
		##remove primary datasets without "proper" distribution
		m0 <- apply(CategoryID2GeneID[,-1], 2, function(p) {
			pi0 <- 0.99
			tryCatch(pi0 <- qvalue(p)$pi0, error = function(e) {pi0 <- 1})
			#print(pi0)
			pi0})
		if(all(m0 >= 1)) warning("None of the reference sets appear to have significant p-values. Please check.")
		else if(any(m0 >= 1)) {
			warning(paste("One or more reference sets appear to not have significant p-values. ", sum(m0 >= 1), "data sets removed from analysis."))
			CategoryID2GeneID <- CategoryID2GeneID[,-c(which(m0 >= 1)+1)]
		}
		l <- primSetSigInClustering(tree = tree, allGenes = allGenes, primRefSetPvalues=CategoryID2GeneID,
	    	refGenes = refGenes, maxSize = maxSize, minSize = minSize, maxNumOfClusters=maxNumOfClusters,
			na.rm=na.rm, verbose=verbose, sigFDR = sigFDR, inBkg=inBkg)
		CLEANscore <- ifelse(allGenes %in% CategoryID2GeneID[,1], l$minFDR, 0)
		if(estimateNullDistribution) {
			CLEANscore_nullDistr <- primSetSigInClustering_nullDistr(allClusters = l$allClusters, 
				allGenes = allGenes[sample(length(allGenes),length(allGenes))], 
				primRefSetPvalues = CategoryID2GeneID, refGenes = refGenes,
				na.rm=na.rm, verbose=verbose, sigFDR = sigFDR, inBkg=inBkg)
		} else {
			CLEANscore_nullDistr <- NULL
		}

	} else {
		l <- geneSigInClustering(tree = tree, allGenes = allGenes, CategoryID2GeneID = CategoryID2GeneID,
			refGenes = refGenes, maxSize = maxSize, minSize = minSize, maxNumOfClusters=maxNumOfClusters, 
			minGenesInCategory = minGenesInCategory, maxGenesInCategory = maxGenesInCategory, 
			verbose=verbose, inBkg = inBkg, sigFDR = sigFDR)
		geneMatrixSig <- matrix(1, dim(l$geneMatrixSigs)[1], dim(l$geneMatrixSigs)[2])
		refCategories <- l$refCategories 
		for(i in 1:length(l$sigCategories)) {
			#clusterGenes <- refCategories$AllGeneIDs[l$allClusters[[i]]]
			clusterGenes <- allGenes[l$allClusters[[i]]]
			if(!is.null(l$sigCategories[[i]])) {
				if(dim(l$sigCategories[[i]])[1] > 0) {
					for(j in 1:(dim(l$sigCategories[[i]])[1])) {
						categoryGenes <- unique(refCategories$CategoriesWithGenes[[as.character(l$sigCategories[[i]][j,1])]])
						index <- intersect(l$allClusters[[i]][!is.na(match(clusterGenes, categoryGenes))], which(geneMatrixSig[, i] > l$sigCategories[[i]][j,3]))
						geneMatrixSig[index, i] <- l$sigCategories[[i]][j,3]
					}
				}
			}
		}
		CLEANscore <- -log10(apply(geneMatrixSig, 1, min))
		if(estimateNullDistribution) {
			CLEANscore_nullDistr <- geneSigInClustering_nullDistr(allClusters = l$allClusters, 
				allGenes = allGenes[sample(length(allGenes),length(allGenes))], 
				CategoryID2GeneID = CategoryID2GeneID, refCategories = l$refCategories,  
				minGenesInCategory = minGenesInCategory, maxGenesInCategory = maxGenesInCategory, 
				verbose = verbose, inBkg = inBkg, sigFDR = sigFDR)
		} else {
			CLEANscore_nullDistr <- NULL
		}
	}

	if(any(l$minFDR > 0)) cwCLEANscore <- pruneMinFDR(l) else cwCLEANscore <- l$sminFDR

	if(!saveDataObjects) {
		l$geneMatrixSigs <- NULL
		l$geneMatrixORs <- NULL
	}
	
	list(c(l, CategoryID2Desc = CategoryID2Desc, allGenes = list(allGenes), 
		CLEANscore = list(CLEANscore), cwCLEANscore = list(cwCLEANscore), 
		CLEANscore_nullDistr = list(CLEANscore_nullDistr)))
}
