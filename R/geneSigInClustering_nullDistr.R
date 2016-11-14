`geneSigInClustering_nullDistr` <-
function(allClusters,allGenes,CategoryID2GeneID,refCategories = NULL,minGenesInCategory=10,maxGenesInCategory=1000,verbose=TRUE,inBkg=FALSE, sigFDR = 0.1){
	
	if(verbose) cat("Number of Clusters\n",length(allClusters),"\n")
	if(is.null(refCategories)) {
		refGenes<-unique(allGenes)
		refCategories<-initializeCategories(CategoryID2GeneID, AllGeneIDs=refGenes,verbose=FALSE, 
			minGenesInCategory=minGenesInCategory, maxGenesInCategory=maxGenesInCategory, inBkg = inBkg)
	}	
	geneMatrixSigs<-matrix(1,length(allGenes),length(allClusters))
	#geneMatrixORs<-matrix(0,length(allGenes),length(allClusters))
	sigCategories<-NULL
	if (verbose) cat("\n")
	for(i in 1:length(allClusters)){
		if(verbose) cat("Cluster ", i, " of ", length(allClusters), "\n")
		sigcategories<-significantCategories(unique(allGenes[allClusters[[i]]]), refCategories=refCategories, verbose = FALSE)
		#if(verbose) print(sigcategories$categories[1:10,])
		#if(verbose) print(allClusters[[i]])
		if(length(sigcategories$categories)>0) {
			geneMatrixSigs[as.integer(allClusters[[i]]), i]<-sigcategories$categories[1,"FisherFDR"]
			#geneMatrixORs[as.integer(allClusters[[i]]), i]<-max(sigcategories$categories[1,"logOR"])
		}
		sigCategories<-c(sigCategories,list(sigcategories$categories[sigcategories$categories[,"FisherFDR"]<sigFDR,]))
	}
	#minFDR<- -log10(apply(geneMatrixSigs,1,min))
	#maxOR<-apply(geneMatrixORs,1,max)
	
	geneMatrixSig <- matrix(1, dim(geneMatrixSigs)[1], dim(geneMatrixSigs)[2])
	for(i in 1:length(sigCategories)) {
		clusterGenes <- refCategories$AllGeneIDs[allClusters[[i]]]
		if(!is.null(sigCategories[[i]])) {
			if(dim(sigCategories[[i]])[1] > 0) {
				for(j in 1:(dim(sigCategories[[i]])[1])) {
					categoryGenes <- unique(refCategories$CategoriesWithGenes[[as.character(sigCategories[[i]][j,1])]])
					index <- intersect(allClusters[[i]][!is.na(match(clusterGenes, categoryGenes))], which(geneMatrixSig[, i] > sigCategories[[i]][j,3]))
					geneMatrixSig[index, i] <- sigCategories[[i]][j,3]
				}
			}
		}
	}
	#list(minFDR_nullDistr=minFDR,maxOR_nullDistr=maxOR)#,geneMatrixSigs=geneMatrixSigs,geneMatrixORs=geneMatrixORs,allClusters=allClusters,sigCategories=sigCategories)
	-log10(apply(geneMatrixSig, 1, min))
}

