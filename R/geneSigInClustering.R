`geneSigInClustering` <-
function(tree,allGenes,CategoryID2GeneID,refGenes = NULL,maxSize=1000,minSize=10,maxNumOfClusters=NULL, minGenesInCategory=10,maxGenesInCategory=1000,minSigsInCategory=2,verbose=TRUE,inBkg=TRUE, sigFDR = 0.1){
	
	allClusters<-getAllClusters(tree=tree,maxSize=maxSize,minSize=minSize,clusterList=NULL, k=maxNumOfClusters)
	if(verbose) cat("Number of Clusters\n",length(allClusters),"\n")
	if(is.null(refGenes)) refGenes<-unique(allGenes)
	refCategories<-initializeCategories(CategoryID2GeneID, AllGeneIDs=refGenes, minGenesInCategory = minGenesInCategory, 
		maxGenesInCategory = maxGenesInCategory, inBkg = inBkg, verbose=verbose)

	geneMatrixSigs<-matrix(1,length(allGenes),length(allClusters))
	geneMatrixORs<-matrix(0,length(allGenes),length(allClusters))
	sigCategories<-NULL
	
	for(i in 1:length(allClusters)){
		if(verbose) cat("\n\nCluster ", i, " of ", length(allClusters), "\n\n")
		sigcategories<-significantCategories(unique(allGenes[allClusters[[i]]]),refCategories=refCategories, verbose = verbose)
		if(verbose) print(sigcategories$categories[1:10,])
		if(verbose) print(allClusters[[i]])
		if(length(sigcategories$categories)>0) {
			geneMatrixSigs[as.integer(allClusters[[i]]), i]<-sigcategories$categories[1,"FisherFDR"]
			geneMatrixORs[as.integer(allClusters[[i]]), i]<-max(sigcategories$categories[1,"logOR"])
		}
		sigCategories<-c(sigCategories,list(sigcategories$categories[sigcategories$categories[,"FisherFDR"]<sigFDR,]))
	}
	minFDR<- -log10(apply(geneMatrixSigs,1,min))
	maxOR<-apply(geneMatrixORs,1,max)
	list(minFDR=minFDR,maxOR=maxOR,geneMatrixSigs=geneMatrixSigs,geneMatrixORs=geneMatrixORs,allClusters=allClusters,sigCategories=sigCategories,refCategories=refCategories)
}

