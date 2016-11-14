`primSetSigInClustering_nullDistr` <-
function(allClusters, allGenes, primRefSetPvalues, refGenes = NULL, verbose=TRUE, sigFDR = 0.1, na.rm=TRUE, inBkg=TRUE){

	## RandomSet function (Newton et al.)
	myAllez <- function(pValues, categoryList.index) {
    	#categoryList <- unique(categoryList)
	    #index <- na.omit(match(categoryList, geneList))
		scores <- -log(pValues, 10)
		index <- which(!is.na(scores[categoryList.index]))
		m <- length(index)
		G <- ifelse(na.rm, length(na.omit(scores)), length(scores))
	    if(m > 0) {
			xBar <- mean(scores[categoryList.index[index]])
			mu <- mean(scores, na.rm=na.rm)
			sigma <- sqrt((sum(scores^2, na.rm=na.rm)/G-mu^2) * (G - m) / (G - 1) / m)
			if(sigma==0) sigma <- 1 ## mu and xBar are the same
			#list(xBar=xBar, mu=mu, sigma=sigma, z=(xBar-mu)/sigma)
			c(zScore=(xBar-mu)/sigma, n.genes=m, n.allGenes=G)
		} else {
			#list(NA)
			c(zScore=NA, n.genes=0, n.allGenes=G)
		}
	}
	
	#allClusters<-getAllClusters(tree=tree, maxSize=maxSize, minSize=minSize, clusterList=NULL, k=maxNumOfClusters)
	if(verbose) cat("Number of Clusters\n",length(allClusters),"\n")

	geneMatrixSigs<-matrix(1,length(allGenes),length(allClusters))
	#sigCategories<-NULL

	if(inBkg) {
		if(is.null(refGenes)) {
			refGenes <- allGenes
		}
		refGenes <- intersect(refGenes, primRefSetPvalues[,1])
		pValues <- primRefSetPvalues[match(refGenes, primRefSetPvalues[,1]), -1]
	} else {
		refGenes <- primRefSetPvalues[,1]
		pValues <- primRefSetPvalues[, -1]
	}
	pValues <- as.matrix(pValues)
	refIndex <- match(allGenes, refGenes)

	for(i in 1:length(allClusters)){
		if(verbose) cat("Cluster ", i, " of ", length(allClusters), "\n")
		currentCluster <- refIndex[allClusters[[i]]]
		sigPrimDatasets <- apply(pValues, 2, myAllez, categoryList.index=currentCluster)
		sigPrimDatasets <- t(sigPrimDatasets)
		sigPrimDatasets[, 1] <- pnorm(sigPrimDatasets[, 1], lower.tail=FALSE)
		sigPrimDatasets[, 1] <- ifelse(sigPrimDatasets[, 1] < .Machine$double.xmin, .Machine$double.xmin, sigPrimDatasets[, 1])
		sigPrimDatasets <- cbind(sigPrimDatasets, RS.FDR=p.adjust(sigPrimDatasets[, 1], method="fdr"))
		sigPrimDatasets <- sigPrimDatasets[order(sigPrimDatasets[,"RS.FDR"]),]
		#if(verbose) print((sigPrimDatasets[order(sigPrimDatasets[, 1]),])[1:10,])
		#if(verbose) print(allClusters[[i]])
		geneMatrixSigs[as.integer(allClusters[[i]]), i]<-sigPrimDatasets[1, "RS.FDR"]
		
		#index <- sigPrimDatasets[, "RS.FDR"] < sigFDR
		#if(sum(index, na.rm=TRUE) > 0) {
		#	sigPrimDatasets <- matrix(sigPrimDatasets[index, ], ncol=4)
		#	sigPrimDatasets <- data.frame(ID=colnames(pValues)[index], FisherPValue=sigPrimDatasets[, 1], 
		#		FisherFDR=sigPrimDatasets[, 4], nGenesInCategory=sigPrimDatasets[, 2], 
		#		nAllGenesInCategory=sigPrimDatasets[, 3], logOR=NA)
		#} else {
		#	sigPrimDatasets <- data.frame(ID=character(0), FisherPValue=numeric(0),
		#		FisherFDR=numeric(0), nGenesInCategory=numeric(0), nAllGenesInCategory=numeric(0), logOR=numeric(0))
		#}
		#sigCategories<-c(sigCategories,list(sigPrimDatasets))
	}
	-log10(apply(geneMatrixSigs,1,min))
}

