`pruneMinFDR` <-
function(clustCategories) {
	clusters <- clustCategories$allClusters
	clustVals <- unlist(lapply(apply(clustCategories$geneMatrixSigs, 2, function(x) unique(x[x < 1])), 
		function(x) if(length(x) < 1) 1 else x))
	clustSize <- unlist(lapply(clusters, length))
	geneMatrixSig <- clustCategories$geneMatrixSigs
	toDo <- (1:length(clusters))[order(clustVals)]
	while(length(toDo) > 0 & clustVals[toDo[1]] < 1) {
		cluster <- toDo[1]
		index <- toDo[-1]
		if(length(index) > 0) {
			index <- index[sapply(index, function(j) length(intersect(clusters[[j]], clusters[[cluster]])) > 0)]
			toDo <- setdiff(toDo, index)
			geneMatrixSig[,index] <- 1
		}
		toDo <- toDo[-1]
	}
	-log10(apply(geneMatrixSig,1, min))
}

