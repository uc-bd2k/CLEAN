`getAllClusters` <-
function(tree, maxSize=1000, minSize=10, clusterList=NULL, k=NULL) {
#recursive version of this function
#deprecated because unstable with some clustering algorithms
#	pos<-c(F,T)
#	for(side in pos) {
#		subTree<-getBinarySubTree(tree,side)
#		genes<-subTree$labels
#		if(length(genes) >= minSize && length(genes)<=maxSize) {
#			clusterList[[length(clusterList)+1]]<-genes
#		}
#		if(length(genes) >= minSize) {
#			clusterList<-getAllClusters(subTree,maxSize,minSize,clusterList)
#		}
#	}
#	clusterList
	if(missing(k) | is.null(k)) k <- Inf
	else { 
		if(k < 0) k <- 0
		if(nrow(tree$merge) < maxSize) k <- k + 1
	}
	treeList <- list(list(node = tree, splitMore = TRUE, isLeaf = FALSE))
	if(is.null(tree$labels)) warning("tree does not have labels, function may not generate any meaningful clusters.")
	if(length(tree$labels) <= maxSize & length(tree$labels) >= minSize) treeList[[1]]$isLeaf = TRUE
	if(minSize < 2) {
		warning(paste("minimum cluster size", minSize, "too small, set to 2."))
		minSize <- 2 
	}
	f <- function(treeListElement, maxSize, minSize) {
		treeList <- as.list(NULL)
		if(treeListElement$splitMore) {
			pos<-c(FALSE, TRUE)
			subTrees <- list(getBinarySubTree(treeListElement$node, FALSE),
				getBinarySubTree(treeListElement$node, TRUE))
			if(max(subTrees[[1]]$height) < max(subTrees[[2]]$height))
				subTrees <- subTrees[2:1]
			for(side in pos) {
				newTree <- NULL
				#subTree<-getBinarySubTree(treeListElement$node,side)
				subTree <- subTrees[[side+1]]
				subTree$call <- "getBinarySubTree(treeListElement$node,side)"
				newTree$node <- subTree
				genes<-subTree$labels
				if(length(genes) >= minSize) {
					newTree$splitMore <- TRUE
					if(length(genes)<=maxSize) {
						newTree$isLeaf <- TRUE
					} else {
						newTree$isLeaf <- FALSE
					}
					treeList <- c(treeList,list(newTree))
				} 
			}
		}
		treeList
	}
	allLeaves <- FALSE
	while(!allLeaves) {
		allLeaves <- TRUE
		trees <- NULL
		toBeRemoved <- NULL
		for(i in 1:length(treeList)) {
			newTrees <- NULL
			if(treeList[[i]]$splitMore) {
				newTrees <- f(treeList[[i]], maxSize, minSize)
				if(!treeList[[i]]$isLeaf) {
					toBeRemoved <- c(toBeRemoved, i)
				} else {
					treeList[[i]]$splitMore <- FALSE
				}
			}
			if(length(newTrees) > 0) {
				allLeaves <- FALSE
				trees <- c(trees, newTrees)
			}
		}
		if(!is.null(toBeRemoved)) treeList <- treeList[-toBeRemoved]
		treeList <- c(treeList, trees)
		if(length(treeList) > k) allLeaves <- TRUE
	}
	res <- lapply(treeList, function(subTree) subTree$node$labels)
	res[1:min(k, length(res))]
}

