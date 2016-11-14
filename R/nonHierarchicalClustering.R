`nonHierarchicalClustering` <-
function(clustering, data = NULL, k = NULL, nRepeats = 1, ...) {
#usage:
# data(gimmOut)
# d <- nonHierarchicalClustering(function(m, k, ...) kmeans(m, k, ...)$cluster, gimmOut$clustData[,-(1:2)], nstart = 10) # all possible cluster numbers (k)
# funcClustAnnot(gimmOut$clustData, d, NA, callTreeView=T)
# d <- nonHierarchicalClustering(function(m, k, ...) kmeans(m, k, ...)$cluster, gimmOut$clustData[,-(1:2)], 4, nstart = 10) # k = 4
# funcClustAnnot(gimmOut$clustData, d, NA, callTreeView=T)
# d2 <- nonHierarchicalClustering(function(m, k, ...) kmeans(m, k, ...)$cluster, t(gimmOut$clustData[,-(1:2)]), 9, nstart = 10) #sample clustering
# funcClustAnnot(gimmOut$clustData, d, d2, callTreeView=T)
	if(is.function(clustering)) {
		if(is.null(k)) {
			k <- 2:round(dim(data)[1]/2)
			k <- rep(k, each = nRepeats)
		}
		if(is.null(data)) {
			warning("data = NULL, no clustering generated.")
			clustering <- NULL
		} else {
			clustering <- t(sapply(k, function(i) clustering(data, i, ...)))
		}
	} else if(is.vector(clustering)) {
		clustering <- t(clustering)
	}
	n <- dim(clustering)[2]
	m <- dim(clustering)[1]
	d <- matrix(0, n, n) 
	for(i in 1:m) {
		cl <- split(1:n, clustering[i,])
		for(j in 1:length(cl)) {
			index <- cl[[j]]
			d[index, index] <- d[index, index] + 1
		}
	}
	#d <- as.dist(1-d/m)
	d <- as.dist(1-d/m + runif(length(d), 0, 0.05))
	#d <- hclust(d, method = "average")
	#d$labels <- 1:n
	d
}

