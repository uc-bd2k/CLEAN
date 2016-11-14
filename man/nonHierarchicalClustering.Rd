\name{nonHierarchicalClustering}
\alias{nonHierarchicalClustering}
\title{Convert a non-hierarchical clustering into hierarchical tree-structure }
\description{
coming soon.
}
\usage{
nonHierarchicalClustering(clustering, data = NULL, k = NULL, nRepeats = 1, ...)
}
\arguments{
  \item{clustering}{A vector assigning genes or samples to clusters, or a function}
  \item{data}{The data to be clustered if clustering is a function.}
  \item{k}{The number of clusters to be generated}
  \item{nRepeats}{The number of runs to be averaged.}
  \item{\dots}{Additional parameters to be passed to clustering}
}
\details{
coming soon.
}
\value{
a dist object.
}
\references{coming soon. }
\author{Johannes Freudenberg }
\note{coming soon.
}
\seealso{ \code{\link{funcClustAnnot}} }
\examples{
data(gimmOut)
require(CLEAN.Rn)
d <- nonHierarchicalClustering(function(m, k, ...) kmeans(m, k, ...)$cluster,
    gimmOut$clustData[,-(1:2)], k = 2:4, nstart = 10)
trt <- sapply(colnames(gimmOut$clustData)[-(1:2)], function(str) strsplit(str, split = "_")[[1]][1])
#not run
#generateTreeViewFiles(gimmOut, rclust = d, cclust = NA, verbose = FALSE, functionalCategories=c("geneRIFs",
#	"CpGislands", "GO", "KEGG"), species = "Rn", callTreeView = TRUE, sampleDesc = trt)
generateTreeViewFiles(gimmOut, rclust = d, cclust = NA, verbose = FALSE, functionalCategories=c("geneRIFs",
	"CpGislands", "GO", "KEGG"), species = "Rn", callTreeView = FALSE, sampleDesc = trt, maxNumOfClusters=2)

# all possible cluster numbers (k)
d <- nonHierarchicalClustering(function(m, k, ...) kmeans(m, k, ...)$cluster, 
	gimmOut$clustData[,-(1:2)], nstart = 10)
# k = 4
d <- nonHierarchicalClustering(function(m, k, ...) kmeans(m, k, ...)$cluster, 
	gimmOut$clustData[,-(1:2)], 4, nstart = 10)
#sample clustering
d2 <- nonHierarchicalClustering(function(m, k, ...) kmeans(m, k, ...)$cluster, t(gimmOut$clustData[,-(1:2)]), 9, nstart = 10) #sample clustering

##not run
#generateTreeViewFiles(gimmOut$clustData, d, d2, callTreeView=TRUE, species = "Rn")
generateTreeViewFiles(gimmOut$clustData, d, d2, callTreeView=FALSE, species = "Rn")
}
\keyword{cluster }
