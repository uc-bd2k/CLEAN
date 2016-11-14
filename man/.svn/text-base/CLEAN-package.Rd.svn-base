\name{CLEAN-package}
\alias{CLEAN-package}
\alias{CLEAN}
\docType{package}
\title{
Clustering Enrichment Analysis
}
\description{
Given an hierarchical gene clustering and a list of functional categories, this package performs functional enrichment analysis of all possible clusters and generates files to simultaneously display gene expression data, gene clustering, sample clustering, and functional annotation of gene clustering.
}
\details{
\tabular{ll}{
Package: \tab CLEAN\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2008-08-23\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}
}
\author{
Johannes Freudenberg, Xiangdong Liu, Mario Medvedovic

Maintainer: Mario Medvedovic <Mario.Medvedovic@uc.edu>
}
\references{

}
\keyword{ package }
\seealso{
 \code{\link[gimmR:gimmR-package]{gimmR}}
}
\examples{
data(gimmOut)
require(CLEAN.Rn)
res <- runCLEAN(gimmOut, species = "Rn")

generateTreeViewFiles(gimmOut, functionalCategories=getFunctionalCategories("geneRIFs", species = "Rn"))
#same as 
generateTreeViewFiles(gimmOut, functionalCategories="geneRIFs", species = "Rn")

#multiple category types
generateTreeViewFiles(gimmOut, functionalCategories=c("geneRIFs", "CpGislands", "GO", "KEGG"), species = "Rn")

trt <- sapply(colnames(gimmOut$clustData)[-(1:2)], function(str) strsplit(str, split = "_")[[1]][1])
#not run
#generateTreeViewFiles(gimmOut, cclust = NA, verbose = FALSE, functionalCategories=c("geneRIFs",
#	"CpGislands", "GO", "KEGG"), species = "Rn", callTreeView = TRUE, sampleDesc = trt)
generateTreeViewFiles(gimmOut, cclust = NA, verbose = FALSE, functionalCategories=c("geneRIFs",
	"CpGislands", "GO", "KEGG"), species = "Rn", callTreeView = FALSE, sampleDesc = trt)

#non-hierarchical clustering
d <- nonHierarchicalClustering(function(m, k, ...) kmeans(m, k, ...)$cluster, 
	gimmOut$clustData[,-(1:2)], k = 2:4, nstart = 10)
#not run
#generateTreeViewFiles(gimmOut, rclust = d, cclust = NA, verbose = FALSE, functionalCategories=c("geneRIFs",
#    "CpGislands", "GO", "KEGG"), species = "Rn", callTreeView = TRUE, sampleDesc = trt)
generateTreeViewFiles(gimmOut, rclust = d, cclust = NA, verbose = FALSE, functionalCategories=c("geneRIFs",
    "CpGislands", "GO", "KEGG"), species = "Rn", callTreeView = FALSE, sampleDesc = trt)

#geneList enrichment
geneList <- gimmOut$clustData[,1]
require(org.Rn.eg.db)
allGenes <- unique(keys(org.Rn.egSYMBOL)) #one should really use the list of 
                                          #genes represented on the microarray instead
res <- geneListEnrichment(geneList, allGenes, functionalCategories = "GO", 
	species = "Rn", sigFDR = 0.01, maxGenesInCategory = 10000)
genesInEnrichedCategories(res[,1], geneList, funcCategories = "GO", species = "Rn")

}

