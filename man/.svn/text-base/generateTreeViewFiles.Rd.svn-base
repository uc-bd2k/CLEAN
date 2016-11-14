\name{funcClustAnnot}
\alias{funcClustAnnot}
\alias{geneListEnrichment}
\alias{generateTreeViewFiles}
\alias{runCLEAN}
\title{Wrapper function to perform functional Clustering Enrichment Analysis (CLEAN)}
\description{
This function takes a gene expression data set, hierarchical clusterings of genes and samples, and a list of gene sets representing functional categories.  It performs hierarchical clustering (if not  provided) and performs the Clustering Enrichment Analysis.  Finally, it generates files to display data, clustering, and functional annotation using tools such as the Java-based fTreeView.
}
\usage{
funcClustAnnot(data = NULL, rclust = NULL, cclust = NULL, funcCategories = NULL, fClustAnnotations = NULL, bkgList = NULL, file = "cluster", maxSize = 1000, minSize = 10, minGenesInCategory = 10, maxGenesInCategory = 1000, minSigsInCategory = 2, inBkg = TRUE, sigFDR = 0.1, addCLEANscore2cdt = TRUE, species = NULL, estimateNullDistribution = FALSE, atr.distance = "euc", gtr.distance = "euc", dec = ".", digits = 5, verbose = TRUE, saveDataObjects = FALSE, sampleDesc = NULL, callTreeView = FALSE)
generateTreeViewFiles(data = NULL, rclust = NULL, cclust = NULL, functionalCategories = NULL, fClustAnnotations = NULL, bkgList = NULL, file = "cluster", maxSize = 1000, minSize = 10, maxNumOfClusters=NULL, minGenesInCategory = 10, maxGenesInCategory = 1000, inBkg = TRUE, sigFDR = 0.1, addCLEANscore2cdt = TRUE, species = NULL, estimateNullDistribution = FALSE, dec = ".", verbose = TRUE, saveDataObjects = FALSE, sampleDesc = NULL, callTreeView = FALSE)
geneListEnrichment(geneList, allGenes, functionalCategories = NULL, species = NULL, minGenesInCategory=10, maxGenesInCategory=1000, inBkg=TRUE, sigFDR = 0.1, verbose=TRUE)
runCLEAN(data = NULL, rclust = NULL, functionalCategories = NULL, bkgList = NULL, file = "cluster", maxSize=1000, minSize=10, maxNumOfClusters=NULL, minGenesInCategory=10, maxGenesInCategory=1000, inBkg=TRUE, sigFDR = 0.1, species = NULL, estimateNullDistribution = FALSE, verbose=TRUE, saveDataObjects = FALSE)
}
\arguments{
  \item{data}{ Gene expression data set.}
  \item{geneList}{A gene list to be tested for enrichment of functional categories.}
  \item{rclust}{Gene clustering. Can be file name, distance matrix, hclust object etc.  If hclust, rclust$labels must be a numeric vector indexing the rows of the data parameter.}
  \item{cclust}{Sample clustering.}
  \item{functionalCategories}{A collection of gene sets, representing functional categories.}
  \item{funcCategories}{A collection of gene sets, representing functional categories.}
  \item{fClustAnnotations}{A list of Clustering Enrichment Analysis results.  If provided, CLEAN is not performed which considerably reduces computing time.}
  \item{allGenes}{A list of genes to be used as background for gene list enrichment analysis.}
  \item{bkgList}{A list of genes to be used as background for CLEAN.}
  \item{file}{Name of the file(s) (without extensions) to be generated.}
  \item{maxSize}{Maximum cluster size to be considered by CLEAN.}
  \item{minSize}{Minimum cluster size to be considered by CLEAN.}
  \item{maxNumOfClusters}{Maximum number of clusters to be considered by CLEAN.  If specified, the top `maxNumOfClusters` from the tree cutting function will be chosen}
  \item{minGenesInCategory}{Minimum number of all genes in the expression data set overlapping with a given functional category.}
  \item{maxGenesInCategory}{Maximum number of all genes in the expression data set overlapping with a given functional category.}
  \item{minSigsInCategory}{Minimum number of genes in a cluster overlapping with a functional category.}
  \item{inBkg}{ -- to be completed -- } 
  \item{sigFDR}{FDR cutoff for a functional category to be overrepresented in a cluster.}
  \item{addCLEANscore2cdt}{If TRUE one or more columns are added to the cdt file indicating the highest significance level observed for a cluster by CLEAN.}
  \item{species}{Two letter description of the species to be used to generate gene ontology categories (e.g. "Hs" for human, "Mm" for mouse). This parameter is used when funcCategories = "GO"}
  \item{estimateNullDistribution}{If TRUE an empirical estimate of the Null-distribution of the CLEAN score is computed.}
  \item{atr.distance}{The distance metric used by the sample clustering -- deprecated.}
  \item{gtr.distance}{The distance metric used by the gene clustering -- deprecated.}
  \item{dec}{The decimal character to be used when writing the treeview files.}
  \item{digits}{Number of significant digits to be used when writing the treeview files -- deprecated. }
  \item{verbose}{If TRUE lots of additional output will be generated.}
  \item{saveDataObjects}{If TRUE, data, clusterings, and CLEAN results will be saved to disc as with filename <file>.RData.}
  \item{sampleDesc}{A vector or factor describing the samples.  Samples with same description will be grouped together if cclust = NA.}
  \item{callTreeView}{If TRUE the Java application fTreeView is called to display the CLEAN results.}
}
\details{
coming soon.
}
\value{
  A list:
  \item{data}{A dataframe containing the gene expression data.}
  \item{rclust, cclust}{Hierarchical clusterings (hclust objects) of genes and samples, respectively.}
  \item{fClustAnnotations}{CLEAN results.}
}
\references{Coming soon.}
\author{Johannes Freudenberg, Xiangdong Liu, Mario Medvedovic.}
\note{
}
\seealso{\code{\link{hclust}}, \code{\link[gimmR:gimmR-package]{gimmR}}, \code{\link[GO:GO-package]{GO}}, \code{\link[KEGG:KEGG-package]{KEGG}}, \code{\link{r2cdt}}, \code{\link{call.treeview}}}
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
#    "CpGislands", "GO", "KEGG"), species = "Rn", callTreeView = TRUE, sampleDesc = trt)
generateTreeViewFiles(gimmOut, cclust = NA, verbose = FALSE, functionalCategories=c("geneRIFs",
    "CpGislands", "GO", "KEGG"), species = "Rn", callTreeView = FALSE, sampleDesc = trt)

#geneList enrichment
geneList <- gimmOut$clustData[,1]
require(org.Rn.eg.db)
allGenes <- unique(keys(org.Rn.egSYMBOL)) #one should really use the list of
                                          #genes represented on the microarray instead
res <- geneListEnrichment(geneList, allGenes, functionalCategories = "GO",
	species = "Rn", sigFDR = 0.01, maxGenesInCategory = 10000)

#using primary geneset
data(cMap)
data(gimmOut)
#download.file("ftp://ftp.ncbi.nih.gov/pub/HomoloGene/current/homologene.data",destfile="homologene.data",mode="wb")
#pValues.rat <- convertGeneTable(cMap$cMapPvalues, fromSpecies="h", toSpecies="r")
#generateTreeViewFiles(gimmOut, functionalCategories=list(cMap=list(pValues.rat, cMap$cMapDescr)), species="Rn", bkgList=NULL)
#call.treeview()
}
\keyword{cluster}

