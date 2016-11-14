\name{RandomSet}
\alias{RandomSet}
\alias{RS}
\title{ Testing GO terms, KEGG pathways, and other categories with Random Set method }
\description{
    This function uses the Random Set method to test for enriched biological
    categories in gene expression data.
}
\usage{
RandomSet(sigvals, geneids, database = "GO", functionalCategories = NULL,
    species = "Hs", min.g = 10, minGenesInCategory = NULL, max.g = NA,
    maxGenesInCategory = NULL, sig.cutoff = 0.1, sigFDR = NULL,
	averageMultipleProbes = TRUE, allGenesInCategoriesAsBackground=TRUE,
	two.sided = FALSE, na.rm=TRUE, verbose = TRUE)
RS(sigvals, geneids, database = "GO", functionalCategories = NULL,
    species = "Hs", min.g = 10, minGenesInCategory = NULL, max.g = NA,
    maxGenesInCategory = NULL, sig.cutoff = 0.1, sigFDR = NULL,
	averageMultipleProbes = TRUE, allGenesInCategoriesAsBackground=FALSE,
	two.sided = FALSE, na.rm=TRUE, verbose = TRUE)
}
\arguments{
  \item{sigvals}{A vector of p-values, same length and order as "geneids"}
  \item{geneids}{A vector of Entrez gene IDs, may contain duplicates and missing values}
  \item{min.g}{Deprecated. Please use `minGenesInCategory` instead.}
  \item{minGenesInCategory}{The minimum number of unique gene IDs analyzed in category to be tested, if NULL it is set to 0.}
  \item{max.g}{Deprecated. Please use `minGenesInCategory` instead.}
  \item{maxGenesInCategory}{The maximum number of unique gene IDs analyzed in category to be tested, if NULL it is set to Inf.}
  \item{sig.cutoff}{Deprecated. Please use `sigFDR` instead.}
  \item{sigFDR}{Categories with FDR <= sigFDR will be returned.  If NULL it is set to 1.}
  \item{database}{Deprecated. Please use `functionalCategories` instead.}
  \item{functionalCategories}{Functional categories to be tested- currently, options include "GO", "KEGG" and various other categories, default = "GO".  Can be provided by function getFunctionalCategories().}
  \item{species}{Species to further specify database, human="Hs", mouse="Mm", rat="Rn", etc. Default ="Hs".}
  \item{averageMultipleProbes}{If multiple probes per geneID, the (geometric) mean is computed.}
  \item{allGenesInCategoriesAsBackground}{If TRUE, all genes in a list of functional categories (e.g. "GO") will be used as background gene list, and computations are limited to intersection of background and genes in "geneids" paramter.}
  \item{two.sided}{If TRUE, the two-sided p-value is computed.}
  \item{na.rm}{If TRUE, potential NAs and NaNs in sigvals are removed before computing the Random Set statistic.}
  \item{verbose}{If TRUE, produces lots of output.}
}
\details{
This function uses the Random Set method (Newton et al., 2007) to test for enriched biological
    categories in gene expression data.

}
\value{
Object is a dataframe with the following columns:
category ID 
category description 
n.genes - number of genes overlapping between gene list and category
zScore - the random set z-score
p-value - the corresponding (one-sided or two-sided) p-value
FDR - False Discovery Rate (Benjamini & Hochberg, 1995)
}
\references{Newton, 2007. Annals App Stat `Random Set Methods Identify Distinct Aspects of the Enrichment Signal in Gene Set Analysis`}
\author{Johannes Freudenberg}
\seealso{\code{\link[CLEAN]{LRpath}}, \code{\link[GO.db:GO.db-package]{GO.db}}, \code{\link[KEGG.db:KEGG.db-package]{KEGG.db}}, \code{\link[gimmR:gimmR-package]{gimmR}}}
\examples{
data(gimmOut)
p <- rbeta(94, 0.5, 2)
res <- RandomSet(sigvals=p, geneids=gimmOut$clustData[,1], functionalCategories=c("GO", "KEGG"), species="Rn")
names(res)
head(res[[1]])
}
\keyword{regression}
\keyword{htest}
