\name{LRpath}
\alias{LRpath}
\title{ Testing GO terms or KEGG with logistic regression }
\description{
    This function uses logistic regression to test for enriched biological
    categories in gene expression data. Our method models the probability
    of a randomly selected gene belonging to a specific category given the
    significance level of that gene.  For categories significantly affected by
    the experimental condition, this probability will increase as the significance
    statistic increases.  Categories with significant p-values and positive slope
    coefficients are enriched with differentially expressed genes.
}
\usage{
LRpath(sigvals, geneids, min.g=10, max.g=NA, sig.cutoff=0.05, database="GO", functionalCategories=NULL, 
	odds.min.max=c(0.001,0.5), species="Hs") 
}
\arguments{
  \item{sigvals}{A vector of p-values, same length and order as "geneids"}
  \item{geneids}{A vector of Entrez gene IDs, may contain duplicates and missing values}
  \item{min.g}{The minimum number of unique gene IDs analyzed in category to be tested, default = 10}
  \item{max.g}{The maximum number of unique gene IDs analyzed in category to be tested, default = NA (99999)}
  \item{sig.cutoff}{Entrez gene IDs in each category with p-values<sig.cutoff will be returned, default = 0.05}
  \item{database}{Deprecated. Please use `functionalCategories` instead.}
  \item{functionalCategories}{Functional categories to be tested- currently, options include "GO", "KEGG" and various other categories, default = "GO".  Can be provided by function getFunctionalCategories().}
  \item{odds.min.max}{Lower and upper p-values to be used for odds ratio calculation, default= c(0.001, 0.5)}
  \item{species}{Species to further specify database, human="Hs", mouse="Mm", rat="Rn", etc. Default ="Hs".}
}
\details{
LRpath: testing GO terms or KEGG with logistic regression
Written by: Maureen Sartor, University of Cincinnati, 2008
   
This function uses logistic regression to test for enriched biological
categories in gene expression data. Our method models the probability
of a randomly selected gene belonging to a specific category given the
significance level of that gene.  For categories significantly affected by
the experimental condition, this probability will increase as the significance
statistic increases.  Categories with significant p-values and positive slope
coefficients are enriched with differentially expressed genes.

Please acknowledge your use of LRpath in publications by referencing the Sartor et al. (2009) paper.

}
\value{
Object is a dataframe with the following columns:
GO or KEGG ID   - posterior t-value for IBMT
GO or KEGG term     - name of category
Ontology (only for GO) - BP, MF, or CC
n.genes - Number of unique Entrez Gene IDs in category
coeff - coefficient of slope (only GO terms with positive values are significant)
odds.ratio -  Odds ratio, as measure of strenth of enrichment
p.value -  P-value that slope does not equal zero (that term is enriched)
FDR - False Discovery Rate (Benjamini & Hochberg, 1995)
sig.genes   - comma separated Entrez gene ids in category with p-value<"sig.cutoff"
}
\references{Sartor MA, Leikauf GD, Medvedovic M. 2009. LRpath: A logistic regression approach for identifying enriched 
  biological groups in gene expression data. Bioinformatics 25(2):211-7.}
\author{Maureen Sartor}
\seealso{\code{\link[stats]{glm}}, \code{\link[GO.db:GO.db-package]{GO.db}}, \code{\link[KEGG.db:KEGG.db-package]{KEGG.db}}, \code{\link[gimmR:gimmR-package]{gimmR}}}
\examples{
data(gimmOut)
p <- rbeta(94, 0.5, 2)
LRpath(sigvals=p, geneids=gimmOut$clustData[,1], functionalCategories="GO", species="Rn")
}
\keyword{regression}
\keyword{htest}
