\name{enrichmentGSEA}
\title{Performing GSEA analysis with the fgsea package}
\description{
    This function uses the GSEA methodology to perform enrichment analysis of transciptional signatures. The statistical significance is calculated by random gene permutations.
}
\usage{
enrichmentGSEA(sigvals,geneids,diffexp=NULL, functionalCategories=c("KEGG","GO","tfacts","dorotheatfs","hallmark"), directions=c("both"),species="Hs", excel=F,sigFDR=0.1,nperm=100000)
}
\arguments{
  \item{sigvals}{A vector of p-values, same length and order as "geneids".}
  \item{geneids}{A vector of Entrez gene IDs, may contain duplicates and missing values.}
  \item{diffexp}{log2 differential expression levels.}
  \item{functionalCategories}{Functional categories to be tested- currently, options include "GO", "KEGG" and various other categories, default = "GO".  Can be provided by function getFunctionalCategories().}
  \item{directions}{Specifies which of the directions are changes is used in the enrichment values can be "both" (default), "up", "down" and the combinations of the three.}
  \item{species}{Species to further specify database, human="Hs", mouse="Mm", rat="Rn", etc. Default ="Hs".}
  \item{excel}{Specifies whether the results should be organize in excell sheets. Default is FALSE. If results are big it can cause Java issues.}
  \item{sigFDR}{Categories with FDR <= sigFDR will be returned.  If NULL it is set to 1.}
  \item{nperm}{Number of random gene samples used for calculating p-values.}
}
\details{
  This function uses the GSEA methodology as implemented in the fgsea package to perform enrichment analysis of transciptional signatures. The statistical significance is calculated by random gene permutations.
}
\value{
  Object is a list of lists of dataframes for each of the functional categories used in the analysis. For functional category a list of dataframes of results corresponding to directions (both, up, down) specified is generated.
}
\references{}
\author{Mario Medvedovic}
\seealso{\code{\link[CLEAN.Hs]{CLEAN.Hs}}, \code{\link[fgsea]{fgsea}}}
\examples{
data(cMap)
de<-rnorm(11911)
res <- enrichmentGSEA(sigvals=cMap$cMapPvalues[,2],diffexp=de, geneids=cMap$cMapPvalues[,1], functionalCategories=c("KEGG","tfacts"), species="Hs",sigFDR=.5,nperm=10000)
names(res)
names(res$both)
head(res$both$KEGG)
}
\keyword{GSEA}
\keyword{gsea}
\keyword{enrichment}
