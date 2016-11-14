\name{cdt2r}
\alias{cdt2r}
\alias{atr2r}
\alias{gtr2r}
\title{Functions to import files in TreeView format into R.}
\description{
  cdt2r() imports a cdt file and converts it to a dataframe where rows represent genes and columns represent samples.  The first column of the dataframe contains the gene identifiers and the second column contains gene names or descriptions.  Column names represent sample names.
}
\usage{
cdt2r(cdtFile)
gtr2r(gtrFile)
atr2r(atrFile)
}
\arguments{
  \item{cdtFile}{A cdt file}
  \item{gtrFile}{A gtr file}
  \item{atrFile}{An atr file}
}
\details{
  cdt2r() reads a tab-delimited file in generalized cdt format where the first row contains sample names, followed by additional optional rows.  The first two or three columns contain gene IDs and descriptions, and optional additional columns, followed by gene expression data.  The function returns a dataframe with gene IDs and descriptions in the first two columns, gene expression data in the remaining columns, and column names representing samples. 
  gtr2r() and atr2r() read tree files representing the hierarchical gene tree or sample tree structure, respectively.  Both functions return an hclust object.
}
\value{
  cdt2r() returns a dataframe (see details)
  gtr2r() and atr2r() each return an hclust objects
}
\references{ coming soon }
\author{ Johannes Freudenberg, Xiangdong Liu, Mario Medvedovic }
\seealso{ \code{\link{funcClustAnnot}}, \code{\link{hclust}}, \code{\link{r2cdt}} }
\examples{
data(gimmOut)
generateTreeViewFiles(gimmOut$clustData, file="cluster") 
exprData <- cdt2r("cluster.cdt")
gclust <- gtr2r("cluster.gtr")
plot(as.dendrogram(gclust))
sclust <- atr2r("cluster.atr")
}
\keyword{ cluster }
\keyword{ IO }
