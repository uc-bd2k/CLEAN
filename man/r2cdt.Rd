\name{r2cdt}
\alias{r2cdt}
\alias{r2gtr}
\alias{r2atr}
\alias{r2fni}
\title{Functions to generate fTreeView files}
\description{
r2cdt, r2gtr, and r2atr are slightly modified functions from the package ctc.  

r2cdt generates a tab-delimited cdt file from a dataframe where the first column contains gene IDs, the second column contains gene names or descriptions, and the remaining columns represent samples.  Column names are assumed to be sample names.  In addition, corresponding hierarchical gene and sample clusterings are required to re-order genes and samples accordingly.  

r2gtr and r2atr generate tab-delimited TreeView files representing the hierarchical gene and sample clustering, respectively. 

r2fni generates the Functional Node Information file for the Java application fTreeView using an hierarchical gene clustering and corresponding functional Cluster Enrichment Analysis results.  This tab-delimited file is then displayed interactively as the user selects a node (i.e. gene cluster) in the gene tree.  
}
\usage{
r2cdt(hr, hc, data, labels = TRUE, description = TRUE, file = "cluster.cdt", dec = ".")
r2gtr(hr, file = "cluster.gtr", dec = ".")
r2atr(hc, file = "cluster.atr", dec = ".")
r2fni(hr, file = "cluster.fni", clustCategories = NULL, dec = ".")
}
\arguments{
  \item{hr}{An hclust object representing the gene clustering}
  \item{hc}{An hclust object representing the sample clustering}
  \item{data}{A dataframe where the first column may contain gene IDs, the second column may contain gene names or descriptions, and the remaining columns represent samples.  Column names are assumed to be sample names. }
  \item{labels}{If TRUE, the first column of data is assumed to contain gene IDs.}
  \item{description}{If TRUE, the second column of data is assumed to contain gene names or descriptions.}
  \item{file}{The name of the file to be generated.}
  \item{dec}{The character to be used for decimal points.}
  \item{clustCategories}{A list of functional Cluster Enrichment Analysis results. See funcClustAnnot for more details.}
}
\details{
}
\value{
Nothing is returned but a file with the specified name is created.
}
\references{ctc package: http://antoinelucas.free.fr/ctc/

fTreeView: http://eh3.uc.edu/ftreeview}
\author{Johannes Freudenberg}
\note{  
}
\seealso{\code{\link[ctc:ctc-package]{ctc}}, \code{\link{funcClustAnnot}}, \code{\link{hclust}}, \code{\link{call.treeview}}}
\examples{
data(gimmOut)
r2cdt(gimmOut$hGClustData, gimmOut$hSClustData, gimmOut$clustData)
r2gtr(gimmOut$hGClustData)
r2atr(gimmOut$hSClustData)
## not run
# call.treeview("cluster.cdt")
}
\keyword{ IO }
\keyword{ cluster }
