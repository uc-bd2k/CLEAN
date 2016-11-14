\name{getBinarySubTree}
\alias{getBinarySubTree}
\title{Get the left or right subtree of a given hierarchical clustering.}
\description{
  Reads a hierarchical clustering and returns the left or right child of the root node.
}
\usage{
getBinarySubTree(tree, right = FALSE)
}
\arguments{
  \item{tree}{An hclust object}
  \item{right}{If TRUE, return the right subtree, i.e. the right child of the root node.}
}
\details{
}
\value{
  The subtree as an hclust object.
}
\references{}
\author{Johannes Freudenberg, Xiangdong Liu, Mario Medvedovic}
\note{ 
}
\seealso{\code{\link{hclust}}}
\examples{
hc <- hclust(dist(matrix(1:20,5)))
hc$labels <- 1:5
plot(hc)
hc.sub <- getBinarySubTree(hc, TRUE)
plot(hc.sub)
}
\keyword{cluster }
