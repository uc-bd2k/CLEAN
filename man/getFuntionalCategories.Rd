\name{getFunctionalCategories}
\alias{getFunctionalCategories}
\title{Retrieve functional categories from species-specific library.}
\description{
  Load the species-specific library/ies and retrieve the functional catogories of type(s) "functionalCategories".
}
\usage{
getFunctionalCategories(functionalCategories, species)
}
\arguments{
  \item{functionalCategories}{A character vector of functional catogory types.}
  \item{species}{Species abreviation such as "Hs", "Mm", or "Rn"}
}
\details{
}
\value{
  A list containing the functional catories and their description (optional).
}
\references{}
\author{Johannes Freudenberg}
\note{ 
}
\seealso{\code{\link{runCLEAN}}}
\examples{
library(CLEAN.Rn)
data(gimmOut)
myFunctionalCategories <- getFunctionalCategories(CLEAN.Rn(), species = "Hs")
generateTreeViewFiles(gimmOut, functionalCategories = myFunctionalCategories)
}
\keyword{data}
