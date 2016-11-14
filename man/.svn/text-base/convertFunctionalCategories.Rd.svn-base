\name{convertFunctionalCategories}
\alias{convertFunctionalCategories}
\alias{convertGeneTable}
\title{Function to convert Entrez gene ID-based functional categories from one species to another using the Homologene database}
\description{
    This function loads the Homologene table and a list of functional categories and converts Entrez gene IDs from one species to another.  Currently implemented species are "h" (Homo sapiens), "m" (Mus musculus), and "r" (Rattus norvegicus).  

	Use the following command to download the Homologene table into the working directory:
	> download.file("ftp://ftp.ncbi.nih.gov/pub/HomoloGene/current/homologene.data",destfile="homologene.data",mode="wb")

	This function was adapted from from convertSmc() of the PGSEA package.
}
\usage{
convertFunctionalCategories(functionalCategories, fromSpecies = "h", toSpecies = "r", homologTable = "./homologene.data")
convertGeneTable(geneTable, fromSpecies = "h", toSpecies = "r", homologTable = "./homologene.data")
}
\arguments{
  \item{functionalCategories}{A list of functional categories.}
  \item{geneTable}{A matrix or data.frame where the first column contains gene IDs and the second column contains p-values. }
  \item{fromSpecies}{Species the functional categories are currently defined for. Currently implemented species are "h" (Homo sapiens), "m"
  (Mus musculus), and "r" (Rattus norvegicus).}
  \item{toSpecies}{Species the functional categories are to be converted to. Currently implemented species are "h" (Homo sapiens), "m"
  (Mus musculus), and "r" (Rattus norvegicus).}
  \item{homologTable}{File name of the Homologene Table.}
}
\author{Johannes Freudenberg}
\seealso{\code{\link{runCLEAN}}, \code{\link{GRS}}}
\examples{
## not run
# download.file("ftp://ftp.ncbi.nih.gov/pub/HomoloGene/current/homologene.data",destfile="homologene.data",mode="wb")
# library(CLEAN.Hs)
# data(geneRIFs.Hs)
# geneRIFs.Rn <- geneRIFs.Hs
# geneRIFs.Rn$RIFterm2Genes <- convertFunctionalCategories(geneRIFs.Hs$RIFterm2Genes)
#
# data(bCa)
# geneTable <- convertGeneTable(bCa$Schmidt.pValues)
}
\keyword{ IO }

