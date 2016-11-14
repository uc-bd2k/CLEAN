\name{genesInEnrichedCategories}
\alias{genesInEnrichedCategories}
\title{Function to intersect a gene list with enriched categories}
\description{
This function takes a gene list and a list of IDs of (enriched) categories and returns a bolean matrix where each column indicates whether or not the gene is member of the corresponding categories.
}
\usage{
genesInEnrichedCategories(categoryIDs, geneList, funcCategories = NULL, species = NULL)
}
\arguments{
  \item{categoryIDs}{A list of functional category IDs.}
  \item{geneList}{A gene list to be tested for enrichment of functional categories.}
  \item{funcCategories}{A collection of gene sets, representing functional categories.}
  \item{species}{Two letter description of the species to be used to generate gene ontology categories (e.g. "Hs" for human, "Mm" for mouse). This parameter is used when funcCategories = "GO"}
}
\details{
coming soon.
}
\value{
  A matrix with rows representing the functional categories and columns representing the genes.
}
\references{Coming soon.}
\author{Johannes Freudenberg, Xiangdong Liu, Mario Medvedovic.}
\note{
}
\seealso{\code{\link{geneListEnrichment}}}
\examples{
data(gimmOut)
require(CLEAN.Rn)
require(org.Rn.eg.db)

geneList <- gimmOut$clustData[,1]
allGenes <- unique(keys(org.Rn.egSYMBOL)) #one should really use the list of
                                          #genes represented on the microarray instead
res <- geneListEnrichment(geneList, allGenes, functionalCategories = "GO",
	species = "Rn", sigFDR = 0.01, maxGenesInCategory = 10000)

#display membership of each gene in top ten categories
genesInEnrichedCategories(res[1:10,1], geneList, funcCategories = "GO", species = "Rn")

}
\keyword{cluster}

