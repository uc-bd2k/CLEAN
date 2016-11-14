\name{cMap}
\alias{cMap}
\docType{data}
\title{
A set of p-values derived from the Connectivity Map dataset.
}
\description{
A set of p-values derived from the Connectivity Map dataset. Intended for CLEAN with primary reference data, and for RS and GRS.

}
\usage{data(cMap)}
\format{
  The format is:
List of 2
 $ cMapPvalues: num [1:11911, 1:100] 10 100 1000 10000 10001 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:11911] "10_at" "100_at" "1000_at" "10000_at" ...
  .. ..$ : chr [1:100] "geneID" "cID_1" "cID_2" "cID_3" ...
 $ cMapDescr  : chr [1:99, 1:2] "cID_1" "cID_2" "cID_3" "cID_4" ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : NULL
  .. ..$ : chr [1:2] "ID" "descr"
}
\details{
}
\source{
}
\references{
Lamb J, Crawford ED, Peck D, Modell JW, Blat IC, Wrobel MJ, Lerner J, Brunet JP, Subramanian A, Ross KN, Reich M, Hieronymus H, Wei G, Armstrong SA, Haggarty SJ, Clemons PA, Wei R, Carr SA, Lander ES, Golub TR.  The Connectivity Map: using gene-expression signatures to connect small molecules, genes, and disease. Science. 2006 Sep 29;313(5795):1929-35.
}
\examples{
data(cMap)
data(gimmOut)
#download.file("ftp://ftp.ncbi.nih.gov/pub/HomoloGene/current/homologene.data",destfile="homologene.data",mode="wb")
#pValues.rat <- convertGeneTable(cMap$cMapPvalues, fromSpecies="h", toSpecies="r")
#generateTreeViewFiles(gimmOut, functionalCategories=list(cMap=list(pValues.rat, cMap$cMapDescr)), species="Rn", bkgList=NULL)
#call.treeview()
}
\keyword{datasets}
