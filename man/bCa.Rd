\name{bCa}
\alias{bCa}
\docType{data}
\title{Sets of p-values derived from two breast cancer data sets.}
\description{
  Sample data set for GRS().  P-values were calculated using the t.test() function and ER status as the discriminitory variable.
}
\usage{data(bCa)}
\format{
  The format is:
List of 2 data.frames where 1st column contains gene IDs and 2nd column contains p-values
}
\source{
   "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE11121"
   "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE3494"
}
\references{
      "Schmidt M, Böhm D, von Törne C, Steiner E et al.: The humoral immune system has a key prognostic impact in node-negative breast cancer. Cancer Res 2008 Jul 1;68(13):5405-13. PMID: 18593943"
	  "Miller LD, Smeds J, George J, Vega VB et al. An expression signature for p53 status in human breast cancer predicts mutation status, transcriptional effects, and patient survival. Proc Natl Acad Sci U S A 2005 Sep 20;102(38):13550-5. PMID: 16141321"
}
\examples{
data(bCa)
## not run:
## geneTable <- convertGeneTable(bCa$Schmidt.pValues)
}
\keyword{datasets}
