\name{chipSeqWeightedSum}
\alias{chipSeqWeightedSum}
\title{Calculate binding strengh by summing weighted chip-seq peaks.}
\description{
  This function is intended to calculate chip-seq binding score and associated probability for each gene. The input chip-seq list contains all peaks found in up/down stream within certain range whose weights contribute to final score is assigned based on the its distance to transcription start site.
}
\usage{
chipSeqWeightedSum(ChipList,file=NULL,verbose=FALSE)
}
\arguments{
  \item{ChipList}{ A matrix contains all peaks for each genes with each row represent one peak and the first column contains ref-seq ID, second column contain entrez ID, third column contains chromosome ID, fourth column contains distance from the center of the peak to transcription start site and the fifth column contains the peak's strength}
  \item{file}{ File name and corresponding path to save output}
}
\details{
  This function is intended to calculate chip-seq binding score and associated probability for each gene.
}
\value{
  This function returns a matrix with first column contains the entrez ID, the second column contains the score and the third column contains the associated probability. 
}
\references{ Publication forthcoming. }
\author{Zhen Hu }
\seealso{}
\examples{
data(nMyc)
RefTable<-createRefGenome("mm8_refGene",genome.id = "mm8", tablename = "refGene")
ChipSeq<-annotateChipSeqPeaks(n.Myc,transcriptDB=RefTable,species="Mm",distanceRange=c(-1e+06,1e+06),file=NULL)
res<-chipSeqWeightedSum(ChipSeq,verbose=T)

}
\keyword{ chip-seq }
\keyword{ EM algorithm }
