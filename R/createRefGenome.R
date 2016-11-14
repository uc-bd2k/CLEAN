`createRefGenome`<-
function(filename=NULL,genome.id="hg18",tablename="refGene")
{
	require(GenomicFeatures)
	genomeDB = makeTranscriptDbFromUCSC(genome=genome.id,tablename=tablename)
	if(!is.null(filename)) saveFeatures(genomeDB, file=paste(filename,".sqlite",sep=""))
	genomeDB
}