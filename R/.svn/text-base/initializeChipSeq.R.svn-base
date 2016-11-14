`initializeChipSeq`<-
function(chip.seq,genome.id="hg18",tablename="refGene",transcriptDB=NULL,species="Hs",distanceRange=c(-1e+06,1e+06),file=NULL)
{
	require(GenomicFeatures)
	colnames(chip.seq)<-c("Chromosome","Start","End","Des","Score")
	require(paste("org",species,"eg.db",sep="."),character.only=T)
	
	if(is.null(transcriptDB))genomeDB = makeTranscriptDbFromUCSC(genome=genome.id,tablename=tablename)
	else genomeDB<-loadFeatures(transcriptDB)
	allchr<-unique(chip.seq[,1])

	refseqs = names(as.list(get(paste("org",species,"egREFSEQ2EG",sep="."))))
	exps = matrix(NA,ncol=5)
	colnames(exps)<-c("RefID","GeneID","Chromosome","TSSDist","Score")

	for(chr in allchr)
	{
	  cat("Annotating ",chr,"\n")
	  genome<-transcripts(genomeDB,val<-list(tx_name=refseqs,tx_chrom=chr))
	  genome.starts=start(genome)
	  genome.ends=end(genome)
	  genome.strands=as.vector(strand(genome))
	  genome.0.starts=genome.starts
	  genome.0.ends=genome.ends
	  genome.0.starts[genome.strands=="+"]<-genome.starts[genome.strands=="+"]+distanceRange[1]
	  genome.0.starts[genome.0.starts<0]<-0
	  genome.0.ends[genome.strands=="+"]<-genome.starts[genome.strands=="+"]+distanceRange[2]
	  genome.0.starts[genome.strands=="-"]<-genome.ends[genome.strands=="-"]+distanceRange[1]
	  genome.0.starts[genome.0.starts<0]<-0
	  genome.0.ends[genome.strands=="-"]<-genome.ends[genome.strands=="-"]+distanceRange[2]
	  gene.id = unlist(mget(as.character(values(genome)[,"tx_name"]),get(paste("org",species,"egREFSEQ2EG",sep="."))))


	  tree <- IntervalTree(IRanges(genome.0.starts,genome.0.ends))
	  temp<-chip.seq[chip.seq$Chromosome==chr,]
	  tables<-matchMatrix(findOverlaps(query=IRanges(temp$Start,temp$End),tree,type="within",select="all"))
	  distances<-genome.starts[tables[,2]]
	  distances[is.element(tables[,2],which(genome.strands=="-"))]<-genome.ends[tables[is.element(tables[,2],which(genome.strands=="-")),2]]
	  starts1<-temp[tables[,1],"Start"]
	  starts2<-temp[tables[,1],"End"]
	  distances<-as.integer(abs((starts1+starts2)/2-distances))
	  tables<-cbind(tables,Dist=distances,Score=temp[tables[,1],"Score"])
	  rowID<-paste(values(genome)[tables[,2],"tx_name"],tables[,3],tables[,1],sep="-")
	  tables<-cbind(values(genome)[tables[,2],"tx_name"],gene.id[tables[,2]],as.character(rep(chr,dim(tables)[1])),tables[,c(3,4)])
	  rownames(tables)<-rowID
	  colnames(tables)<-c("RefID","GeneID","Chromosome","TSSDist","Score")
	  tables<-data.frame(tables,stringsAsFactors=FALSE)
	  exps<-rbind(exps,tables)
	  rm(tables)
	}
	exps<-exps[!is.na(exps[,1]),]
	exps
}
