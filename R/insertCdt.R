`insertCdt` <-
function(cdtFile, insertVector,tree=NULL,columnName="FunctionalClustering",newCdtFile=NULL){
	if(length(tree)!=0) insertVector<-insertVector[tree$order]
	cdtData<-read.table(file=cdtFile,header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE, check.names=FALSE)
	if(length(insertVector) == 1) {
		if(insertVector == "") cdtData <- cbind(cdtData, "")
	} else {
		i<-1
		while(length(grep("GENE",cdtData[i,1]))==0) {insertVector<-c(1, insertVector); i<-i+1}
		cdtData$FunctionalClustering<-insertVector
		if(columnName!="FunctionalClustering") names(cdtData)[dim(cdtData)[2]]<-columnName
	}
	if(length(newCdtFile)!=0) write.table(cdtData,file=newCdtFile,sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
	else write.table(cdtData,file=cdtFile,sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
	invisible(cdtData)
}

