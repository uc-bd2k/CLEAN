`gtr2r` <-
function(gtrFile){
	if (!file.exists(gtrFile)) {
		warning(paste("File not found: ", gtrFile, ". Trying ", gtrFile, ".gtr", sep=""))
		gtrFile <- paste(gtrFile, "gtr", sep = ".")
	}
	colc<-c("character","character","character","numeric")
	gtr<-read.table(gtrFile,colClasses=colc)
	height<-1-gtr[,4]
	clusterNum<-dim(gtr)[[1]]
	merge<-matrix(0,clusterNum,2)
	clusterSize<-integer(clusterNum)
	for(i in 1:clusterNum) {
		for( j in 1:2) {
			tempStr<-gtr[i,(j+1)]
			positive<-substr(tempStr,1,4) == "NODE"
			id<-substr(tempStr,5,(nchar(tempStr)-1))
			if(positive) merge[i,j]<-as.integer(id)
			else  merge[i,j]<- -as.integer(id)
			if(merge[i,j] > 0) {
				clusterSize[i]<-clusterSize[i]+clusterSize[merge[i,j]]
			}else {
				clusterSize[i]<-clusterSize[i]+1
			}
		}
	}
	perm<-clusterNum
	allLeafs <- FALSE
	while(!allLeafs) {
		allLeafs <- all(perm < 0)
		if(!allLeafs) {
			index <- which(perm>=0)[1]
			newPerm <- as.vector(merge[perm[index], ])
			if(index > 1) newPerm <- c(perm[1:(index-1)], newPerm)
			if(index < length(perm)) newPerm <- c(newPerm, perm[(index+1):length(perm)])
			perm <- newPerm
		}
	}
	res <- list(merge=merge, height=height, order=-perm, labels = 1:length(perm), method = "average")
	class(res) <- "hclust"
	res
}

