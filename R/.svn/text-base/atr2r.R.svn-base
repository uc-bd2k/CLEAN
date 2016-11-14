`atr2r` <-
function(atrFile){
	if (!file.exists(atrFile)) {
		warning(paste("File not found: ", atrFile, ". Trying ", atrFile, ".atr", sep=""))
		atrFile <- paste(atrFile, "atr", sep = ".")
	}
	colc<-c("character","character","character","numeric")
	atr<-read.table(atrFile,colClasses=colc)
	height<-1-atr[,4]
	clusterNum<-dim(atr)[[1]]
	merge<-matrix(0,clusterNum,2)
	clusterSize<-integer(clusterNum)
	for(i in 1:clusterNum) {
		for( j in 1:2) {
			tempStr<-atr[i,(j+1)]
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

