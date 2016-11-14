`cdt2r` <-
function(cdtFile) {
	if (!file.exists(cdtFile)) {
		warning(paste("File not found: ", cdtFile, ". Trying ", cdtFile, ".cdt", sep=""))
		cdtFile <- paste(cdtFile, "cdt", sep = ".")
	}
	cdt <- read.table(cdtFile, sep = "\t", stringsAsFactors=FALSE, header = TRUE, quote = "")
	data.start.row <- which(cdt[,1] == "EWEIGHT") + 1
	data.start.col <- which(colnames(cdt) == "GWEIGHT") + 1
	if(length(data.start.row) < 1) { #there was no "EWEIGHT" row
		data.start.row <- 1
	}
	if(length(data.start.col) < 1) { #there was no "GWEIGHT" column
		if(colnames(cdt)[1] == "GID") data.start.col <- 4
		else data.start.col <- 3
	}
	if(colnames(cdt)[1] == "GID"){
		ID.col <- 2
		NAME.col <- 3
	} else {
		ID.col <- 1
		NAME.col <- 2
	}
	data.last.row <- dim(cdt)[1]
	data.last.col <- dim(cdt)[2]
	cdt <- cdt[data.start.row:data.last.row, c(ID.col,NAME.col,data.start.col:data.last.col)]
	for(i in 3:dim(cdt)[2]) cdt[, i]<-as.numeric(as.character(cdt[, i]))
	cdt
}

