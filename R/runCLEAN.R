`runCLEAN` <-
function(data = NULL, 
		rclust = NULL, 
		functionalCategories = NULL, 
		bkgList = NULL,
		file = "cluster", 
		maxSize=1000, 
		minSize=10, 
		maxNumOfClusters=NULL,
		minGenesInCategory=10,
		maxGenesInCategory=1000,
		inBkg=TRUE,
		sigFDR = 0.1,
		species = NULL,
		estimateNullDistribution = FALSE,
		verbose=TRUE, 
		saveDataObjects = FALSE) {
	if(is.null(data)) { 
		if(file.access(paste(file, "RData", sep = "."))==0) 
			data <- paste(file, "RData", sep = ".")
		else if(file.access(paste(file, "cdt", sep = "."))==0)
			data <- paste(file, "cdt", sep = ".")
		else 
			data <- paste(file, "txt", sep = ".")
	}
	fromCdtFile <- FALSE
	if(is.character(data)) { #is a file name
		if(substr(data,nchar(data)-4,nchar(data)) == "RData") {
			e1 <- new.env()
			load(file = data, envir = e1)
			data <- get(ls(envir = e1, pattern = "gimmOut")[1], envir = e1)
		}
		else if(substr(data,nchar(data)-2,nchar(data)) == "cdt") {
			data <- cdt2r(data)
			fromCdtFile <- TRUE
		} else {
			data <- read.table(file = data, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
		}
	}
	if(is.matrix(data)) {
		cdt <- cbind(ID=rownames(data), NAME=as.character(rownames(data)), data) #rownames are assumed to be gene IDs
	} else if(is.data.frame(data)) {
		if(is.numeric(data[,2]) & #second column contains data
		length(grep("descr", colnames(data)[2], ignore.case=TRUE)) <= 0 &
		length(grep("name", colnames(data)[2], ignore.case=TRUE)) <= 0) {
			if(length(grep("id",colnames(data)[1], ignore.case=TRUE)) > 0 |   #if first column contains gene id
				length(grep("acc",colnames(data)[1], ignore.case=TRUE)) > 0 |
				!is.numeric(data[, 1])) {
				cdt <- cbind(ID=data[,1], NAME=as.character(data[,1]), data[,-1])
				warning("data[, 2] is numeric.  Assuming data[, 1] = ID column, data[, 2] data column")
			} else { #first column contains data
				cdt <- cbind(ID=rownames(data), NAME=as.character(rownames(data)), data) #rownames are assumed to be gene IDs
				warning("data[, 1:2] are numeric.  Assuming all columns to contain data.")
			}
	
		} else {
			cdt <- data
		}
	} else { #is gimmOut object
		cdt <- data$clustData
		if(is.null(rclust)) rclust <- data$hGClustData
	}
	if(is.null(rclust)) rclust <- "euclidean"  #default is hclust with euclidean distance
	run.CLEAN = TRUE
	if(length(rclust) == 1) {
		if(is.na(rclust)) {                        #no gene clustering specified
			run.CLEAN <- FALSE
		}
	}
	if(is.character(rclust)) {  #is filename or dist method
		if(substr(rclust,nchar(rclust)-2,nchar(rclust)) == "gtr") { #import gtr file
			rclust <- gtr2r(rclust)
		} else if(substr(rclust,nchar(rclust)-1,nchar(rclust)) == "zm") { #import zm file
			rclust <-data.matrix(scan(rclust))
			rclust <- 1 - matrix(rclust,sqrt(length(rclust)))
		} else {
			rclust <- dist(cdt[,-(1:2)], method = rclust) #compute distance matrix from data
			fromCdtFile <- FALSE
		}
	}
	if(is.numeric(rclust)) {
		gtr <- hclust(as.dist(rclust), method="average") 
		gtr$labels <- 1:nrow(cdt) 
	}
	else gtr <- rclust
	if(fromCdtFile) cdt <- cdt[order(gtr$order),]
	if(is.null(gtr$labels)) gtr$labels <- 1:nrow(cdt)
	if(is.character(gtr$labels)) {
		gtr$labels <- 1:nrow(cdt)
		warning("rclust$labels are character, reassigned as rclust$labels<-1:nrow(cdt)")
	}
	if(length(intersect(gtr$labels, 1:nrow(cdt))) < nrow(cdt)) {
		gtr$labels <- 1:nrow(cdt)
		warning("rclust$labels incompatible, reassigned as rclust$labels<-1:nrow(cdt)")
	}
	if(is.null(functionalCategories)) functionalCategories <- "GO"  #default is GO categories
	if(length(functionalCategories) == 1) {
		if(is.na(functionalCategories)) {  #no functional categories specified
			run.CLEAN <- FALSE
		}
	}
	if(run.CLEAN) {
		if(!is.list(functionalCategories)) functionalCategories <- as.list(unique(functionalCategories))
		CategoryID2GeneID <- list()
		CategoryID2Desc <- list()
		run.CLEAN <- rep(TRUE, length(functionalCategories))
		for(i in 1:length(functionalCategories)) {
			if(is.character(functionalCategories[[i]])) {
				if(substr(functionalCategories[[i]],nchar(functionalCategories[[i]])-4,nchar(functionalCategories[[i]])) == "RData") {
					e1 <- new.env()
					load(file = functionalCategories[[i]], envir = e1)
					if(length(ls(envir = e1, pattern = "2[G,g][E,e][N,n][E,e]")) > 0) {
						l <- list(get(ls(envir = e1, pattern = "2[G,g][E,e][N,n][E,e]")[1], envir = e1))
						names(l) <- names(functionalCategories)[i]
						if(is.null(names(l))) names(l) <- substr(functionalCategories[[i]], 1, nchar(functionalCategories[[i]])-6)
						CategoryID2GeneID <- c(CategoryID2GeneID, l)
					} else {
						run.CLEAN[i] <- FALSE
						CategoryID2GeneID <- c(CategoryID2GeneID, list(NA))
						warning(paste("No functional categories found in file", functionalCategories[[i]], "-  No functional clustering annotation generated."))
					}
					if(length(ls(envir = e1, pattern = "2[D,d][E,e][S,s][C,c]")) > 0) {
						l <- list(get(ls(envir = e1, pattern = "2[D,d][E,e][S,s][C,c]")[1], envir = e1))
						names(l) <- names(functionalCategories)[i]
						if(is.null(names(l))) names(l) <- substr(functionalCategories[[i]], 1, nchar(functionalCategories[[i]])-6)
						CategoryID2Desc <- c(CategoryID2Desc, l)
					} else CategoryID2Desc <- c(CategoryID2Desc, list(NA))
				} else {
					l1 <- getFunctionalCategories(functionalCategories[[i]], species)
					if(length(l1) > 0) {
						l <- l1[[1]][1]
						names(l) <- functionalCategories[[i]]
						CategoryID2GeneID <- c(CategoryID2GeneID, l)
						if(length(l1[[1]]) > 1) {
							l <- l1[[1]][2]
							names(l) <- functionalCategories[[i]]
							CategoryID2Desc <- c(CategoryID2Desc, l)
						} 
						else CategoryID2Desc <- c(CategoryID2Desc, list(NA))
					} else {
						run.CLEAN[i] <- FALSE
						CategoryID2GeneID <- c(CategoryID2GeneID, list(NA))
						CategoryID2Desc <- c(CategoryID2Desc, list(NA))
						warning(paste("functional categories", functionalCategories[[i]], "not found.  No functional clustering annotation generated."))
					}
				}
			} else {   #functionalCategories[[i]] is a list
				l <- functionalCategories[[i]][1]
				names(l) <- names(functionalCategories)[i]
				CategoryID2GeneID <- c(CategoryID2GeneID, l)
				if(length(functionalCategories[[i]]) > 1) {
					l <- functionalCategories[[i]][2]
					names(l) <- names(functionalCategories)[i]
					CategoryID2Desc <- c(CategoryID2Desc, l)
				}
				else CategoryID2Desc <- c(CategoryID2Desc, list(NA))
			}
		}
	}
	if(any(run.CLEAN)) {
		fClustAnnotations <- list()
		for(i in 1:length(run.CLEAN)) {
			if(run.CLEAN[i]) {
				l <- CLEAN(tree = gtr, allGenes=cdt[,1], CategoryID2GeneID=CategoryID2GeneID[[i]],
					CategoryID2Desc = CategoryID2Desc[i], refGenes = unique(bkgList), maxSize=maxSize,
					minSize=minSize, maxNumOfClusters=maxNumOfClusters, minGenesInCategory=minGenesInCategory, 
					maxGenesInCategory=maxGenesInCategory, verbose=verbose,inBkg=inBkg, sigFDR=sigFDR, 
					saveDataObjects=saveDataObjects, estimateNullDistribution = estimateNullDistribution)
				names(l) <- names(CategoryID2GeneID)[i]
				fClustAnnotations <- c(fClustAnnotations, l)
			}
		}
	}
	if(saveDataObjects) { 
		save(fClustAnnotations, file = paste(file,"RData",sep ="."))
	}
	invisible(fClustAnnotations)
}

