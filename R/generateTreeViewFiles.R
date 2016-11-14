`generateTreeViewFiles` <-
function(data = NULL, 
		rclust = NULL, 
		cclust = NULL, 
		functionalCategories = NULL, 
		fClustAnnotations = NULL,
		bkgList = NULL,
		file = "cluster", 
		maxSize=1000, 
		minSize=10, 
		maxNumOfClusters=NULL,
		minGenesInCategory=10,
		maxGenesInCategory=1000,
		inBkg=TRUE,
		sigFDR = 0.1,
		addCLEANscore2cdt = TRUE,
		species = NULL,
		estimateNullDistribution = FALSE,
		dec='.',
		verbose=TRUE, 
		saveDataObjects = FALSE,
		sampleDesc = NULL,
		callTreeView = FALSE) {
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
			if(is.null(cclust)) cclust <- data$hSClustData
	}
	make.gtr <- TRUE
	if(is.null(rclust)) rclust <- "euclidean"  #default is hclust with euclidean distance
	if(length(rclust) == 1) {
		if(is.na(rclust)) {                        #no gtr file desired but need dummy file
			rclust <- hclust(dist(1:(nrow(cdt))))
			rclust$order<-rclust$order[order(rclust$order)]
			make.gtr <- FALSE
			fromCdtFile <- FALSE
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
		warning("rclust$labels incompatible, reassigned to rclust$labels<-1:nrow(cdt)")
	}
	make.atr <- TRUE
	if(is.null(cclust)) cclust <- 1-cor(cdt[,-(1:2)]) #default: hclust w/ pearson corr as similarity measure
	if(length(cclust) == 1) {
		if(is.na(cclust)) {                               #no atr file desired but need dummy file
				cclust <- hclust(dist(1:(ncol(cdt)-2))) 
				if(is.null(sampleDesc)) {
					cclust$order<-cclust$order[order(cclust$order)]
				} else {
					cclust$order<-(1:length(cclust$order))[order(sampleDesc)]
				}
				make.atr <- FALSE
		}
	}
	if(is.character(cclust)) {  #is filename or dist method
		if(substr(cclust,nchar(cclust)-2,nchar(cclust)) == "atr") { #import atr file
			cclust <- atr2r(cclust)
		} else if(substr(cclust,nchar(cclust)-1,nchar(cclust)) == "zm" | 
					substr(cclust,nchar(cclust)-2,nchar(cclust)) == "zm2") { #import zm file
			cclust <-data.matrix(scan(cclust))
			cclust <- 1 - matrix(cclust,sqrt(length(cclust)))
		} else {
			cclust <- dist(t(cdt[,-(1:2)]), method = cclust) #compute distance matrix from data
		}
	}
	if(is.numeric(cclust)) atr <- hclust(as.dist(cclust), method = "average")
	else atr <- cclust

	runCLEAN <- TRUE
	if(!is.null(fClustAnnotations)){
		runCLEAN <- FALSE
		if(is.character(fClustAnnotations)){ #is vector of (a) file name(s) for .RData object
		    fClustAnnotFiles <- fClustAnnotations
			fClustAnnotations <- list()
			for(i in 1:length(fClustAnnotFiles)) {
				e1 <- new.env()
				if(file.exists(fClustAnnotFiles[i])) {
					load(file = fClustAnnotFiles[i], envir = e1)
					fClustAnnotations<- c(fClustAnnotations, get("fClustAnnotations", envir = e1))
				} else {
					warning("File ", fClustAnnotFiles[i], " not found. Functional Cluster Annotations not loaded.")
				}
			}
		}
	}
	make.fni <- TRUE
	if(is.null(functionalCategories)) functionalCategories <- "GO"  #default is GO categories
	if(length(functionalCategories) == 1) {
		if(is.na(functionalCategories)) {  #no fni file desired
			runCLEAN <- FALSE
			make.fni <- FALSE
		}
	}
	if(runCLEAN) {
		if(!is.list(functionalCategories)) functionalCategories <- as.list(unique(functionalCategories))
		CategoryID2GeneID <- list()
		CategoryID2Desc <- list()
		runCLEAN <- rep(TRUE, length(functionalCategories))
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
						runCLEAN[i] <- FALSE
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
						} else CategoryID2Desc <- c(CategoryID2Desc, list(NA))
					} else {
						runCLEAN[i] <- FALSE
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
	r2cdt(gtr,atr,cdt,labels=TRUE, description=TRUE,file=paste(file, "cdt",sep="."),dec=dec)
    if(!make.atr & !is.null(sampleDesc)) {
		cdtFromFile <- read.table(file=paste(file, "cdt",sep="."), header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE, check.names=FALSE)
		data.start.col <- which(cdtFromFile[1,] == "GWEIGHT") + 1
		if (length(data.start.col) < 1) {
			if (colnames(cdt)[1] == "GID")
				data.start.col <- 4
			else data.start.col <- 3
		}
		newCdt <- cdtFromFile[,1:(data.start.col-1)]
		l <- data.start.col
		lSampleDesc <- levels(as.factor(sampleDesc))	
		for(f in lSampleDesc[-length(lSampleDesc)]) {
			l <- l + sum(sampleDesc == f)
			newCdt <- cbind(newCdt, cdtFromFile[,data.start.col:(l-1)], "")
			data.start.col <- l
		}
		l <- l + sum(sampleDesc == lSampleDesc[length(lSampleDesc)])
		newCdt <- cbind(newCdt, cdtFromFile[,data.start.col:(l-1)])
		write.table(newCdt, file=paste(file, "cdt",sep="."), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
	}
	if(make.atr) 
		r2atr(atr,file=paste(file,"atr",sep ="."), dec=dec)
	if(make.gtr) {
		r2gtr(gtr,file=paste(file,"gtr",sep ="."), dec=dec)
		if(make.fni) {
			if(any(runCLEAN)) {
				fClustAnnotations <- list()
				for(i in 1:length(runCLEAN)) {
					if(runCLEAN[i]) {
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
			if(length(fClustAnnotations) > 0) {
				if(addCLEANscore2cdt) {
					cwCLEANscore.all <- rep(0, dim(cdt)[1])
					for(i in 1:length(fClustAnnotations)) {
						#if(any(fClustAnnotations[[i]]$minFDR > 0)) 
						#	minFDR <- pruneMinFDR(fClustAnnotations[[i]])
						#else 
						#	minFDR <- fClustAnnotations[[i]]$sminFDR
						cwCLEANscore <- fClustAnnotations[[i]]$cwCLEANscore
						if(!is.null(cwCLEANscore)) {
							cwCLEANscore[cwCLEANscore > 300] <- 300
							sigLevels <- c(1, 3, 5, 10, 15, 20, 50, 100, 150, 200, 250, 299)
							if(max(cwCLEANscore) > 0) {
								newCdt<-insertCdt(cdtFile=paste(file, "cdt", sep = "."), insertVector = "", tree=gtr,columnName="", newCdtFile=NULL) #insert blank column
								for(sigLevel in sigLevels[sigLevels < max(cwCLEANscore)]) {
									newCdt<-insertCdt(cdtFile=paste(file, "cdt", sep = "."), insertVector=cwCLEANscore-sigLevel, 
										tree=gtr,columnName=paste(names(fClustAnnotations)[i], ": CLEANscore-", sigLevel, sep = ""), newCdtFile=NULL)
								}
							}
							cwCLEANscore.all <- apply(cbind(cwCLEANscore.all,cwCLEANscore), 1, max)
						}
					}
					if(length(fClustAnnotations) > 1 & max(cwCLEANscore.all) > 0) {
						sigLevels <- c(1, 3, 5, 10, 15, 20, 50, 100, 150, 200, 250, 299)
						newCdt<-insertCdt(cdtFile=paste(file, "cdt", sep = "."), insertVector = "", tree=gtr,columnName="", newCdtFile=NULL) #insert blank column
						for(sigLevel in sigLevels[sigLevels < max(cwCLEANscore.all)]) {
							newCdt<-insertCdt(cdtFile=paste(file, "cdt", sep = "."), insertVector=cwCLEANscore.all-sigLevel, 
								tree=gtr,columnName=paste("Overall CLEANscore-", sigLevel, sep = ""), newCdtFile=NULL)
						}
					}
				}
				r2fni(gtr, clustCategories = fClustAnnotations, file=paste(file, "fni", sep = "."))
			}
		}
	}
	if(callTreeView) call.treeview(paste(file,"cdt",sep ="."))
	if(saveDataObjects) { 
		gimmOut <- NULL
		gimmOut$clustData <- cdt
		gimmOut$hGClustData <- gtr
		gimmOut$hSClustData <- atr
		save(gimmOut, fClustAnnotations, file = paste(file,"RData",sep ="."))
	}
	invisible(list(data = cdt, rclust = gtr, cclust = atr, fClustAnnotations= fClustAnnotations))
}

