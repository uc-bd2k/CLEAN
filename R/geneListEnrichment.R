`geneListEnrichment` <- 
function(geneList,
						allGenes,
						functionalCategories = NULL, 
						species = NULL,
						minGenesInCategory=10,
						maxGenesInCategory=1000,
						inBkg=TRUE,
						sigFDR = 0.1,
						verbose=TRUE) {
	if(is.null(functionalCategories)) {
		functionalCategories <- "GO" #default is GO categories
	}
	if(is.character(functionalCategories)) {
		if(substr(functionalCategories,nchar(functionalCategories)-4,nchar(functionalCategories)) == "RData") {
			e1 <- new.env()
			load(file = functionalCategories, envir = e1)
			if(length(ls(envir = e1, pattern = "2[G,g][E,e][N,n][E,e]")) > 0) {
				l <- get(ls(envir = e1, pattern = "2[G,g][E,e][N,n][E,e]")[1], envir = e1)
				#names(l) <- names(functionalCategories)
				#if(is.null(names(l))) names(l) <- substr(functionalCategories, 1, nchar(functionalCategories)-6)
				CategoryID2GeneID <- l
			} else {
				CategoryID2GeneID <- NA
				warning(paste("No functional categories found in file", functionalCategories, "-  No functional clustering annotation generated."))
			}
			if(length(ls(envir = e1, pattern = "2[D,d][E,e][S,s][C,c]")) > 0) {
				l <- get(ls(envir = e1, pattern = "2[D,d][E,e][S,s][C,c]")[1], envir = e1)
				#names(l) <- names(functionalCategories)
				#if(is.null(names(l))) names(l) <- substr(functionalCategories, 1, nchar(functionalCategories)-6)
				CategoryID2Desc <- l
			} 
			else CategoryID2Desc <- NA
		} else {
			functionalCategories <- getFunctionalCategories(functionalCategories, species)
			CategoryID2GeneID <- functionalCategories[[1]][[1]]
			if(length(functionalCategories[[1]]) > 1) 
				CategoryID2Desc <- functionalCategories[[1]][[2]]
			else CategoryID2Desc <- NA
		}
	} else {   #functionalCategories is a list
		l <- functionalCategories[[1]]
		CategoryID2GeneID <- l
		if(length(functionalCategories) > 1) 
			CategoryID2Desc <- functionalCategories[[2]]
		else CategoryID2Desc <- NA
	}
	categories <- as.data.frame(geneSigInList(geneList=geneList, allGenes=allGenes, CategoryID2GeneID=CategoryID2GeneID,  
		minGenesInCategory=minGenesInCategory, maxGenesInCategory=maxGenesInCategory, verbose=verbose, inBkg=inBkg, sigFDR = sigFDR),stringsAsFactors=F)
	data.frame(categories[,1], Description = CategoryID2Desc[match(categories[,1], CategoryID2Desc[,1]),2], categories[,-1],stringsAsFactors=F)
}

