`genesInEnrichedCategories` <- 
function(categoryIDs, geneList, funcCategories = NULL, species = NULL) {
	if(is.null(funcCategories)) funcCategories <- "GO"  #default is GO categories
	if(is.character(funcCategories)) {
		if(substr(funcCategories,nchar(funcCategories)-4,nchar(funcCategories)) == "RData") {
			e1 <- new.env()
			load(file = funcCategories, envir = e1)
			if(length(ls(envir = e1, pattern = "2[G,g][E,e][N,n][E,e]")) > 0) {
				l <- get(ls(envir = e1, pattern = "2[G,g][E,e][N,n][E,e]")[1], envir = e1)
				#names(l) <- names(funcCategories)
				#if(is.null(names(l))) names(l) <- substr(funcCategories, 1, nchar(funcCategories)-6)
				CategoryID2GeneID <- l
			} else {
				CategoryID2GeneID <- NA
				warning(paste("No functional categories found in file", funcCategories, "-  No functional clustering annotation generated."))
			}
			if(length(ls(envir = e1, pattern = "2[D,d][E,e][S,s][C,c]")) > 0) {
				l <- get(ls(envir = e1, pattern = "2[D,d][E,e][S,s][C,c]")[1], envir = e1)
				#names(l) <- names(funcCategories)
				#if(is.null(names(l))) names(l) <- substr(funcCategories, 1, nchar(funcCategories)-6)
				CategoryID2Desc <- l
			} else 
				CategoryID2Desc <- list(NA)
		} else {
			l <- getFunctionalCategories(funcCategories, species)
			if(length(l) > 0) {
				CategoryID2GeneID <- l[[1]][[1]]
				if(length(l[[1]]) > 1) {
					CategoryID2Desc <- l[[1]][[2]]
				} else {
					CategoryID2Desc <- list(NA)
				}
			} else {
				CategoryID2GeneID <- NA
				CategoryID2Desc <- NA
				warning(paste("functional categories", funcCategories, "not found.  No functional clustering annotation generated."))
			}
		}
	} else {   #funcCategories is a list
		l <- funcCategories[[1]]
		CategoryID2GeneID <- l
		if(length(funcCategories) > 1)
			CategoryID2Desc <- funcCategories[[2]]
		else CategoryID2Desc <- list(NA)
	}
	genesInCategories <- sapply(categoryIDs, function(id) intersect(CategoryID2GeneID[[as.character(id)]], geneList))
	enrichedGenes <- matrix(rep(FALSE, length(geneList)*length(categoryIDs)), length(categoryIDs))
	for(i in 1:length(genesInCategories)) {
		enrichedGenes[i,match(genesInCategories[[i]], geneList)] <- TRUE
	}
	colnames(enrichedGenes) <- geneList
	data.frame(CategoryID2Desc[match(categoryIDs, CategoryID2Desc[,1]),], enrichedGenes,stringsAsFactors = F)
}

