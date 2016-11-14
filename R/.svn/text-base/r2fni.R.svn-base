`r2fni` <-
function(hr, file = "cluster.fni", clustCategories = NULL, dec = ".") {
	get.genes <- function(merge,component) {
		#recursive version deprecated
		#if(component > 0) {
		#	c(get.genes(merge,merge[component,2]), get.genes(merge,merge[component,1]))
		#}  
		#else -component
	    toDoList <- component   
		doneList <- NULL
		while(length(toDoList)>0) {
			if(toDoList[1] > 0) {
				toDoList <- c(merge[toDoList[1],], toDoList[-1])
			} else {
				doneList <- c(-toDoList[1], doneList)
				toDoList <- toDoList[-1]
			}
		}
		doneList
	}

	CategoryCols <- c("nGenesInCategory","nAllGenesInCategory","FisherFDR", "logOR")
	fni <- NULL
	if(!is.null(clustCategories)) {
		fni <- rbind(fni, c("NodeID","CategoryType", "CategoryID", "CategoryDescription", CategoryCols))
		emptyCols <- rep("", length(CategoryCols)+ 2)
		i <- 1
		if(is.null(names(clustCategories))) {
			warning("Names of functional category not specified.  Using default")
			names(clustCategories) <- paste("CategoryType", 1:length(clustCategories))
		}
		while(i <= dim(hr$merge)[1]) {
			NodeID <- paste("NODE", i, "X", sep = "")
			genes <- get.genes(hr$merge, i)
			index <- which(unlist(lapply(clustCategories[[1]]$allClusters, length)) == length(genes))
			if(length(index) > 0) {
				cluster <- index[sapply(index, function(j) 
					setequal(clustCategories[[1]]$allClusters[[j]], genes))]
				if(length(cluster) > 0) {
					sigCategories <- NULL
					for(j in 1:length(clustCategories)) {
						type <- names(clustCategories)[j]
						if(!is.null(clustCategories[[j]]$sigCategories[[cluster]])) {
							if(dim(clustCategories[[j]]$sigCategories[[cluster]])[1] > 0) {
								ID <- as.character(clustCategories[[j]]$sigCategories[[cluster]][,1])
								#if(is.null(clustCategories[[j]]$CategoryID2Desc) | is.na(clustCategories[[j]]$CategoryID2Desc)) {
								if(length(clustCategories[[j]]$CategoryID2Desc) < 2) {
									ID <- cbind(ID, Description = "")
								} else {
									ID <- cbind(ID, Description = clustCategories[[j]]$CategoryID2Desc[match(ID, clustCategories[[j]]$CategoryID2Desc[,1]), 2])
								}
								sigCategories <- rbind(sigCategories, cbind(NodeID, type, ID, clustCategories[[j]]$sigCategories[[cluster]][,CategoryCols]))
							}
						}
					}	
					if(is.null(sigCategories)) {
						fni <- rbind(fni, c(NodeID, "cluster is not significantly enriched for any of the considered categories", emptyCols))
					} else {
						sigCategories <- sigCategories[order(sigCategories[, "FisherFDR"]),]
						sigCategories[,"NodeID"] <- as.character(sigCategories[,"NodeID"])
						sigCategories[,"type"] <- as.character(sigCategories[,"type"])
						sigCategories[,"ID"] <- as.character(sigCategories[,"ID"])
						sigCategories[,"Description"] <- as.character(sigCategories[,"Description"])
						sigCategories[,"nGenesInCategory"] <- as.character(sigCategories[,"nGenesInCategory"])
						sigCategories[,"nAllGenesInCategory"] <- as.character(sigCategories[,"nAllGenesInCategory"])
						sigCategories[,"FisherFDR"] <- format(sigCategories[,"FisherFDR"], digits = 3, scientific = TRUE)
						sigCategories[,"logOR"] <- format(sigCategories[,"logOR"], digits = 3)
						fni <- rbind(fni, as.matrix(sigCategories))
					}
				} else {
					fni <- rbind(fni, c(NodeID, "cluster was not found in functional annotation file", emptyCols))
				}
			} else {
				fni <- rbind(fni, c(NodeID, "cluster size is outside the range considered for functional annotation", emptyCols))
			}
			i <- i+1
		}
		write.table(fni, file = file, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
	}
}

