`RandomSet` <-
function(sigvals, geneids, 
    database="GO", functionalCategories=NULL, species="Hs", 
	min.g=10, minGenesInCategory=NULL, 
	max.g=NA, maxGenesInCategory=NULL, 
	sig.cutoff=0.1, sigFDR=NULL,
	averageMultipleProbes=TRUE,
	allGenesInCategoriesAsBackground=TRUE,
	two.sided=FALSE,
	na.rm=TRUE, verbose=TRUE) {

## converting sigvals to -log10 scale if necessary
if(any(is.na(sigvals))) {
	if(na.rm) {
		geneids <- geneids[!is.na(sigvals)]
		sigvals <- sigvals[!is.na(sigvals)]
	} else {
		warning("NA or NaN detected, but na.rm=FALSE. Check results.")
	}
}
if(max(sigvals, na.rm=TRUE) > 1) {
	warning("sigvals > 1 detected.  Assuming that sigvals are on -log10 scale.")
} else {
	sigvals <- -log(sigvals, 10)
	if(verbose) cat("converted sigvals to -log10 scale.\n")
}

## sanity checks etc.
if(!(is.factor(geneids) | is.numeric(geneids) | is.character(geneids))) {
	warning("Parameter `geneids` should be a character or numeric vector, or a factor.")
}
if(is.factor(geneids)) {
	geneids <- as.character(geneids)
	warning("Converted geneids from factor to character.")
}

if(!missing(database)) {
	warning("Argument `database` deprecated.  Please use `functionalCategories` instead.")
	if(missing(functionalCategories)) functionalCategories <- database
}
if (is.null(functionalCategories)) functionalCategories <- "GO"

if(!(missing(min.g) & missing(max.g))) {
	warning("Argument `min.g`, `max.g` deprecated.  Please use `minGenesInCategory`, `maxGenesInCategory` instead.")
	if(missing(minGenesInCategory)) minGenesInCategory <- min.g
	if(missing(maxGenesInCategory)) maxGenesInCategory <- max.g
}
if(is.null(minGenesInCategory)) {
	minGenesInCategory <- 0
	if(verbose) cat("minGenesInCategory set 0.\n")
}
if(is.null(maxGenesInCategory)) maxGenesInCategory <- NA
if(is.na(maxGenesInCategory)) {
	maxGenesInCategory <- Inf
	if(verbose) cat("maxGenesInCategory set Inf.\n")
}

if(!missing(sig.cutoff)) {
	warning("Argument `sig.cutoff` deprecated.  Please use `sigFDR` instead.")
	if(missing(sigFDR)) sigFDR <- sig.cutoff
}
if(is.null(sigFDR)) {
	sigFDR <- 1
	if(verbose) cat("sigFDR set to 1.\n")
}

## prepare database/functionalCategories, etc.
if(!is.list(functionalCategories)) functionalCategories <- as.list(unique(functionalCategories))
CategoryID2GeneID <- list()
CategoryID2Desc <- list()
runRS <- rep(TRUE, length(functionalCategories))
for(i in 1:length(functionalCategories)) {
    if (functionalCategories %in% unique(c(names(getFunctionalCategories(CLEAN.Hs(), species = "Hs")),
                                           names(getFunctionalCategories(CLEAN.Mm(), species = "Mm")),
                                           names(getFunctionalCategories(CLEAN.Rn(), species = "Rn"))))) {
		if(substr(functionalCategories[[i]],nchar(functionalCategories[[i]])-4,nchar(functionalCategories[[i]])) == "RData") {
			e1 <- new.env()
			load(file = functionalCategories[[i]], envir = e1)
			if(length(ls(envir = e1, pattern = "2[G,g][E,e][N,n][E,e]")) > 0) {
				l <- list(get(ls(envir = e1, pattern = "2[G,g][E,e][N,n][E,e]")[1], envir = e1))
				names(l) <- names(functionalCategories)[i]
				if(is.null(names(l))) names(l) <- substr(functionalCategories[[i]], 1, nchar(functionalCategories[[i]])-6)
				CategoryID2GeneID <- c(CategoryID2GeneID, l)
			} else {
				runRS[i] <- FALSE
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
			species <- paste(toupper(substr(species, 1, 1)), substr(species, 2, 2), sep = "")
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
				runRS[i] <- FALSE
				CategoryID2GeneID <- c(CategoryID2GeneID, list(NA))
				CategoryID2Desc <- c(CategoryID2Desc, list(NA))
				warning(paste("functional categories", functionalCategories[[i]], "not found.  No functional clustering annotation generated."))
			}
		}
	} else {   
	    l1 <- get(functionalCategories[[i]])
	        l <- l1[[1]][1]
	        names(l) <- functionalCategories[[i]]
	        CategoryID2GeneID <- c(CategoryID2GeneID, l)
	        if(length(l1[[1]]) > 1) {
	            l <- l1[[1]][2]
	            names(l) <- functionalCategories[[i]]
	            CategoryID2Desc <- c(CategoryID2Desc, l)
	        } else CategoryID2Desc <- c(CategoryID2Desc, list(NA))
	}
}

## RandomSet function (Newton et al.)
myAllez <- function(scores, geneList, categoryList) {
	categoryList <- unique(categoryList)
	index <- na.omit(match(categoryList, geneList))
	if(length(index) > 0) {
		m <- length(index)
		G <- length(scores)
		xBar <- mean(scores[index], na.rm=TRUE)
		mu <- mean(scores, na.rm=TRUE)
		sigma <- sqrt((sum(scores^2, na.rm=TRUE)/G-mu^2) * (G - m) / (G - 1) / m)
		#list(xBar=xBar, mu=mu, sigma=sigma, z=(xBar-mu)/sigma)
		c(n.genes=length(index), zScore=(xBar-mu)/sigma)
	} else {
		#list(NA)
		c(n.genes=0, zScore=NA)
	}
}

## run RandomSet method for each set of functional categories
if(any(runRS)) {
	
	#deal with multiple probes for the same geneID 
	if(averageMultipleProbes & length(geneids) > length(unique(geneids))) {
		if(verbose) cat("Averaging multiple probes per gene ID where applicable.\n")
		sigvals <- split(sigvals, geneids)
		geneids <- names(sigvals)
		sigvals <- unlist(lapply(sigvals, mean, na.rm=na.rm))
	}

	RSresults <- list()
	for(ii in 1:length(runRS)) {
		if(runRS[ii]) {
			if(verbose) cat(paste("Analyzing", names(CategoryID2GeneID)[ii], "categories.\n"))
			
			#remove categories that are smaller than minGenesInCategory and bigger than maxGenesInCategory
			categories <- CategoryID2GeneID[[ii]]
			categoryLengths <- unlist(lapply(categories, function(x) length(unique(x))))
			categories <- categories[categoryLengths <= maxGenesInCategory & categoryLengths >= minGenesInCategory]
			if(any(is.na(CategoryID2Desc[[ii]]))) { 
				categoryDescr <- NA
			} else {
				categoryDescr <- CategoryID2Desc[[ii]][match(names(categories), CategoryID2Desc[[ii]][,1]), 2]
			}
			if(allGenesInCategoriesAsBackground) {
				if(verbose) cat(paste("  .. Restricting gene list to", names(CategoryID2GeneID)[ii], "background gene list.\n"))
				geneids2 <- intersect(geneids, unique(unlist(categories)))
				sigvals2 <- sigvals[match(geneids2, geneids)]
			} else {
				geneids2 <- geneids
				sigvals2 <- sigvals
			}
			if(verbose) cat(paste("  .. Analzying", length(categories), "categories.\n"))
			res <- matrix(unlist(lapply(categories, function(categoryList) myAllez(sigvals2, geneids2, categoryList))), byrow=TRUE, ncol=2)
			if(two.sided) {
				pvalues <- pnorm(-abs(res[,2])) * 2
			} else {
				if(verbose) cat("  .. Computing one-sided p-values.\n")
				pvalues <- pnorm(-res[,2])
			}
			res <- data.frame(categoryID=names(categories), 
				description=categoryDescr, 
				nGenes=res[,1],
				zScore=res[,2], 
				pValue=pvalues, 
				FDR=p.adjust(pvalues, method="fdr"))
			res <- res[order(res[, "pValue"]), ]
			res <- res[which(res[,"FDR"] <= sigFDR), ]

			RSresults <- c(RSresults, list(res))
			names(RSresults)[length(RSresults)] <- names(CategoryID2GeneID)[ii]
		}
	}
}
invisible(RSresults)
}

