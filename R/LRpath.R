`LRpath` <-function(sigvals, geneids, min.g=10, max.g=NA, sig.cutoff=0.05,
functionalCategories=NULL, odds.min.max=c(0.001,0.5), species="Hs") {

## prepare database/functionalCategories, etc.
if(!(is.factor(geneids) | is.numeric(geneids) | is.character(geneids))) {
	warning("Parameter `geneids` should be a character or numeric vector, or a factor.")
}
if(is.factor(geneids)) {
	geneids <- as.character(geneids)
	warning("Converted geneids from factor to character.")
}

if (is.null(functionalCategories)) {
  warning("functionalCategories is NULL, use `GO` as default.")
	functionalCategories <- "GO"
}
if(!is.list(functionalCategories)) functionalCategories <- as.list(unique(functionalCategories))\
CategoryID2GeneID <- list()
CategoryID2Desc <- list()
runLRpath <- rep(TRUE, length(functionalCategories))
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
				runLRpath[i] <- FALSE
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
				runLRpath[i] <- FALSE
			    CategoryID2GeneID <- c(CategoryID2GeneID, list(NA))
				CategoryID2Desc <- c(CategoryID2Desc, list(NA))
				warning(paste("functional categories", functionalCategories[[i]], "not found.  No functional clustering annotation generated."))
			}
		}
    } else { #functionalCategories is a list
        l1 <- functionalCategories[i]
        if(length(l1) > 0) {
            l <- l1[[1]][1]
            names(l) <- names(functionalCategories[i])
            CategoryID2GeneID <- c(CategoryID2GeneID, l)
            if(length(l1[[1]]) > 1) {
                l <- l1[[1]][2]
                names(l) <- names(functionalCategories[i])
                CategoryID2Desc <- c(CategoryID2Desc, l)
            } else CategoryID2Desc <- c(CategoryID2Desc, list(NA))
        } else {
            runLRpath[i] <- FALSE
            CategoryID2GeneID <- c(CategoryID2GeneID, list(NA))
            CategoryID2Desc <- c(CategoryID2Desc, list(NA))
            warning(paste("functional categories", functionalCategories[[i]], "not found.  No functional clustering annotation generated."))
        }
    }
## run LRpath
if(any(runLRpath)) {
	LRresults <- list()
	for(ii in 1:length(runLRpath)) {
		if(runLRpath[ii]) {
			
			cat(paste("Analyzing", names(CategoryID2GeneID)[ii], "categories"))

			xx <- CategoryID2GeneID[[ii]]
			#xx <- lapply(xx, unique)
			ENTREZ.in.DB <- unique(unlist(xx))
			if(length(intersect(ENTREZ.in.DB, unique(geneids))) < min.g) 
				warning(paste("min.g=", min.g, " is larger than number of gene IDs overlapping with ", names(CategoryID2GeneID)[ii], ". Please consider lower min.g.", sep=""))

			sigvals[sigvals==0]<- 10^(-15)
			geneids2<-geneids[!is.na(geneids)&!is.na(sigvals)&geneids %in% ENTREZ.in.DB]
			sigvals2<-sigvals[!is.na(geneids)&!is.na(sigvals)&geneids %in% ENTREZ.in.DB]

			uniqids<-unique(geneids2)
			numuniq<-length(uniqids)
			LOR.mult<-log(odds.min.max[2])-log(odds.min.max[1])

			################ Average -log(p-values) for duplicate locusids
			nlp<-(-1)*log(sigvals2)
			newp<-NA
			for (i in (1:numuniq)){
				current<-nlp[geneids2==uniqids[i]]  #Get all values for i-th geneid
				#numrep<-length(current)
				newp[i]<-mean(current)
			}
			#plot(c(1,2,3),c(1,2,3))
			### Limit testing to those categories with at least min.g genes TOTAL
			catsizes<-sapply(xx,length)
			yy<-xx[catsizes>=min.g]
			ncats<-length(yy)

			#siggenes <- uniqids[exp(-newp)<0.05]  ##  If averaged p-value < 0.05
			siggenes <- uniqids[exp(-newp)<sig.cutoff]

			######################################################################
			###  Test all categories to which at least min ANALYZED genes belong

			LRcoeff<-NA; LRpval<-NA; cattots<-NA; yyind<-NA; GOnums<-NA
			catsigIDs<-NA
			ind<-0
			if (is.na(max.g)) {
				max.g<-99999
			}
			for (i in (1:ncats)){
				catgenes<-as.character(yy[[i]])
				catgenes<-unique(catgenes)
				catpoprows<-match(catgenes, uniqids)
				catpoprows<-catpoprows[!is.na(catpoprows)]
					cattot<-length(catpoprows)	#Number of analyzed genes in category
				if (cattot>=min.g&cattot<=max.g) {
					ind<-ind+1
					cat<-rep(0,numuniq)
					cat[catpoprows]<- 1
					forLR<- as.data.frame(cbind(cat,newp))
						names(forLR)<-c("cat","nlogp")

					withCallingHandlers(
						glm.lrGO <- glm(cat ~ nlogp, family=binomial(link="logit"),forLR),
						warning = function(w) {
							if(length(grep("fitted probabilities numerically 0 or 1 occurred", as.character(w)))==1)
							    warning("fitted probabilities numerically 0 or 1 occurred: functional category appears to induce linear separation of p-values.")})
					lrGO<-summary(glm.lrGO)
					###########  Extract p-value and coefficient
					LRcoeff[ind]<-lrGO$coefficients["nlogp","Estimate"]
					LRpval[ind]<-lrGO$coefficients["nlogp","Pr(>|z|)"]
					cattots[ind]<-cattot
					#yyind[ind]<-i
					GOnums[ind]<-names(yy[i])

					catsig.IDlist<-intersect(siggenes,catgenes)
					catsigIDs[ind]<-paste(catsig.IDlist[order(catsig.IDlist)],collapse=", ")
   				}
				if (i%%100 == 0) {
					cat(".")
				}
			}
			cat("\n")
			odds.ratio<-exp(LOR.mult*LRcoeff)
			BenjFDR<-p.adjust(LRpval,"BH")
			kterms <- CategoryID2Desc[[ii]][,2]
			k.ids <- CategoryID2Desc[[ii]][,1]
			keggrows<- match(GOnums,k.ids)
			kegg.annot<-kterms[keggrows]
			if(all(is.na(GOnums))) {
				warning(paste("No results found for ", names(CategoryID2GeneID)[ii], ".", sep=""))
				LRres <- NA
			} else {
				LRres <-cbind(GOnums,kegg.annot,cattots,LRcoeff, odds.ratio,as.data.frame(LRpval),BenjFDR,catsigIDs)
				names(LRres)<- c("category.ID", "category.description", "n.genes", "coeff", "odds.ratio", "p.value", "FDR", "sig.genes")
			}
			LRresults <- c(LRresults, list(LRres))
			names(LRresults)[length(LRresults)] <- names(CategoryID2GeneID)[ii]
		}
	}
}
}
LRresults <- LRresults[[1]]
invisible(LRresults)
}#end LRpath function