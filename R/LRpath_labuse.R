LrPath<-function (sigvals, geneids, min.g = 10, max.g = NA, sig.cutoff = 0.05, functionalCategories = NULL, odds.min.max = c(0.001, 0.5), species = "Hs") ## this function is modified by Lixia Zhang based on maureen's LRpath function
{
  if (!(is.factor(geneids) | is.numeric(geneids) | is.character(geneids))) {
    warning("Parameter `geneids` should be a character or numeric vector, or a factor.")
  }
  if (is.factor(geneids)) {
    geneids <- as.character(geneids)
    warning("Converted geneids from factor to character.")
  }
  xx <- functionalCategories
  ENTREZ.in.DB <- unique(unlist(xx))
  if (length(intersect(ENTREZ.in.DB, unique(geneids))) < min.g) 
    warning(paste("min.g=", min.g, " is larger than number of gene IDs overlapping with ", 
                  names(CategoryID2GeneID)[ii], ". Please consider lower min.g.", 
                  sep = ""))
  sigvals[sigvals == 0] <- 10^(-15)
  geneids2 <- geneids[!is.na(geneids) & !is.na(sigvals) ]
  sigvals2 <- sigvals[!is.na(geneids) & !is.na(sigvals) ]
  uniqids <- unique(geneids2)
  numuniq <- length(uniqids)
  LOR.mult <- log(odds.min.max[2]) - log(odds.min.max[1])
  nlp <- (-1) * log(sigvals2)
  newp <- NA
  for (i in (1:numuniq)) {
    current <- nlp[geneids2 == uniqids[i]]
    newp[i] <- mean(current)
  }
  catsizes <- sapply(xx, length)
  yy <- xx[catsizes >= min.g]
  ncats <- length(yy)
  siggenes <- uniqids[exp(-newp) < sig.cutoff]
  LRcoeff <- NA
  LRpval <- NA
  cattots <- NA
  catsigIDs <- NA
  ind <- 0
  if (is.na(max.g)) {
    max.g <- 99999
  }
  for (i in (1:ncats)) {
    catgenes <- as.character(yy[[i]])
    catgenes <- unique(catgenes)
    catpoprows <- match(catgenes, uniqids)
    catpoprows <- catpoprows[!is.na(catpoprows)]
    cattot <- length(catpoprows)
    if (cattot >= min.g & cattot <= max.g) {
      ind <- ind + 1
      cat <- rep(0, numuniq)
      cat[catpoprows] <- 1
      forLR <- as.data.frame(cbind(cat, newp))
      names(forLR) <- c("cat", "nlogp")
      withCallingHandlers(glm.lrGO <- glm(cat ~ 
                                            nlogp, family = binomial(link = "logit"), 
                                          forLR), warning = function(w) {
                                            if (length(grep("fitted probabilities numerically 0 or 1 occurred", 
                                                            as.character(w))) == 1) 
                                              warning("fitted probabilities numerically 0 or 1 occurred: functional category appears to induce linear separation of p-values.")
                                          })
      lrGO <- summary(glm.lrGO)
      LRcoeff[ind] <- lrGO$coefficients["nlogp", 
                                        "Estimate"]
      LRpval[ind] <- lrGO$coefficients["nlogp", 
                                       "Pr(>|z|)"]
      cattots[ind] <- cattot
      catsig.IDlist <- intersect(siggenes, catgenes)
      catsigIDs[ind] <- paste(catsig.IDlist[order(catsig.IDlist)], 
                              collapse = ", ")
    }
    if (i%%100 == 0) {
      cat(".")
    }
  }
  cat("\n")
  odds.ratio <- exp(LOR.mult * LRcoeff)
  BenjFDR <- p.adjust(LRpval, "BH")
  LRres <- cbind(cattots, 
                 LRcoeff, odds.ratio, as.data.frame(LRpval), 
                 BenjFDR, catsigIDs)
  names(LRres) <- c(
    "n.genes", "coeff", "odds.ratio", "p.value", 
    "FDR", "sig.genes")
  LRresults <-  list(LRres)
  names(LRresults)[length(LRresults)] <- names(functionalCategories)
  LRresults <- LRresults[[1]]
  invisible(LRresults)
}

