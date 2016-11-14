GRS <- function(query.p, reference.p, query.gL=NULL, reference.gL=NULL, na.rm=TRUE, 
	estimateNullDistr=TRUE, nullDistrQuantiles=c(0.9, 0.95, 0.99), nullDistrN=100, 
	tolerateWarnings=TRUE, pAdjust.method.query=NULL, pAdjust.method.reference=NULL, lambda=0.005) {

	pAdjust <- function (p, method="top", pAdjustN=100, sigLevel=0.1) {
		if(method == "top" & !is.null(pAdjustN)) {
			index <- rank(p, ties.method="random") > pAdjustN
			p[index] <- 1
			p[!index] <- sigLevel
			p
		} else {
			p.adjust(p, method)
		}
	}
	postP <- function(p, scale=1) {
		B <- -exp(1) * p * log(p)
		B[p > exp(-1)] <- 1
		1-1/(1+1/(max(1, scale)*B))
	}
	estimateH0 <- function(p, lambda) {
		pi0 <- 1
		tryCatch(pi0 <- qvalue(p, lambda=lambda)$pi0, error = function(e) {pi0 <- 1})
		pi0 <- min(pi0, 0.999)
		pi0/(1-pi0)
	}
	RSCM <- function(p1, p2, s1, s2) {

		E <- 1/2 * (sum(p1*s2)/sum(p1) + sum(p2*s1)/sum(p2))

		mu.p1 <- mean(p1)
		mu.p2 <- mean(p2)
		mu.s1 <- mean(s1)
		mu.s2 <- mean(s2)
		sigma.p1 <- sd(p1)
		sigma.p2 <- sd(p2)
		sigma.s1 <- sd(s1)
		sigma.s2 <- sd(s2)
		gamma1 <- sum((p1 - mu.p1)*(s1 - mu.s1))/(length(s1) - 1)
		gamma2 <- sum((p2 - mu.p2)*(s2 - mu.s2))/(length(s2) - 1)

		mu.E <- (mu.s1 + mu.s2) / 2
	
		delta <- matrix(c(1/mu.p1, -(mu.s2/mu.p1), 1/mu.p2, -(mu.s1/mu.p2)), 4, 1)
		Sigma.11 <- sigma.p1^2 * sigma.s2^2 + sigma.p1^2 * mu.s2^2 + sigma.s2^2 * mu.p1^2
		Sigma.12 <- mu.s2 * sigma.p1^2
		Sigma.13 <- gamma1 * gamma2 + gamma1 * mu.p2 * mu.s2 + gamma2 * mu.p1 * mu.s1
		Sigma.14 <- gamma2 * mu.p1
		Sigma.22 <- sigma.p1^2
		Sigma.23 <- gamma1 * mu.p2
		Sigma.24 <- 0
		#Sigma.24 <- var(p1, p2)
		Sigma.33 <- sigma.p2^2 * sigma.s1^2 + sigma.p2^2 * mu.s1^2 + sigma.s1^2 * mu.p2^2
		Sigma.34 <- mu.s1 * sigma.p2^2
		Sigma.44 <- sigma.p2^2
		Sigma <- matrix(c(Sigma.11, Sigma.12, Sigma.13, Sigma.14, 
		Sigma.12, Sigma.22, Sigma.23, Sigma.24, 
		Sigma.13, Sigma.23, Sigma.33, Sigma.34, 
		Sigma.14, Sigma.24, Sigma.34, Sigma.44), 4, 4)
		sigma.E <- 0.5 * sqrt(((t(delta) %*% Sigma) %*% delta)/length(p1)) 
		list(p.value=(pnorm((E-mu.E)/sigma.E, lower.tail=FALSE)), 
			z.score=(E-mu.E)/sigma.E, 
			E.gene=length(p1)/2 * ((p1*s2)/sum(p1) + (p2*s1)/sum(p2)))
	}
		
	if(!is.null(dim(query.p)) | !is.null(query.gL)) {
		if(!is.null(dim(query.p))) {
			query.gL <- query.p[,1]
			query.p <- query.p[,2:dim(query.p)[2]]
			if(is.null(dim(reference.p))) { 
				if(tolerateWarnings) warning("parameter reference.p should be a matrix or data.frame")
				else return("parameter reference.p should be a matrix or data.frame")
			} else {
				reference.gL <- reference.p[,1]
				reference.p <- reference.p[,2:dim(reference.p)[2]]
			}
		}

		pValue.reference<-FALSE
		pValue.query<-FALSE
		if(is.null(dim(reference.p))) pValue.reference<-TRUE
		if(is.null(dim(query.p))) pValue.query<-TRUE
		if(is.null(query.gL)) { 
			if(tolerateWarnings) warning("parameter query.gL should not be NULL")
			else return("parameter query.gL should not be NULL")
		} else {
			if(pValue.query) NAindex <- is.na(query.p)
			else NAindex <- is.na(query.p[,1])
			if(any(NAindex)) {
				if(na.rm) {
					query.gL <- query.gL[!NAindex]
					if(pValue.query) query.p <- query.p[!NAindex]
					else query.p <- query.p[!NAindex,]
				} else {
					if(tolerateWarnings) warning("query.p has NA's but na.rm is FALSE")
					else return("query.p has NA's but na.rm is FALSE")
				}
			}
			uniqueGenes <- unique(query.gL)
			if(length(uniqueGenes) != length(query.gL)) {
			      if(pValue.query){
				averagePvalues <- split(query.p, query.gL)
				names(averagePvalues) <- as.character(names(averagePvalues))
				averagePvalues <- unlist(lapply(averagePvalues, function(x) exp(mean(log(x), na.rm=TRUE))))
				query.gL <- uniqueGenes
				query.p <- averagePvalues[match(query.gL, names(averagePvalues))]
			      }
			      else
			      {
				averageScore <- split(query.p[,1], query.gL)
				averageProb <- split(query.p[,2], query.gL)
				names(averageScore) <- as.character(names(averageScore))
				names(averageProb) <- as.character(names(averageProb))
				averageScore <- unlist(lapply(averageScore, function(x) mean(as.numeric(x), na.rm = TRUE)))
				averageProb <- unlist(lapply(averageProb, function(x) mean(as.numeric(x), na.rm = TRUE)))
				query.gL <- uniqueGenes
				query.p <- averageScore[match(query.gL, names(averageScore))]
				query.p <- cbind(query.p,averageProb[match(query.gL, names(averageProb))])
				colnames(query.p)<-c("Score","Prob")
			      }
			}
		}
		if(is.null(reference.gL)) { 
			if(tolerateWarnings) warning("parameter reference.gL should not be NULL")
			else("parameter reference.gL should not be NULL")
		} else {
			if(pValue.reference) NAindex <- is.na(reference.p)
			else NAindex <- is.na(reference.p[,1])
			
			if(any(NAindex)) {
				if(na.rm) {
					reference.gL <- reference.gL[!NAindex]
					if(pValue.reference) reference.p <- reference.p[!NAindex]
					else reference.p <- reference.p[!NAindex,]
				} else {
					if(tolerateWarnings) warning("reference.p has NA's but na.rm is FALSE")
					else return("reference.p has NA's but na.rm is FALSE")
				}
			}
			uniqueGenes <- unique(reference.gL)
			if(length(uniqueGenes) != length(reference.gL)) {
				if(pValue.reference){
				  averagePvalues <- split(reference.p, reference.gL)
				  names(averagePvalues) <- as.character(names(averagePvalues))
				  averagePvalues <- unlist(lapply(averagePvalues, function(x) exp(mean(log(x), na.rm=TRUE))))
				  reference.gL <- uniqueGenes
				  reference.p <- averagePvalues[match(reference.gL, names(averagePvalues))]
				}
				else
				{
				  averageScore <- split(reference.p[,1], reference.gL)
				  averageProb <- split(reference.p[,2], reference.gL)
				  names(averageScore) <- as.character(names(averageScore))
				  names(averageProb) <- as.character(names(averageProb))
				  averageScore <- unlist(lapply(averageScore, function(x) mean(as.numeric(x), na.rm = TRUE)))
				  averageProb <- unlist(lapply(averageProb, function(x) mean(as.numeric(x), na.rm = TRUE)))
				  reference.gL <- uniqueGenes
				  reference.p <- averageScore[match(reference.gL, names(averageScore))]
				  reference.p <- cbind(reference.p,averageProb[match(reference.gL, names(averageProb))])
				  colnames(reference.p)<-c("Score","Prob")
				}
			}
		}
		commonGenes <- intersect(query.gL, reference.gL)
		if(length(commonGenes) < 2)  {
			if(tolerateWarnings) warning("number of common gene IDs is < 2")
			else return("number of common gene IDs is < 2")
		}
		
		if(pValue.query) query.p <- query.p[match(commonGenes, query.gL)]
		else query.p <- query.p[match(commonGenes, query.gL),]
		if(pValue.reference) reference.p <- reference.p[match(commonGenes, reference.gL)]
		else reference.p <- reference.p[match(commonGenes, reference.gL),]
	}
	
	if(pValue.query) s1 <- -log(query.p, 10)
	else  s1 <- query.p[,"Score"]
	if(pValue.reference) s2 <- -log(reference.p, 10)
	else s2 <- reference.p[,"Score"]

	if(is.null(pAdjust.method.query)) {
		if(pValue.query){
		  query.m0 <- estimateH0(query.p, lambda)
		  p1 <- postP(query.p, scale=query.m0)
		  query.m0 <- query.m0/(query.m0 + 1)
		}
		else  
		{
		  p1 <- query.p[,"Prob"]
		  query.m0 <- estimateH0(1-p1, lambda)
		  query.m0 <- query.m0/(query.m0 + 1)
		}
	} else {
		p1 <- 1 - p.adjust(query.p, method=pAdjust.method.query)
		query.m0 <- NULL
	}
	if(is.null(pAdjust.method.reference)) {
		if(pValue.reference){
		  reference.m0 <- estimateH0(reference.p, lambda)
		  p2 <- postP(reference.p, scale=reference.m0)
		  reference.m0 <- reference.m0 / (reference.m0 + 1)
		}
		else
		{
		      p2 <- reference.p[,"Prob"]
		      reference.m0 <- estimateH0(1-p2, lambda)
		      reference.m0 <- reference.m0/(reference.m0 + 1)
		}
	} else {
		p2 <- 1 - p.adjust(reference.p, method=pAdjust.method.reference)
		reference.m0 <- NULL
	}

	res <- RSCM(p1, p2, s1, s2)
	geneTable <- data.frame(geneID=commonGenes, E.value=res$E.gene)
	ret <- list(p.value=res$p.value, z.score=res$z.score, geneTable=geneTable[order(res$E.gene, decreasing=TRUE), ], 
					query.m0=query.m0, reference.m0=reference.m0)
	if(estimateNullDistr) {
		f <- function(p1, p2, s1, s2) {
			index1 <- sample(length(p1))
			index2 <- sample(length(p2))
			RSCM(p1[index1], p2[index2], s1[index1], s2[index2])$E.gene
		}
		geneTable<-geneTable[order(res$E.gene, decreasing=TRUE), ]
		geneTable<-cbind(geneTable,Escore.pValue=rep(0,dim(geneTable)[1]),Escore.qValue=rep(0,dim(geneTable)[1]))
		temp <- sapply(1:nullDistrN, function(i) list(f(p1, p2, s1, s2)))
		distr<-sapply(temp,function(x)quantile(x,probs=nullDistrQuantiles))
		temp<-unlist(temp)
		temp<-cbind(temp,rep(1,length(temp)))
    		junk<-rbind(as.matrix(geneTable[,c(2,3)]),temp)
    		junk<-junk[order(junk[,1],decreasing=T),]
   		counts<-which(junk[,2]==0)
    		for(i in 1:length(counts)) counts[i]<-counts[i]-i
        	geneTable[,3]<-counts/(nullDistrN*length(counts))
		geneTable[which(geneTable[,3]==0),3]<-1/(dim(temp)[1]+1)
    		geneTable[,4]<-qvalue(geneTable[,3])$qvalues

		if(length(nullDistrQuantiles) == 1) distr <- t(distr)
		distr <- rowMeans(distr)
		
		ret$geneTable<-geneTable
		ret <- c(ret, list(EvalueNullDistrQ=distr))
	}
	ret

}

