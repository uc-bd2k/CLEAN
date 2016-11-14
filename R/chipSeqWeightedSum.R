`chipSeqWeightedSum`<-
function(ChipList,file=NULL,verbose=FALSE)
{
	require(GenomicFeatures)
	emExpUniform<-function(data,p=0.1,lamda=0.1,steps=1000,stopCond=0.01)
	{
		maxValue<-max(data)
		minValue<-min(data)
		likelihood1<-sum(log(p*lamda*exp(-lamda*data)+(1-p)/(maxValue-minValue)))
		likelihood2<-9999
		if(verbose)cat("Log likelihood ",likelihood1, "\nIteration\tLikelihood\tWeights\tLamda\n")
		iter<-1
		while(abs(likelihood2-likelihood1) > stopCond & iter < steps)
		{
			if(iter > 1) likelihood1 <-likelihood2
			tao<-p*lamda*exp(-lamda*data)/(p*lamda*exp(-lamda*data)+(1-p)/(maxValue-minValue))
			p<-sum(tao)/length(data)
			lamda<-sum(tao)/sum(tao*data)
			likelihood2<-sum(log(p*lamda*exp(-lamda*data)+(1-p)/(maxValue-minValue)))
			iter<-iter+1
			if(verbose)cat(iter,likelihood2,p,lamda,"\n")
		}
		res<-list(p=p,lamda=lamda)
		res
	} 

    emExpNormal<-function(data,p=0.1,lamda=0.1,mean=2.0,var=1.0,steps=1000,stopCond=0.01)
    {
	      likelihood1<-sum(log(p*lamda*exp(-lamda*data)+(1-p)*dnorm(data,mean=mean,sd=sqrt(var))))
	      likelihood2<-9999
	     if(verbose) cat("Log likelihood ",likelihood1, "\nIteration\tLikelihood\tWeights\tLamda\tMean\tVariance\n")
	      iter<-1
	      while(abs(likelihood2-likelihood1) > stopCond & iter < steps)
	      {
		      if(iter > 1) likelihood1 <-likelihood2
		      tao1<-p*lamda*exp(-lamda*data)/(p*lamda*exp(-lamda*data)+(1-p)*dnorm(data,mean=mean,sd=sqrt(var)))
		      tao2<-1-tao1
		      p<-sum(tao1)/length(data)
		      lamda<-sum(tao1)/sum(tao1*data)
		      mean<-sum(tao2*data)/sum(tao2)
		      var<-sum(tao2*(data-mean)^2)/sum(tao2)
		      likelihood2<-sum(log(p*lamda*exp(-lamda*data)+(1-p)*dnorm(data,mean=mean,sd=sqrt(var))))
		      iter<-iter+1
		      if(verbose)cat(iter,likelihood2,p,lamda,mean,var,"\n")
	      }
	      res<-list(p=p,lamda=lamda,mean=mean,var=var)
	      res
      } 

    emExpNormal.prior<-function(data,p=0.1,lamda=0.1,mean=2.0,var=1.0,steps=20,stopCond=0.01,background=3,probs=0.5)
    {
	    if(background <= 0) background <- 0.1
	    lamda<--log(probs)/background
	    mean<-background+1
	    likelihood1<-sum(log(p*lamda*exp(-lamda*data)+(1-p)*dnorm(data,mean=mean,sd=sqrt(var))))
	    likelihood2<-9999
	    if(verbose) cat("Log likelihood ",likelihood1, "\nIteration\tLikelihood\tWeights\tLamda\tMean\tVariance\n")
	    iter<-1
	    while(abs(likelihood2-likelihood1) > stopCond & iter < steps)
	    {
		    if(iter > 1) likelihood1 <-likelihood2
		    tao1<-p*lamda*exp(-lamda*data)/(p*lamda*exp(-lamda*data)+(1-p)*dnorm(data,mean=mean,sd=sqrt(var)))
		    tao2<-1-tao1
		    p<-sum(tao1)/length(data)
		    mean<-sum(tao2*data)/sum(tao2)
		    var<-sum(tao2*(data-mean)^2)/sum(tao2)
		    likelihood2<-sum(log(p*lamda*exp(-lamda*data)+(1-p)*dnorm(data,mean=mean,sd=sqrt(var))))
		    iter<-iter+1
		    if(verbose) cat(iter,likelihood2,p,lamda,mean,var,"\n")
	    }
	    res<-list(p=p,lamda=lamda,mean=mean,var=var)
	    res
    } 

    estimate.background<-function(data,para)
    {
	    minDist<-min(as.numeric(data[,"TSSDist"]))
	    maxDist<-max(as.numeric(data[,"TSSDist"]))
	    
	    if(is.element("EntrezID",colnames(data))) windowsize<-table(data[,"EntrezID"])
	    else windowsize<-table(data[,"GeneID"])
	    windowsize<-max(windowsize)
	    windowsize<-as.integer((maxDist-minDist)/windowsize)
	    res<-0
	    iter<-minDist
	    bin.start<-seq(minDist,maxDist-windowsize,by=windowsize)
	    bin.end<-seq(windowsize,maxDist,by=windowsize)
	    bins<-IRanges(start=bin.start,end=bin.end)
	    query<-IRanges(as.numeric(data[,"TSSDist"]),as.numeric(data[,"TSSDist"])+0.5)
	    tables<-matchMatrix(findOverlaps(query,bins))
	    bin<-unique(tables[,2])
	    peaks<-unlist(sapply(bin,function(x) mean(as.numeric(data[tables[tables[,2]==x,1],"Score"]))))
	    Weights<-para$p*para$lamda*exp(-para$lamda*(bin.start+windowsize/2))
	    Weights<-Weights/(Weights+(1-para$p)/(maxDist-minDist))
	    res<-Weights*peaks
	    res<-sum(res)
	    
	    res
    }

    WeightedSum<-function(data,TSSDist.Bound=10000)
    {
	    data<-data[as.numeric(data[,"TSSDist"])<=TSSDist.Bound,]
	    para<-emExpUniform(as.numeric(data[,"TSSDist"]),p=0.1,lamda=0.1)
	    if(verbose) cat("Estimate background...")
	    uniform<-estimate.background(data,para)
	    if(verbose) cat(log(uniform+1),"\n")
	    Weights<-para$p*para$lamda*exp(-para$lamda*as.numeric(data[,"TSSDist"]))
	    Weights<-Weights/(Weights+(1-para$p)/(max(as.numeric(data[,"TSSDist"]))-min(as.numeric(data[,"TSSDist"]))))
	    data<-cbind(data,WeightScore=Weights*as.numeric(data[,"Score"]))
	    Weighted.Sum<-matrix(0,nrow=length(unique(data[,"RefID"])),ncol=3)
	    colnames(Weighted.Sum)<-c("RefID","EntrezID","Score")
	    if(verbose) cat("Summarize Chip.Seq peaks...\n")
	    temp<-split(data,data[,"RefID"])
	    Weighted.Sum[,1]<-names(temp)
	    Weighted.Sum[,2]<-sapply(temp,function(x)x[1,"GeneID"])
	    Weighted.Sum[,3]<-sapply(temp,function(x) sum(x[,"WeightScore"]))
	    rownames(Weighted.Sum)<-names(temp)
	    Weighted.Sum<-data.frame(Weighted.Sum,stringsAsFactors=F)

	    res<-list(Score=Weighted.Sum,Uniform=uniform)
	    res
    }
  
    WeightedSum.Prob<-function(data,file=NULL,prob=0.005,prior=TRUE,log=TRUE)
    {
	data.2<-data[[1]]
	if(log==TRUE)
	{
		data.2[,"Score"]<-log(as.numeric(data.2[,"Score"])+1)
		backgrounds<-log(data[[2]]+1)
	}
	else backgrounds<-data[[2]]
	if(prior==TRUE) para<-emExpNormal.prior(as.numeric(data.2[,"Score"]),background=backgrounds,probs=prob)
	else para<-emExpNormal(as.numeric(data.2[,"Score"]))
	Prob<-(1-para$p)*dnorm(as.numeric(data.2[,"Score"]),mean=para$mean,sd=sqrt(para$var))
	Prob<-Prob/(Prob+para$p*para$lamda*exp(-para$lamda*as.numeric(data.2[,"Score"])))
 	data.2<-cbind(data.2,Prob=Prob)
	rownames(data.2)<-data.2[,1]
	data.2<-data.2[,-1]
	
	if(!is.null(file)) save(data.2,file=file)
	data.2
      }

    if(verbose) cat("Sum Chip-seq peaks...\n")
    Chip.Weighted<-WeightedSum(ChipList,TSSDist.Bound=1000000)
    if(verbose) cat("Get binding probs...\n")
    Chip.Weighted.Prob<-WeightedSum.Prob(Chip.Weighted,file=file,prob=0.001,prior=TRUE)

    Chip.Weighted.Prob
} 