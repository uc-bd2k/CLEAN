'enrichmentGSEA'<-  function(sigvals=NULL,geneids,diffexp=NULL, functionalCategories=c("KEGG","GO","tfacts","dorotheatfs","hallmark"), directions=c("up","down","both","combined"),species="Hs", excel=F,sigFDR=0.1,nperm=100000,minSize = 15, maxSize = 500){
  require(CLEAN)
  print(paste0("CLEAN.",species),character.only=T)
  require(paste0("CLEAN.",species),character.only=T)
  require(fgsea)

  if(is.list(functionalCategories)) gs<-functionalCategories
  else gs<-getFunctionalCategories(functionalCategories,species)
  outlist=as.list(NULL)
# print(paste0("is.null(sigvals)=",is.null(sigvals)))
# print(paste0("is.null(diffexp)=",is.null(diffexp)))
  if("both" %in% directions){
    if(is.null(diffexp)){
      stats<-sigvals
      names(stats)<-geneids
      stats=sort(stats)
    }else if(is.null(sigvals)){
      stats<-diffexp
      names(stats)<-geneids
      stats<-sort(stats)
    }else {
      stats<-ifelse(diffexp>0,-log10(sigvals),log10(sigvals))
      names(stats)<-geneids
      stats=sort(stats)
    }
    fgseaRes <- lapply(gs, function(x){
      if(is.null(nperm)) {
        fgseaOut<-fgsea(pathways = x[[1]],stats = stats, minSize=minSize, maxSize=maxSize,eps=0)
      } else{
        fgseaOut<-fgsea(pathways = x[[1]],stats = stats, minSize=minSize, maxSize=maxSize, nperm=nperm)
      }
      fgseaOut<-data.frame(ID=fgseaOut$pathway,Name=x[[2]][match(fgseaOut$pathway,x[[2]][,1]),2],fgseaOut[,-1],stringsAsFactors = F)
      fgseaOut[order(fgseaOut$pval),]
    })
    names(fgseaRes)<-names(gs)
    outlist=c(outlist,list(both=fgseaRes))
  }

  if("up" %in% directions){
    stats=-log10(ifelse(diffexp>0,sigvals/2,1-sigvals/2))
    names(stats)<-geneids
    stats=sort(stats)
    fgseaRes <- lapply(gs, function(x){
      if(is.null(nperm)) {
        fgseaOut<-fgsea(pathways = x[[1]],stats = stats, minSize=minSize, maxSize=maxSize)
      } else{
        fgseaOut<-fgsea(pathways = x[[1]],stats = stats, minSize=minSize, maxSize=maxSize, nperm=nperm)
      }
      fgseaOut<-data.frame(ID=fgseaOut$pathway,Name=x[[2]][match(fgseaOut$pathway,x[[2]][,1]),2],fgseaOut[,-1],stringsAsFactors = F)
      fgseaOut[order(fgseaOut$pval),]
    })
    names(fgseaRes)<-names(gs)
    outlist=c(outlist,list(up=fgseaRes))
  }

  if("down" %in% directions){
    stats=-log10(ifelse(diffexp<0,sigvals/2,1-sigvals/2))
    names(stats)<-geneids
    stats=sort(stats)
    fgseaRes <- lapply(gs, function(x){
      if(is.null(nperm)) {
        fgseaOut<-fgsea(pathways = x[[1]],stats = stats, minSize=minSize, maxSize=maxSize)
      } else{
        fgseaOut<-fgsea(pathways = x[[1]],stats = stats, minSize=minSize, maxSize=maxSize, nperm=nperm)
      }
      fgseaOut<-data.frame(ID=fgseaOut$pathway,Name=x[[2]][match(fgseaOut$pathway,x[[2]][,1]),2],fgseaOut[,-1],stringsAsFactors = F)
      fgseaOut[order(fgseaOut$pval),]
    })
    names(fgseaRes)<-names(gs)
    outlist=c(outlist,list(down=fgseaRes))
  }

  if("combined" %in% directions){
    statsUp=-log10(ifelse(diffexp>0,sigvals/2,1-sigvals/2))
    names(statsUp)<-geneids
    statsUp=sort(statsUp)
    print(head(statsUp))
    statsDown=-log10(ifelse(diffexp<0,sigvals/2,1-sigvals/2))
    names(statsDown)<-geneids
    statsDown=sort(statsDown)
    print(head(statsDown))
    fgseaRes <- lapply(gs, function(x){
      if(is.null(nperm)) {
        fgseaOutUp<-fgsea(pathways = x[[1]],stats = statsUp, minSize=minSize, maxSize=maxSize)
          print(head(fgseaOutUp))
        fgseaOutDown<-fgsea(pathways = x[[1]],stats = statsDown, minSize=minSize, maxSize=maxSize)
        print(head(fgseaOutDown))
      } else{
        fgseaOutUp<-fgsea(pathways = x[[1]],stats = statsUp, minSize=minSize, maxSize=maxSize, nperm=nperm)
        fgseaOutDown<-fgsea(pathways = x[[1]],stats = statsDown, minSize=minSize, maxSize=maxSize, nperm=nperm)
      }
      print("fgseaOutUp")
      print(fgseaOutUp)
      print("fgseaOutDown")
      print(fgseaOutDown)
      fgseaOutUp<-fgseaOutUp[match(fgseaOutDown$pathway,fgseaOutUp$pathway),]
      posGreater<-fgseaOutUp$ES>abs(fgseaOutDown$ES)
      fgseaOut<-data.frame(ID=fgseaOutDown$pathway,Name=x[[2]][match(fgseaOutDown$pathway,x[[2]][,1]),2],fgseaOutDown[,-1],stringsAsFactors = F)
      fgseaOut$ES<-(-fgseaOut$ES)
      fgseaOut$NES<-(-fgseaOut$NES)
      fgseaOut[posGreater,-c(1:2)]<-data.frame(fgseaOutUp[posGreater,-1],stringsAsFactors = F)
      fgseaOut[order(fgseaOut$pval),]
    })
    names(fgseaRes)<-names(gs)
    outlist=c(outlist,list(combined=fgseaRes))
  }

## fgsea calls not working
# if("upnew" %in% directions){
#   # stats=-log10(ifelse(diffexp>0,sigvals/2,1-sigvals/2))
#   stats<-ifelse(diffexp>0,-log10(sigvals),log10(sigvals))
#   names(stats)<-geneids
#   stats=sort(stats)
#   fgseaRes <- lapply(gs, function(x){
#     if(is.null(nperm)) {
#       fgseaOut<-fgsea(pathways = x[[1]],stats = stats, minSize=minSize, maxSize=maxSize,scoreType="pos")
#     } else{
#       fgseaOut<-fgsea(pathways = x[[1]],stats = stats, minSize=minSize, maxSize=maxSize, nperm=nperm,scoreType="pos")
#     }
#     print(head(fgseaOut))
#     fgseaOut<-data.frame(ID=fgseaOut$pathway,Name=x[[2]][match(fgseaOut$pathway,x[[2]][,1]),2],fgseaOut[,-1],stringsAsFactors = F)
#     fgseaOut[order(fgseaOut$pval),]
#   })
#   names(fgseaRes)<-names(gs)
#   outlist=c(outlist,list(up=fgseaRes))
# }
#
# if("downnew" %in% directions){
#   # stats=-log10(ifelse(diffexp<0,sigvals/2,1-sigvals/2))
#   stats<-ifelse(diffexp>0,-log10(sigvals),log10(sigvals))
#   names(stats)<-geneids
#   stats=sort(stats)
#   fgseaRes <- lapply(gs, function(x){
#     if(is.null(nperm)) {
#       fgseaOut<-fgsea(pathways = x[[1]],stats = stats, minSize=minSize, maxSize=maxSize,scoreType="neg")
#     } else{
#       fgseaOut<-fgsea(pathways = x[[1]],stats = stats, minSize=minSize, maxSize=maxSize, nperm=nperm,scoreType="neg")
#     }
#     fgseaOut<-data.frame(ID=fgseaOut$pathway,Name=x[[2]][match(fgseaOut$pathway,x[[2]][,1]),2],fgseaOut[,-1],stringsAsFactors = F)
#     fgseaOut[order(fgseaOut$pval),]
#   })
#   print(head(fgseaOut))
#   names(fgseaRes)<-names(gs)
#   outlist=c(outlist,list(down=fgseaRes))
# }



  if(excel){
    require(xlsx)
    upexcel <- createWorkbook()
    for (geneset in genesets){
      assign(paste0("up",geneset),createSheet(wb=upexcel,sheetName=paste0("up",geneset)))
    }
    for (geneset in genesets){
      # addDataFrame(x=get(paste0("rsup$",geneset)),sheet=paste0("up",geneset))
      addDataFrame(x=rsup[geneset][[1]],sheet=get(paste0("up",geneset)))
    }
    outlist=c(outlist,list(upexcel=upexcel))
    downexcel <- createWorkbook()
    for (geneset in genesets){
      assign(paste0("down",geneset),createSheet(wb=downexcel,sheetName=paste0("down",geneset)))
    }
    for (geneset in genesets){
      # addDataFrame(x=get(paste0("rsup$",geneset)),sheet=paste0("up",geneset))
      addDataFrame(x=rsdown[geneset][[1]],sheet=get(paste0("down",geneset)))
    }
    outlist=c(outlist,list(downexcel=downexcel))
    bothexcel <- createWorkbook()
    for (geneset in genesets){
      assign(paste0("both",geneset),createSheet(wb=bothexcel,sheetName=paste0("both",geneset)))
    }
    for (geneset in genesets){
      # addDataFrame(x=get(paste0("rsup$",geneset)),sheet=paste0("up",geneset))
      addDataFrame(x=rsboth[geneset][[1]],sheet=get(paste0("both",geneset)))
    }
    outlist=c(outlist,list(bothexcel=bothexcel))
  }
  return(outlist)
}
