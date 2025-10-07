### estimation
# if mutation-sample mapping is available, example command:
# - Rscript estimateBeta.R workdir=Ycount sampleMut=mutDat.txt out=MAX_isoform_expression.RData
# if mutation-sample mapping is NOT available, example command:
# - Rscript estimateBeta.R workdir=Ycount out=MAX_isoform_expression.RData
# 06April2024/Nghia: improve codes to run with multiple genes

# default settings
sampleMutFn=NULL
workdir="Ycount_MAX2"
fout="MAX2_isoform_expression.RData"

maxiter.X=1
maxiter.b=100
modify=FALSE
keepEst=FALSE
core=4

args = commandArgs(trailingOnly=TRUE)
cat("\nNumber of arguments: ",length(args))
cat("\nList of arguments: ",args,"\n")

for (i in 1:length(args)){
  res=unlist(strsplit(args[i],"="))
  if (res[1]=="workdir") workdir=res[2]
  if (res[1]=="sampleMut") sampleMutFn=res[2]
  if (res[1]=="out") fout=res[2] 
  if (res[1]=="maxiterX") maxiter.X=as.integer(res[2])
  if (res[1]=="maxiterb") maxiter.b=as.integer(res[2])
  if (res[1]=="core") core=as.integer(res[2])
  if (res[1]=="modify") modify=as.logical(res[2])
  if (res[1]=="keepEst") keepEst=as.logical(res[2])
}

source("/cfs/klemming/projects/snic/snic2020-6-4/Meghana/MAX-binary-2.0.1/R/Rsource.R")

#set parallel
library(foreach)
library(doParallel)
ncores = detectCores()
nc = min(ncores,core)     # use 8 or 16 as needed!!
cl <- makePSOCKcluster(nc)   #
registerDoParallel(cl)

flist = list.files(workdir,pattern="_Ycount.RData",recursive=TRUE,full.names = TRUE)
estDataFull=list()
if (!is.null(sampleMutFn)){
  cat("\n Estimation of isoform expression (AEM)")
  ### using AEM
  flist2 = list.files(workdir,pattern="_Ycount.RData",recursive=TRUE,full.names = FALSE)
  sNameID=gsub("_Ycount.RData","",flist2)
  
  
  #Nghia /05Feb2022:
  mut.list=read.csv(sampleMutFn,header=TRUE, sep="\t",stringsAsFactors = FALSE)
  s=mut.list$Samples
  names(s)=mut.list$MutID
  sMap=sapply(s,function(x) sort(unique(trimws(strsplit(x,",")[[1]]))))
  sAll=unique(unlist(sMap))
  
  fileMutDat=NULL
  for (s in sAll){
    myMut=NULL
    for (i in 1:length(sMap)) if (s %in% sMap[[i]]){
      myMut=c(myMut,names(sMap)[i])
    }
    fileMutDat=rbind(fileMutDat,c(s,paste(sort(myMut), collapse=";")))
  }
  colnames(fileMutDat)=c("V1","V2")
  #  fileMutDat=read.csv(sampleMutFn,header=FALSE, sep="\t",stringsAsFactors = FALSE)
  sGroup=tapply(fileMutDat[,1],fileMutDat[,2],c)
} else {
  #consider single sample at once, so will use EM instead of AEM
  flist2 = list.files(workdir,pattern="_Ycount.RData",recursive=TRUE,full.names = FALSE)
  sNameID=gsub("_Ycount.RData","",flist2)
  sGroup=lapply(sNameID, function(x) x)
  names(sGroup)=sNameID
}
  
  maxList=list()
  txAll=NULL
  for (g in 1:length(sGroup)){
    pick=which(sNameID %in% sGroup[[g]])
    flist3=flist[pick]
    gSnames=sNameID[pick]
    
    ## merging Y count
    X.y=NULL
    sample1m=NULL
    for(id in 1:length(flist3)){
      cat("Merging results from sample ",flist3[id],' ...\n')
      load(flist3[id])
      sample1m=c(sample1m,samplename1)
      if(id==1){
        X.y = Y
      }
      if(id>1){
        for(i in 1:length(Y))
          X.y[[i]] = cbind(X.y[[i]],sample1=Y[[i]][,'sample1'])
      }
    }
    
    
    fun = function(crpdat, maxiter.X=5, maxiter.b=100, modify=TRUE){
      xloc = which(colnames(crpdat) != "sample1")
      X0 = matrix(crpdat[,xloc], ncol=length(xloc))
      Ymat = crpdat[,-xloc, drop=FALSE]
      est = AEM(X0, Ymat, maxiter.X=maxiter.X,maxiter.b=maxiter.b, modify=modify)
      return(est)
    }

    EST <- foreach(i= 1:length(X.y)) %dopar% fun(X.y[[i]], maxiter.X=maxiter.X,maxiter.b=maxiter.b, modify=modify)
    
    names(EST) = names(X.y)
    
    th2=NULL
    for(i in 1:length(X.y))
    {
      x.y =  X.y[[i]]
      xloc = which(colnames(x.y) != "sample1")
      X0 = x.y[,xloc]
      Ymat = x.y[,-xloc,drop=FALSE]
      estRes=EST[[i]]
      
      th1=estRes$BETA
      colnames(th1)=colnames(X0)
      rownames(th1)=sample1m
      th2=cbind(th2,th1)
    }
    

    # fn=flist3[1]
    # load(fn)
    # crpdat=Y[[1]]
    # xloc = which(colnames(crpdat) != "sample1")
    # X0 = crpdat[,xloc]
    # txAll=unique(c(txAll,colnames(X0)))
    # Ymat=NULL
    # for (fn in flist3){
    #   load(fn)
    #   crpdat=Y[[1]] #problem here, we use only the first one
    #   xloc = which(colnames(crpdat) != "sample1")
    #   X0 = crpdat[,xloc]
    #   Ymat = cbind(Ymat,crpdat[,-xloc,drop=FALSE])
    # }
    # colnames(Ymat)=gSnames
    # 
    # estRes = AEM(X0, Ymat, maxiter.X=maxiter.X, maxiter.b=maxiter.b, modify=modify)
    #  convergeVec=c(convergeVec,estRes$beta.conv)
    # th2=estRes$BETA
    # colnames(th2)=colnames(X0)
    # rownames(th2)=gSnames
    
    maxList[[g]]=th2
    
    ### save by group
    #estDataFull[[names(sGroup)[g]]]=list(X0=X0,Ymat=Ymat)
    estDataFull[[names(sGroup)[g]]]=list(X.y=X.y,sample1m=sample1m,EST=EST)    
    
  }
  
  #create the final expression matrixt
  txAll=lapply(maxList, colnames)
  txAll=sort(unique(unlist(txAll)))
  maxMat=matrix(0,length(txAll),length(sNameID))
  dim(maxMat)
  rownames(maxMat)=txAll
  colnames(maxMat)=sNameID
  for (g in 1:length(maxList)){
    thMat=maxList[[g]]
    for (i in 1:nrow(thMat)){
      x=thMat[i,drop=FALSE,]
      maxMat[colnames(x),rownames(x)]=x
    }
  }

# }else{ #consider single sample at once, so will use EM instead of AEM
#   ### Using EM
#   cat("\n Estimation of isoform expression (EM)")
#   maxList=list()
#   txAll=NULL
#   convergeVec=NULL
#   
#   
#   for (fn in flist){
#     load(fn)
#     crpdat=Y[[1]]
#     xloc = which(colnames(crpdat) != "sample1")
#     X0 = crpdat[,xloc]
#     Ymat = crpdat[,-xloc,drop=FALSE]
#     estRes = AEM(X0, Ymat, maxiter.X=1, maxiter.b=maxiter.b, modify=modify)
#     convergeVec=c(convergeVec,estRes$beta.conv)
#     th2=estRes$BETA[1,]
#     names(th2)=colnames(X0)
#     txAll=unique(c(txAll,colnames(X0)))
#     maxList[[samplename1]]=th2
#     estDataFull[[fn]]=list(X0=X0,Ymat=Ymat)
#   }
#   maxMat=matrix(0,length(txAll),length(maxList))
#   dim(maxMat)
#   colnames(maxMat)=names(maxList)
#   rownames(maxMat)=txAll
#   for (i in 1:ncol(maxMat)){
#     x=maxList[[colnames(maxMat)[i]]]
#     maxMat[names(x),i]=x
#   }
# }

isoform_count=maxMat
save(isoform_count, file=fout)
if (keepEst) save(estDataFull, file=paste0("estData",fout))
