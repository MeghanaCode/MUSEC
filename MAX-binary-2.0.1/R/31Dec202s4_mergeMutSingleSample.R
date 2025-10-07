### extract and merge mutant isoforms into mut_all
# Run for sample data: 
# - Rscript mergeMutSingleSample.R sampleEq=$inputFn sampleMut=$sampleMutFn sampleID=$sampleFn
# Run for X-matrix generation
# - Rscript mergeMutSingleSample.R sampleEq=$inputFn

### nghiavtr/05April2024: speed up with data.table

inputFn="sample_10/eqClass.txt"
sampleMutFn="MAX_mut_list_keepMutID.txt"
sampleFn="sample_10"
sampleMutFn="NotUse"
sampleFn="NotUse"

args = commandArgs(trailingOnly=TRUE)
cat("\nNumber of arguments: ",length(args))
cat("\nList of arguments: ",args,"\n")

for (i in 1:length(args)){
  res=unlist(strsplit(args[i],"="))
  if (res[1]=="sampleMut") sampleMutFn=res[2]
  if (res[1]=="sampleID") sampleFn=res[2]
  if (res[1]=="sampleEq") inputFn=res[2]
}

library(data.table)

source("/nfs/NGS/Meghana/MAX-binary-2.0.1/R/Rsource.R")

doc1=inputFn
cat("\n Processing ",doc1)
a = read.table(doc1,header=TRUE, stringsAsFactors=FALSE)

mut.list=read.csv(sampleMutFn,header=TRUE, sep="\t",stringsAsFactors = FALSE)
if (!("Samples" %in% colnames(mut.list))){
  cat("\n There is no column 'Samples' containing the information of sample list for each mutation in the mutation list file!!!")
  sampleMutFn=="NotUse"
}

condition=sampleFn=="NotUse" | sampleMutFn=="NotUse"

if (!condition){
  cat("\n Create final mutant isoforms relevant to the mutations of the sample.")

  #Nghia /05Feb2022:
  s=mut.list$Samples
  names(s)=mut.list$MutID
  #Nghia/29Nov2023: fixed the sMap names
  #sMap=sapply(s,function(x) sort(unique(trimws(strsplit(x,",")[[1]]))))
  sMap=lapply(s,function(x) sort(unique(trimws(strsplit(x,",")[[1]]))))
  sAll=unique(unlist(sMap))

  myMut=NULL
  for (i in 1:length(sMap)) if (sampleFn %in% sMap[[i]]){
    myMut=c(myMut,names(sMap)[i])
  }

  txAll=unique(a$Transcript)
  pick=grep("mut",txAll)
  txMut=txAll[pick]
  txWt=txAll[-pick]

  selectedTxMut=sapply(myMut,function(x) txMut[grep(x,txMut)])
  selectedTxMut=unique(unlist(c(selectedTxMut)))
  txSelected=c(txWt,selectedTxMut)
  #refine a, keep only relevant mutant isoforms
  pick=a$Transcript %in% txSelected
  a1=a[pick,]
  #need to merge the same eqclasses after reduction
  x=a1
  x_reduced=NULL
  for(i in unique(x$eqClass)){
    x1 = NULL
    x1 = subset(x,eqClass==i)
    x1 = x1[!duplicated(x1$Transcript),]
    x1 = x1[order(x1$Transcript),]
    x_reduced=rbind(x_reduced,x1)
  }
  x=x_reduced
  x_reduced=NULL
  x_merge = NULL
  for(i in unique(x$eqClass)){ 
    x1 = NULL
    x1 = subset(x,eqClass==i)
    x1 = x1[order(x1$Transcript),]
    x1 = cbind(x1,Isoform=paste(x1$Transcript,collapse='_'))
    x_merge = rbind(x_merge,x1)
  }
  z_merge = NULL
  allIsoforms=unique(x_merge$Isoform)
  for(i in 1:length(allIsoforms)){ 
    b = NULL
    b = subset(x_merge,Isoform==allIsoforms[i])[,-7]
    b1 = rbind(tapply(b$Weight,b$Transcript,sum),tapply(b$Count,b$Transcript,sum),tapply(b$EffLength,b$Transcript,mean),tapply(b$RefLength,b$Transcript,mean))
    b1 = t(b1)
    b2 = data.frame(Transcript=rownames(b1),b1,eqClass=i)
    colnames(b2) = colnames(b)
    z_merge = rbind(z_merge,b2)
  }
  rownames(z_merge) = c(1:nrow(z_merge))
  a=z_merge

#### nghiavtr/05April2024: speed up with data.table - does not work correctly, need to revise if needed
#  x_merge=data.table(a1)
#  x_merge=unique(x_merge,by=c("eqClass","Transcript"))
#  x_merge[,Isoform:=paste(Transcript,collapse='_'), by = c("eqClass")]
#
#  x_merge[,Weight2:=sum(Weight), by=c("Isoform","Transcript")]
#  x_merge[,Count2:=sum(Count), by=c("Isoform","Transcript")]
#  x_merge[,EffLength2:=mean(EffLength), by=c("Isoform","Transcript")]
#  x_merge[,RefLength2:=mean(RefLength), by=c("Isoform","Transcript")]
#
#  x_merge$Weight=x_merge$Weight2
#  x_merge$Count=x_merge$Count2
#  x_merge$EffLength=x_merge$EffLength2
#  x_merge$RefLength=x_merge$RefLength2
#  a=as.data.frame(x_merge[,c(1:6)])


}else{
  cat("\n Create final mutant isoforms from all mutations.")
}

doc2 = gsub('eqClass','Mutated_eqClass',doc1)
allTx=sort(unique(as.character(a$Transcript)))
tt=allTx[grep("mut",allTx,invert=TRUE)]
z_merge=mergeMut(a,tt)
write.table(z_merge,file=doc2,quote=F,row.names=F,sep='\t')