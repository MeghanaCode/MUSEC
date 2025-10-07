## Take the arguments
input="WT_Mut_Ref/WT_Mut_reference.fa"
k=10

args = commandArgs(trailingOnly=TRUE)
cat("\nNumber of arguments: ",length(args))
cat("\nList of arguments: ",args,"\n")

for (i in 1:length(args)){
        res=unlist(strsplit(args[i],"="))
        if (res[1]=="in") input=as.character(res[2])
        if (res[1]=="k") k=as.integer(res[2])

}

suppressMessages(suppressWarnings(require(Biostrings)))

myfasta=readDNAStringSet(input)
makeFolds<-function(k=n, n, seedVal=123456){
  # to do LOO-CV, just set k equal to n
  set.seed(seedVal)
  s=trunc(n/k)
  foldId=rep(c(1:k),s+1)[1:n]
  #sampleIDx=sample(seq(n))
  sampleIDx=seq(n) #not random anymore
  kFolds=split(sampleIDx,foldId)
  return(kFolds)
}
n=length(myfasta)
myfolds=makeFolds(k,n)

for (i in 1:k){
        dir.create(paste0("fold_",i))
        foldFa=myfasta[myfolds[[i]]]
        writeXStringSet(foldFa, paste0("fold_",i,"/subbatch.fa"))
}
