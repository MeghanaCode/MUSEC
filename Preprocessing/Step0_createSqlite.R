### create sqlite data from gtf file
gtfFile="/cfs/klemming/projects/supr/snic2020-6-4/Nghia/referenceDB/hg38.refGene.gtf.gz"
out="/cfs/klemming/projects/supr/snic2020-6-4/Nghia/referenceDB/hg38.refGene.sqlite"

args = commandArgs(trailingOnly=TRUE)
cat("\nNumber of arguments: ",length(args))
cat("\nList of arguments: ",args,"\n")

for (i in 1:length(args)){
  res=unlist(strsplit(args[i],"="))
  if (res[1]=="gtfFile") gtfFile=as.character(res[2])
  if (res[1]=="out") out=as.character(res[2])
}


require(GenomicFeatures)
gtfTxdb <- makeTxDbFromGFF(file=gtfFile,
                           format="gtf",
                           dataSource=paste("Link to the source",sep=""),
                           organism="Homo sapiens")
saveDb(gtfTxdb,file=out)
cat("\nSqlite file generated\n")
