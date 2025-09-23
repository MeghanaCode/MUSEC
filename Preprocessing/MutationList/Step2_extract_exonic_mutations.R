### filter the mutation data to keep only EXONIC mutations
# Rscript Step2_extract_exonic_mutations.R gtfSqliteFn="/nfs/NGS/Meghana/refs/hg38.refGene.sqlite" mutDataIn="breast_Mutect2.RData" out="Breast_exonicMut.RData"


gtfSqliteFn="/nfs/NGS/Meghana/refs/hg38.refGene.sqlite"
mutDataIn="breast_Mutect2.RData"
out="Breast_exonicMut.RData"

args = commandArgs(trailingOnly=TRUE)
cat("\nNumber of arguments: ",length(args))
cat("\nList of arguments: ",args,"\n")

for (i in 1:length(args)){
  res=unlist(strsplit(args[i],"="))
  if (res[1]=="gtfSqliteFn") gtfSqliteFn=as.character(res[2])
  if (res[1]=="mutDataIn") mutDataIn=as.character(res[2])
  if (res[1]=="out") out=as.character(res[2])
}

library(data.table)
require(GenomicFeatures)

#load annotation data
anntxdb <- loadDb(gtfSqliteFn)
genes.all = genes(anntxdb, single.strand.genes.only = FALSE )
#tx.all = transcripts(anntxdb)
#exon_all=exons(anntxdb)
genes.tx.all = suppressMessages(suppressWarnings(select(anntxdb, keys=names(genes.all), columns=c("GENEID", "TXNAME","EXONID","EXONSTART","EXONEND","EXONCHROM","EXONSTRAND"), keytype = "GENEID")))

#load mutation data
load(mutDataIn)
length(mutData)

exonicMut=NULL
#now do filtering for individual patients
patient_num=1
while (patient_num <=length(mutData$MUTGENE)){
  cat("\n processing patient th",patient_num)
  mutgene_col=as.vector(unlist(mutData[["MUTGENE"]][[patient_num]]))
  ref_allele_col=as.vector(unlist(mutData[["Reference_allele"]][[patient_num]]))
  mut_allele_col=as.vector(unlist(mutData[["Mutated_allele"]][[patient_num]]))
  mut_pos_col=as.vector(unlist(mutData[["mutPos"]][[patient_num]]))
  
  patient_1=data.frame(pos=mut_pos_col,ref_allele=ref_allele_col,mut_allele=mut_allele_col,gene=mutgene_col)
  #create mut
  #pos
  patient_1$pos=as.character(patient_1$pos)
  #chr
  patient_1$chr=sapply(patient_1$pos,function(x){strsplit(x,split = "_")[[1]][1]})
  
  patient_1$chr=substr(patient_1$chr,4,nchar(patient_1$chr))
  #start
  patient_1$start=sapply(patient_1$pos,function(x){as.numeric(strsplit(x,split = "_")[[1]][2])})
  #refer_allele
  patient_1$ref_allele=as.character(patient_1$ref_allele)
  #mut_allele
  patient_1$mut_allele=as.character(patient_1$mut_allele)
  #end
  patient_1$end=patient_1$start

  myend=apply(patient_1, 1, function(x){
    x1=as.numeric(x[6])+max(nchar(x[2]),nchar(x[3]))-1
    return(x1)
  })

  patient_1$end=myend


  patient_1$allgenes=patient_1$gene
  gene.list=lapply(patient_1$gene,function(x){unique(unlist(strsplit(x,split=";")))})
  gene.list=lapply(gene.list,function(x){x=x[x!=""]})
  patient_1$gene=sapply(gene.list,function(x)paste(x,collapse=";"))


m1 = GRanges(paste0("chr",patient_1$chr), IRanges(start = as.numeric(patient_1$start), end = as.numeric(patient_1$end)),  id=patient_1$pos)

anntxdb_range=GRanges(genes.tx.all$EXONCHROM, IRanges(start = as.numeric(genes.tx.all$EXONSTART), end = as.numeric(genes.tx.all$EXONEND)), id=genes.tx.all$geneID)

overlap = findOverlaps(m1,anntxdb_range, type = "within")
qh=queryHits(overlap)
sh=subjectHits(overlap)
res=data.frame(pos=m1$id[qh],GENEID=genes.tx.all$GENEID[sh])
res$id=paste(res$pos,res$GENEID, sep="__")
res=res[!duplicated(res$id),]

dupID=tapply(res$GENEID,res$pos,function(x) paste(x,collapse=";"))

mut.list=patient_1[patient_1$pos %in% names(dupID),]
mut.list$gene=dupID[match(mut.list$pos,names(dupID))]

#  for(i in 1:nrow(patient_1)){
#    g=gene.list[[i]]
#    mygene=genes.tx.all[genes.tx.all$GENEID %in% g,] ##mygene strand info
#    #mygene_txlist=unique(mygene$TXNAME)
#
#    p=(patient_1$start[i] >= mygene$EXONSTART) & (patient_1$end[i] <= mygene$EXONEND)
#    if (sum(p)>0){
#      exonicGenes=unique(mygene[p,]$GENEID)
#      patient_1$gene[i]=paste(exonicGenes,collapse=";")
#    }else{
#      patient_1$gene[i]="NotExonic"
#    }
#  }
  
  if (nrow(mut.list)>0){
    mut.list$patientID=names(mutData$mutPos)[patient_num]
    exonicMut=rbind(exonicMut,mut.list)
  }else{
    cat("\n no exonic mutations in ", patient_num)
  }

  patient_num=patient_num+1
}

save(exonicMut, file=out)




