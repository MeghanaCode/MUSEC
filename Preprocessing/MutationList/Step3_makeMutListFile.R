# Command
 Rscript ../preprocess/Step3_makeMutListFile.R \
  gtfFn="/nfs/NGS/Meghana/refs/hg38.refGene.gtf.gz" 
  gtfSqliteFn="/nfs/NGS/Meghana/refs/hg38.refGene.sqlite" 
  fastaFn="/nfs/NGS/Meghana/refs/refMrna.fa.gz" 
  exonicMutIn="mutations/Breast_exonicMut_B1.RData" 
  geneIn="/nfs/NGS/Meghana/refs/10x_Annotated_Human_Pan_Cancer_Panel.txt" 
  out="breast_mutation_list_batch_1.txt" 
  faOut="breast_tx_batch_1.fa"

### filter the mutation data to keep only mutations from the selected gene list
# process the cases when a mutations belonging to two genes

gtfFn="/nfs/NGS/Meghana/refs/hg38.refGene.gtf.gz"
gtfSqliteFn="/nfs/NGS/Meghana/refs/hg38.refGene.sqlite"
fastaFn="/nfs/NGS/Meghana/refs/refMrna.fa.gz"
exonicMutIn="Breast_exonicMut.RData"
geneIn="/nfs/NGS/Meghana/refs/10x_Annotated_Human_Pan_Cancer_Panel.txt"
out="breast_mutation_list.txt"
faOut="breast_tx.fa"

args = commandArgs(trailingOnly=TRUE)
cat("\nNumber of arguments: ",length(args))
cat("\nList of arguments: ",args,"\n")

for (i in 1:length(args)){
  res=unlist(strsplit(args[i],"="))
  if (res[1]=="gtfFn") gtfFn=as.character(res[2])
  if (res[1]=="gtfSqliteFn") gtfSqliteFn=as.character(res[2])
  if (res[1]=="fastaFn") fastaFn=as.character(res[2])
  if (res[1]=="exonicMutIn") exonicMutIn=as.character(res[2])
  if (res[1]=="geneIn") geneIn=as.character(res[2])
  if (res[1]=="out") out=as.character(res[2])
  if (res[1]=="faOut") faOut=as.character(res[2])
}

library(data.table)
require(GenomicFeatures)
require(Biostrings)
require(R.utils)

##load annotation data
anntxdb <- loadDb(gtfSqliteFn)
genes.all = genes(anntxdb, single.strand.genes.only = FALSE )
#tx.all = transcripts(anntxdb)
#exon_all=exons(anntxdb)
genes.tx.all = suppressMessages(suppressWarnings(select(anntxdb, keys=names(genes.all), columns=c("GENEID", "TXNAME","EXONID","EXONSTART","EXONEND","EXONCHROM","EXONSTRAND"), keytype = "GENEID")))
chrList=paste0("chr",c(1:22,"X","Y","M"))
genes.tx.all=genes.tx.all[genes.tx.all$EXONCHROM %in% chrList,]
allGenes=unique(genes.tx.all$GENEID )


load(exonicMutIn)
dim(exonicMut)
geneList=fread(geneIn, header=FALSE)

## if one mutation has a single gene, it is simple
pick1=exonicMut$gene %in% geneList$V2 #consider single gene
exonicMut1=exonicMut[pick1,]
dim(exonicMut1)

## if one mutation has more than one gene
d=strsplit(exonicMut$gene, ";")
p1=which(lengths(d)>1)
d1=d[p1]
d2=lapply(d1, function(x){
  pick=which(x %in% allGenes)
  x1=x[pick]
  if (sum(x1 %in% geneList$V2)>0) return(x1) else return(NULL)
})
p2=which(lengths(d2)>0)
d3=d2[p2]
exonicMut2=exonicMut[p1[p2],]
d4=sapply(d3, function(x){
  pick=which(x %in% geneList$V2)
  x1=x[pick[1]]#select the first one if more than 1 belonging to the gene panel
  return(x1)
})
exonicMut2$normGenes=exonicMut2$gene #genes from standard annotation
exonicMut2$gene=d4 #genes will be used to make mutation list
exonicMut1$normGenes=exonicMut1$gene
exonicMutFinal=rbind(exonicMut1,exonicMut2)


## now, we extract fasta sequence for those genes
g1=unique(exonicMut1$gene)
g2=unique(unlist(d3))
genesForFasta=unique(c(g1,g2)) #include both genes of interest and overlapping genes
#txForFasta=genes.tx.all$TXNAME[genes.tx.all$GENEID %in% genesForFasta & genes.tx.all$EXONCHROM %in% chrList]
txForFasta=genes.tx.all$TXNAME[genes.tx.all$GENEID %in% genesForFasta ] # genes.tx.all was filterd by chrList before
txForFasta=unique(txForFasta)

allFasta=readDNAStringSet(fastaFn)
allFasta_txname=sapply(names(allFasta),function(x){strsplit(x," ")[[1]][1]})
p=allFasta_txname%in% txForFasta
selectedFasta=allFasta[p]
writeXStringSet(selectedFasta,faOut)

## create gtf file for those genes. 
# NOTE: for exonicMut2, information of multiple genes will be assigned to only genes of interest
gtf=fread(gtfFn)
gname=sapply(gtf$V9,function(x){
    x1=strsplit(x,"gene\\_name")[[1]][2]
    x2=trimws(x1, whitespace = "[ \t\r\n\";]")
    return(x2)
  })
gname=strsplit(gtf$V9,"gene\\_name")
gname=setDT(gname)
gname=unlist(as.vector(as.data.frame(gname[2,])))
names(gname)=NULL
gname=trimws(gname, whitespace = "[ \t\r\n\";]")
pick=which(gname %in% genesForFasta)
gtf1=gtf[pick,]
gname1=gname[pick]

#now revise the information for multi-gene mutations
gtf2=gtf1
for (i in 1:nrow(exonicMut2)){
  g=unlist(strsplit(exonicMut2$normGenes[i],";"))
  g0=exonicMut2$gene[i]
  g1=setdiff(g,g0)
  p=which(gname1 %in% g1)
  mygtf=gtf2[p,]
  #replace gene names by the gene of interest
  for (gg in g1) mygtf$V9=gsub(gg,g0,mygtf$V9)
  gtf2$V9[p]=mygtf$V9
}

#remove data not from chrList
#chrList=paste0("chr",c(1:22,"X","Y","M"))
gtf1=gtf1[gtf1$V1 %in% chrList,]
gtf2=gtf2[gtf2$V1 %in% chrList,]
#export raw gtf of the genes of interest
fwrite(gtf1, file=paste0(faOut,"_raw.gtf"), sep="\t",row.names = FALSE, col.names = FALSE)
#export revised gtf
fwrite(gtf2, file=paste0(faOut,"_fixedNames.gtf"), sep="\t",row.names = FALSE, col.names = FALSE)


## make the mutation list file as the input for MAX
mutation_list=exonicMutFinal[,c("chr","start","end","ref_allele","mut_allele","gene","patientID")]
mutID=paste(exonicMutFinal$pos,exonicMutFinal$ref_allele,exonicMutFinal$mut_allele,sep="_")
x1=table(mutID)
x2=tapply(exonicMutFinal$patientID,mutID,function(x){
  paste(x,collapse=",")
})
exonicMutFinal$mutID=mutID
exonicMutFinal$Samples=x2[match(exonicMutFinal$mutID,names(x2))]
p=!duplicated(mutID)
mutation_list=exonicMutFinal[p,c("chr","start","end","ref_allele","mut_allele","gene","Samples")]
colnames(mutation_list)=c("Chromosome", "Start_pos", "End_pos", "REF", "ALT", "Gene",  "Samples")
write.table(mutation_list,file=out,row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

