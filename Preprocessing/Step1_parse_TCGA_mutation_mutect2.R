### parse the mutation output from mutect2
# command: Rscript Step1_parse_TCGA_mutation_mutect2.R manifestFile=gdc_manifest_mutectbreast_filtered.txt mutPath=/cfs/klemming/projects/supr/snic2020-6-4/nobackup/AML/TCGA-LAML/AML_testset out=AMLTestSet_Mutect2.RData
#
# paramter setting

manifestFile="/nfs/NGS/Meghana/breastcancer/gdc_manifest/batches/mutect_batch_1.txt"
mutPath="/nfs/NGS/Meghana/breastcancer/mutations/mutect"
out="/nfs/NGS/Meghana/breastcancer/mutations/mutect/breast_Mutect2_batch1.RData"

#ml nano
#ml PDC/23.12 # to load R
#ml R/4.4.1-cpeGNU-23.12

args = commandArgs(trailingOnly=TRUE)
cat("\nNumber of arguments: ",length(args))
cat("\nList of arguments: ",args,"\n")

for (i in 1:length(args)){
  res=unlist(strsplit(args[i],"="))
  if (res[1]=="manifestFile") manifestFile=as.character(res[2])
  if (res[1]=="mutPath") mutPath=as.character(res[2])
  if (res[1]=="out") out=as.character(res[2])
}

library(TCGAutils)

#All fields in INFO column (8)
field_names=c("Allele","Consequence","IMPACT","SYMBOL","Gene","Feature_type","Feature","BIOTYPE","EXON","INTRON","HGVSc","HGVSp","cDNA_position","CDS_position","Protein_position","Amino_acids","Codons","Existing_variation","ALLELE_NUM","DISTANCE","STRAND","FLAGS","VARIANT_CLASS","SYMBOL_SOURCE","HGNC_ID","CANONICAL","TSL","APPRIS","CCDS","ENSP","SWISSPROT","TREMBL","UNIPARC","RefSeq","GENE_PHENO","SIFT","PolyPhen","DOMAINS","HGVS_OFFSET","GMAF","AFR_MAF","AMR_MAF","EAS_MAF","EUR_MAF","SAS_MAF","AA_MAF","EA_MAF","ExAC_MAF","ExAC_Adj_MAF","ExAC_AFR_MAF","ExAC_AMR_MAF","ExAC_EAS_MAF","ExAC_FIN_MAF","ExAC_NFE_MAF","ExAC_OTH_MAF","ExAC_SAS_MAF","CLIN_SIG","SOMATIC","PHENO","PUBMED","MOTIF_NAME","MOTIF_POS","HIGH_INF_POS","MOTIF_SCORE_CHANGE","ENTREZ","EVIDENCE")


###### results from mutect2
cat("\n Collect Mutect2.. \n")
manifestInfo=read.table(manifestFile,sep="\t",header=TRUE)


folderPath=mutPath
folderList=list.files(folderPath, pattern="*.gz$", recursive = TRUE)

patientID=MUTGENE=mutPos=Reference_allele=Mutated_allele=list()
for (i in 1:length(folderList)){
  cat(", ",i)
  matchID=match(unlist(strsplit(folderList[i],"/"))[2],as.character(manifestInfo$filename))
  #patientBarcode=as.character(manifestInfo$id[matchID])
  #use the function from TCGAutils to get patient ID
  patientBarcode=UUIDtoBarcode(manifestInfo$id[matchID], from_type = c("file_id"))
  patientBarcode=patientBarcode[1,2] #keep the TCGA or the first line
  patientBarcode=substring(patientBarcode,1,12)

  fn=paste(folderPath,folderList[i],sep="/")
  fmDat=read.table(gzfile(fn),sep="\t",header=FALSE)

  #parsing
  res=lapply(as.character(fmDat$V8), function(x){
    x=unlist(strsplit(x,";"))
    y=x[grep("CSQ=",x)]
    y=substring(y,5)
    z=unlist(strsplit(y,","))
    t=strsplit(z,"\\|")
    t=do.call(rbind,t)
    colnames(t)=field_names[1:ncol(t)]
    
    mut=paste(t[,4],collapse = ",")
    mut=unlist(mut)
    mut=mut[!duplicated(mut)]
    mutgene=list(mutgene=mut)
    MUTGENEdata=list(MUTGENE=mut)

    outList=list(MUTGENE=mut)
    return(outList)
  })
  
  res1=sapply(res,function(x) x[1])
  res2=matrix(unlist(res1),ncol=length(res1[[1]]),byrow=TRUE)
  colnames(res2)=names(res1[[1]])
  MUTGENE[[i]]=res2


  res=paste(fmDat$V1,"_",fmDat$V2,sep="")
  mutPos[[i]]=res
  
  refer=paste(fmDat$V4)
  Reference_allele[[i]]=refer
  
  mutation=paste(fmDat$V5)
  Mutated_allele[[i]]=mutation
   
  patientID[[i]]=patientBarcode
  
}

Mutect2=list(MUTGENE=MUTGENE,Reference_allele=Reference_allele, Mutated_allele=Mutated_allele,mutPos=mutPos)
names(Mutect2$MUTGENE)<-patientID
names(Mutect2$Reference_allele)<-patientID
names(Mutect2$Mutated_allele)<-patientID
names(Mutect2$mutPos)<-patientID

length(MUTGENE)

mutData=Mutect2
save(mutData,file=out)