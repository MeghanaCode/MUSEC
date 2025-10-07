#!/bin/bash
mut="";gtf="";ref="";workdir=$PWD;hg="hg19";CPUNUM=8
while getopts ":m:g:r:k:d:v:t:" opt
do
    case $opt in
                v)
        echo "hg version is:$OPTARG"
                hg=$OPTARG
        ;;
        m)
        echo "Mutation file is:$OPTARG"
                mut=$OPTARG
        ;;
        g)
        echo "GFT file is:$OPTARG"
                gtf=$OPTARG
        ;;
        r)
        echo "Wild-type reference is:$OPTARG"
                ref=$OPTARG
        ;;
        k)
        echo "Kmer length is:$OPTARG"
        kmer=$OPTARG
        ;;
                d)
        echo "Working directory is:$OPTARG"
                workdir=$OPTARG
        ;;
        t)
        echo "Number of threads:$OPTARG"
        CPUNUM=$OPTARG
        ;;
                :)
                echo "$varname"
                echo "the option -$OPTARG requires an arguement"
                exit 1
                ;;
        ?)
        echo "Unkonwn parameter"
        exit 1;;
    esac
done

cd $workdir

Rscript /path/to/R/MAX.R mut=$mut gtf=$gtf ref=$ref workdir=$workdir hg=$hg


## do not use $ref, it could contain the whole transcriptome
#cat Mutant3.fa $ref > WT_Mut_tx_ref.1.fa

## use only wild-type transcripts, so smaller
cat Mutant3.fa wild_type_cdna_sequences.fa > WT_Mut_tx_ref.1.fa


awk '{if($1~/^>/)printf("%s Gene_Num=G%s\n",$1,NR);if($1!~/^>/)print $0;}' WT_Mut_tx_ref.1.fa > WT_Mut_tx_ref.2.fa

awk '/^>/{f=!d[$1];d[$1]=1}f' WT_Mut_tx_ref.2.fa > WT_Mut_tx_ref.final.fa

rm WT_Mut_tx_ref.2.fa Mutant3.fa Mutant2.fa WT_Mut_tx_ref.1.fa  tmp_gtf.sqlite
echo -e "\n WT_Mut.fa is generated. \n"

mv WT_Mut_tx_ref.final.fa WT_Mut_reference.fa

######################################
################### Generate X matrix
######################################
# export LD_LIBRARY_PATH=/proj/snic2020-6-4/nobackup/private/wenjiang/ASE2020/XAEM-binary-0.1.1/lib:$LD_LIBRARY_PATH
# export PATH=/proj/snic2020-6-4/nobackup/private/wenjiang/ASE2020/XAEM-binary-0.1.1/bin:$PATH
# XAEM_home=/proj/snic2020-6-4/nobackup/private/wenjiang/ASE2020/XAEM-binary-0.1.1
# module load bioinfo-tools
# module load R_packages/4.0.0

chmod +x /path/to/bin/*

cd $workdir
mkdir WT_Mut_Ref
mv WT_Mut_reference.fa WT_Mut_Ref
# Index
/path/to/bin/TxIndexer -t WT_Mut_Ref/WT_Mut_reference.fa -o Index_reference -k $kmer -p $CPUNUM --force

## generate RNA-seq using polyester
mkdir xSimDat
#Rscript /path/to/R/genPolyesterSimulation.R WT_Mut_Ref/WT_Mut_reference.fa xSimDat
#/path/to/bin/GenTC -i Index_reference -l IU -1 xSimDat/sample_01_1.fasta -2 xSimDat/sample_01_2.fasta -p $CPUNUM -o $workdir

### generate RNA-seq using a fast tool
#/path/to/bin/simx_pair WT_Mut_Ref -o xSimDat -t $CPUNUM -d 2

nbatch=10
Rscript /path/to/R/splitFasta.R in=WT_Mut_Ref/WT_Mut_reference.fa k=$nbatch
for ((i=1; i <= $nbatch; i++))
do
    echo $i
    simFn=$(echo "xSimDat_"$i)
    foldFn=$(echo "fold_"$i)
 #   mkdir $simFn
    /path/to/bin/simx_pair $foldFn -o $simFn -t $CPUNUM -d 2
    cat $simFn/subbatch_sim_1.fasta.gz  >> xSimDat/WT_Mut_reference_sim_1.fasta.gz
    cat $simFn/subbatch_sim_2.fasta.gz  >> xSimDat/WT_Mut_reference_sim_2.fasta.gz

    rm -r $simFn
    rm -r $foldFn
done


/path/to/bin/GenTC -i Index_reference -l IU -1 <(gunzip -c xSimDat/WT_Mut_reference_sim_1.fasta.gz) -2 <(gunzip -c xSimDat/WT_Mut_reference_sim_2.fasta.gz) -p $CPUNUM -o $workdir

rm -r xSimDat

Rscript /path/to/R/buildCRP.R in=Mutated_Combined_eqclass.txt out=$workdir/X_matrix.RData workdir=$workdir

#keep eqclass of X_matrix
mv eqClass.txt raw_Xmatrix.eq

rm tmp*RData
rm *fasta
rm Mutated_Combined_eqclass.txt fragmentInfo.txt gene_tx_tmp.RData wild_type_subset_genes.fa
rm -rf LogDir
######################################
################### Generate X matrix
######################################
