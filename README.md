# MUSEC
Quantification of mutant-allele expression at isoform level for multiple genes

## What is MAX?
______________________________________________________________________________________________________________________________________________
MAX is a novel method to quantify the Mutant-Allele eXpression (MAX) at isoform level from RNA-seq data. Devaloped by Wenjiang Deng. This method was futher optimised to qualtify the Mutant allele expression of multiple genes for multiple samples. 


## Pipeline of running MAX
______________________________________________________________________________________________________________________________________________

<img width="313" height="318" alt="image" src="https://github.com/user-attachments/assets/3aa7836e-b8e1-489b-bd5f-39b69d85ae4f" />

Before running MAX preparation needs to made so all the needed parameters and data is correct for input.

1. First the manifest files are obtained through Genomic Data Commons (GDC) data portal (https://portal.gdc.cancer.gov/analysis_page?app=CohortBuilder&tab=general) through this you can choose the project and cancer. Then [filter](./Preprocessing/Manifest_files_filtering) so only the samples with both RNA data and mutect data are kept.

2. Start the [download](./Preprocessing/Download) of RNA and Mutect files

3. Once mutect files are downloaded, prepare mutation list and fasta files through running the step 1-3 preprocessing Rscipts.

4. Once the RNA bamfiles have completed downloading convert them to [fastq](./Preprocessing/Download). 

6. Prepare parameter file which where the file paths are for input 

7. Run MAX



