# MUSEC
Quantification of mutant-allele expression at isoform level for multiple genes

## What is MAX?
______________________________________________________________________________________________________________________________________________
MAX is a novel method to quantify the Mutant-Allele eXpression (MAX) at isoform level from RNA-seq data. Devaloped by Wenjiang Deng. This method was futher optimised to qualtify the Mutant allele expression of multiple genes for multiple samples. 

MAX requires two essential data components: the individual’s somatic mutation data and corresponding RNA-seq data. To quantify isoform-level expression of genes, both data types must be provided as input. Prior to running MAX, the mutation and RNA-seq data are preprocessed to ensure they are in the correct format for analysis. This preprocessing step prepares the data for accurate isoform-level quantification.


## Main Pipeline of running MAX
______________________________________________________________________________________________________________________________________________

<img width="313" height="318" alt="image" src="https://github.com/user-attachments/assets/3aa7836e-b8e1-489b-bd5f-39b69d85ae4f" />

Before running MAX preparation needs to made so all the needed parameters and data is correct for input. 
- Each step is detailed script is written in the respective folder/code.

1. First the manifest files are obtained through Genomic Data Commons (GDC) data portal (https://portal.gdc.cancer.gov/analysis_page?app=CohortBuilder&tab=general) through this you can choose the project and cancer. Then [filter](./Preprocessing/Manifest_files_filtering) so only the samples with both RNA data and mutect data are kept.
![Mutect_manifest](https://github.com/user-attachments/assets/34bc1ec1-6a93-4a0f-b677-c2d3a39b94c9)

2. Start the [download](./Preprocessing/Download) of RNA and Mutect files. 

3. Once mutect files are downloaded, prepare mutation list and fasta files through running the [step 1-3 preprocessing](Preprocessing/MutationList) Rscipts.

4. Once the RNA bamfiles have completed downloading convert them to [fastq](./Preprocessing/Download). 

5. Prepare parameter file which where the file paths are for input 

6. Run MAX



