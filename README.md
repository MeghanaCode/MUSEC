# MUSEC
MUSEC (Mutant Isoform Explorer for Cancers) is a comprehensive resource of mutant and wild-type expression profiles across 23 cancer types from The Cancer Genome Atlas (TCGA). The web application of MUSEC for interactive exploration is accessible at https://www.meb.ki.se/shiny/truvu/MUSEC/.

We developed MAX2pan, an enhanced implementation of MAX designed for scalable, multi-gene analysis, and applied the tool to quantify the mutant and wild-type isoform expression in MUSEC. This github page provides details of MAX2pan and how to use this tool.

## What is MAX2pan?
______________________________________________________________________________________________________________________________________________
MAX is a method developed by Wenjiang Deng et al. to quantify the Mutant-Allele eXpression (MAX) at isoform level from RNA-seq data. MAX2pan was futher optimised from MAX to qualtify the Mutant allele expression of multiple genes for multiple samples. 

MAX requires two essential data components: the individual’s somatic mutation data and corresponding RNA-seq data. To quantify isoform-level expression of genes, both data types must be provided as input. Prior to running MAX, the mutation and RNA-seq data are preprocessed to ensure they are in the correct format for analysis. This preprocessing step prepares the data for accurate isoform-level quantification.


## Main Pipeline of running MAX2pan
______________________________________________________________________________________________________________________________________________

<img width="313" height="318" alt="image" src="https://github.com/user-attachments/assets/3aa7836e-b8e1-489b-bd5f-39b69d85ae4f" />

Before running MAX2pan preparation needs to made so all the needed parameters and data is correct for input. 
- Each step in detailed script is written in the respective folder/code.
- Folders and scripts migth labelled as MAX this is the same as MAX2pan, "MAX" labelling was used for simplicity 

1. First the manifest files are obtained through Genomic Data Commons (GDC) data portal (https://portal.gdc.cancer.gov/analysis_page?app=CohortBuilder&tab=general) through this you can choose the project and cancer. How to do this in found in detail on [GDC_manifest](GDC_manifest.md)

2. Then [filter](./Preprocessing/Manifest_files_filtering) so only the samples with both RNA data and mutect data are kept.

3. Start the [download](./Preprocessing/Download) of RNA and Mutect files. 

4. Once mutect files are downloaded, prepare mutation list and fasta files through running the [step 1-3 preprocessing](Preprocessing/MutationList) R scripts.

5. Once the RNA bamfiles have completed downloading convert them to [fastq](./Preprocessing/Download). 

6. Prepare parameter file which where the file paths are for input 

7. Run MAX2pan, download the MAX-binary-2.0.1 for all the necessary files. 


## Batch wise running for large data sets
_______________________________________________________________________________________________________________________________________________________

* Cancer types with total sample data exceeding 1 - 1.5tb was split into batches as the MAX2pan tool can not handlle larger data.

* The pipeline for batch wise running is the same as main pipeline. however there is extra step after 1.

1.2 split the manifest files into apporpriate number of batches for example total sample of 5tb --> split to 5 batches
    here you can find template for batch wise split

* After making the batches follow steps 2-7 of the main pipeline use the mainfest files for the specific batch.

8. Stitch the all the Rdata of MAX2pan into one file


## Running MAX2pan
_____________________________________________________________________________________________________________________________________________

* MAX2pan is better excuted in 4 parts
* Before running make a [parameterfile](parameterfile), adjust all paths to suite yours
* Run MAX2pan scripts 1-4 though with [RunningMAX](RunningMAX) script, change the step at the last line of code MAX-binary-2.0.1/runMAX_step
 






