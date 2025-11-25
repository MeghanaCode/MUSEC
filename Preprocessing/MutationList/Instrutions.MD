# Getting the nesscerary input files

Below is script that runs all three steps to create the mutation list
  1. create a folder and store the step0 - step3 rscripts 
  2. create appropriate folders to store the outputs
- Take the paths in script below as an example 
Important to change the paths!!

#!/bin/bash
#SBATCH --job-name=preprocessing_workflow    # Job name
#SBATCH --output=preprocessing_%A.out       # Standard output log
#SBATCH --error=preprocessing_%A.err        # Standard error log
#SBATCH -p core
#SBATCH -n 16
#SBATCH --time=01:00:00

# Load required modules
module load R/4.4.1

# Step 1: Parse TCGA mutations
echo "Starting Step 1: Parsing TCGA mutations..."
Rscript /nfs/NGS/Meghana/preprocess/Step1_parse_TCGA_mutation_mutect2.R \
  manifestFile="/nfs/NGS/Meghana/ESCA/gdc_manifest/batches/mutect_batch_2.txt" \
  mutPath="/nfs/NGS/Meghana/ESCA/Batch_2/mutations/mutect" \
  out="/nfs/NGS/Meghana/ESCA/Batch_2/mutations/ESCA_Mutect2.RData"

if [ $? -ne 0 ]; then
  echo "Error in Step 1. Exiting..."
  exit 1
fi
echo "Step 1 completed successfully."

# Step 2: Extract exonic mutations
echo "Starting Step 2: Extracting exonic mutations..."
Rscript /nfs/NGS/Meghana/preprocess/Step2_extract_exonic_mutations.R \
  gtfSqliteFn="/nfs/NGS/Meghana/refs/hg38.refGene.sqlite" \
  mutDataIn="/nfs/NGS/Meghana/ESCA/Batch_2/mutations/ESCA_Mutect2.RData" \
  out="/nfs/NGS/Meghana/ESCA/Batch_2/mutations/ESCA_exonicMut.RData"

if [ $? -ne 0 ]; then
  echo "Error in Step 2. Exiting..."
  exit 1
fi
echo "Step 2 completed successfully."

# Step 3: Create mutation list file
echo "Starting Step 3: Creating mutation list file..."
Rscript /nfs/NGS/Meghana/preprocess/Step3_makeMutListFile.R \
  gtfFn="/nfs/NGS/Meghana/refs/hg38.refGene.gtf.gz" \
  gtfSqliteFn="/nfs/NGS/Meghana/refs/hg38.refGene.sqlite" \
  fastaFn="/nfs/NGS/Meghana/refs/refMrna.fa.gz" \
  exonicMutIn="/nfs/NGS/Meghana/ESCA/Batch_2/mutations/ESCA_exonicMut.RData" \
  geneIn="/nfs/NGS/Meghana/refs/10x_Annotated_Human_Pan_Cancer_Panel.txt" \
  out="/nfs/NGS/Meghana/ESCA/Batch_2/mutations/ESCA_mutation_list.txt" \
  faOut="/nfs/NGS/Meghana/ESCA/Batch_2/mutations/ESCA_tx.fa"

if [ $? -ne 0 ]; then
  echo "Error in Step 3. Exiting..."
  exit 1
fi
echo "Step 3 completed successfully."

# Workflow completed
echo "All steps completed successfully."
