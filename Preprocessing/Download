# Donwloading from GDC to server

For downloading the data first need to aquire GDC data transfer tool (https://gdc.cancer.gov/access-data/gdc-data-transfer-tool)
Also a token from GDC which gives the user authetication to download the data

Then with this type of bash script you can download the data 
#!/bin/bash
#SBATCH -J gdc_download              # Job name
#SBATCH --output=gdc_download_%j.out # Output file name with job ID
#SBATCH --error=gdc_download_%j.err  # Error file name with job ID
#SBATCH -t 12:00:00                  # Time limit
#SBATCH -p shared                    # Partition
#SBATCH -n 4                         # Number of cores
#SBATCH --mem=2G                     # Memory

# === User-defined variables ===
GDC_CLIENT="/path/to/gdc-client"          # Path to gdc-client executable
TOKEN="/path/to/gdc-user-token.json"      # Path to your GDC token file
MANIFEST="/path/to/gdc_manifest.txt"      # Path to your manifest file
OUTDIR="/path/to/output_directory"        # Directory where files will be saved

# === Run GDC download ===
$GDC_CLIENT download \
    -t $TOKEN \
    -m $MANIFEST \
    -d $OUTDIR

