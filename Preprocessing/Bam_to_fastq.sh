#!/bin/bash
#SBATCH --job-name=bam_to_fastq                # Job name
#SBATCH --output=bam_to_fastq_%j.out           # Output file name with job ID
#SBATCH --error=bam_to_fastq_%j.err            # Error file name with job ID
#SBATCH --time=35:00:00
#SBATCH -p core
#SBATCH -n 16

ml samtools 
# Edit this as per your paths
input_dir="Path/to/RNA_BamFiles"
output_dir="/Path/to/RNA_fastq"
id_file="/path/to/Manifest.RNA_pIDs.txt"


echo "Starting BAM to FASTQ conversion on SLURM..."

# Create the output directory if it does not exist
mkdir -p "$output_dir"

# Skip the header in the RNA_test_setIDs.txt file and process each line
tail -n +2 "$id_file" | while IFS=$'\t' read -r id filename md5 size state patientID; do
  echo "processing BAM file for patient: $patientID"

  # Define the folder corresponding to the patient ID (same name as 'id' in RNA_test_setIDs.txt)
  patient_folder="$input_dir/$id"

  # Check if the folder exists
  if [ ! -d "$patient_folder" ]; then
    echo "ERROR: Folder for patient $patientID not found!"
    continue
  fi

  
expected_bam_suffix="rna_seq.transcriptome.gdc_realn.bam"

# Check if a valid BAM file exists in the patient folder
bam_file=$(find "$patient_folder" -type f -name "*$expected_bam_suffix" ! -name "*.partial" 2>/dev/null)

  # Check if the BAM file exists inside the patient folder
  bam_file="$patient_folder/$filename"
  if [ ! -f "$bam_file" ]; then
    echo "ERROR: BAM file not found for patient: $patientID"
    sleep 360 
    continue
  fi

  # Convert BAM to FASTQ using samtools, with output files labeled by patientID
 processed_log="$output_dir/processed_samples.log"

# Ensure log file exists
touch "$processed_log"

# process BAM files
tail -n +2 "$id_file" | while IFS=$'\t' read -r id filename md5 size state patientID; do
  bam_file="$input_dir/$id/$filename"

  # Skip if already processed
  if grep -q "$patientID" "$processed_log"; then
    echo "Skipping $patientID (already processed)."
    continue
  fi

  # Convert BAM to FASTQ
  echo "processing: $patientID"
  samtools fastq "$bam_file" \
    -1 "$output_dir/${patientID}_1.fastq.gz" \
    -2 "$output_dir/${patientID}_2.fastq.gz"

  # Verify success
  if [[ -f "$output_dir/${patientID}_1.fastq.gz" && -f "$output_dir/${patientID}_2.fastq.gz" ]]; then
    echo "$patientID $(date)" >> "$processed_log"
  else
    echo "ERROR: Conversion failed for $patientID."
  fi
done

# Check if FASTQ files were successfully created and verify their integrity
  if [ -f "$output_dir/${patientID}_1.fastq.gz" ] && [ -f "$output_dir/${patientID}_2.fastq.gz" ]; then
    echo "FASTQ files successfully created for patient: $patientID. Checking integrity..."

    # Use zcat to verify the integrity of the FASTQ files
    if zcat "$output_dir/${patientID}_1.fastq.gz" >/dev/null && zcat "$output_dir/${patientID}_2.fastq.gz" >/dev/null; then
      echo "FASTQ files for patient $patientID passed integrity check. Deleting BAM file..."
      rm -r "$patient_folder"
    else
      echo "ERROR: FASTQ integrity check failed for patient $patientID. Retaining BAM file for troubleshooting."
    fi
  else
    echo "ERROR: FASTQ conversion failed for patient: $patientID. BAM file retained."
  fi
done

echo "BAM to FASTQ conversion complete! check the log files to ensure everything went well "
