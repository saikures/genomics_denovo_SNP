#!/bin/bash
#PBS -N all_read_mapping
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=20:mem=512gb:scratch_local=500gb
#PBS -j oe

# Activate conda environment
source /storage/brno12-cerit/home/saikures/miniconda3/bin/activate samtools_env

# Set working directory
cd "$PBS_O_WORKDIR" || exit 1
WORKDIR="$PBS_O_WORKDIR"

# Reference fasta
REFERENCE="filtered_scaffolds_200.fasta"

# Move to scratch directory
cd "$SCRATCHDIR" || exit 1
echo "[INFO] Running in scratch directory: $SCRATCHDIR"

# Copy necessary input files to scratch
echo "[INFO] Copying files to scratch..."
cp "$WORKDIR"/"$REFERENCE" "$SCRATCHDIR" || { echo "[ERROR] Failed to copy reference"; exit 1; }
cp "$WORKDIR"/trimmed_R1_*.fastq.gz "$SCRATCHDIR" || { echo "[ERROR] Failed to copy R1 files"; exit 1; }
cp "$WORKDIR"/trimmed_R2_*.fastq.gz "$SCRATCHDIR" || { echo "[ERROR] Failed to copy R2 files"; exit 1; }

# Index the reference if not already done
if [[ ! -f "${REFERENCE}.bwt" ]]; then
    echo "[INFO] Indexing reference..."
    bwa index "$REFERENCE" || { echo "[ERROR] BWA indexing failed"; exit 1; }
else
    echo "[INFO] Reference index found, skipping indexing"
fi

# List of sample suffixes
SAMPLES=("AVMp6" "Avp101" "AVp79" "AVSp7" "Vp3")

# Process each sample
for SAMPLE in "${SAMPLES[@]}"; do
    echo "[INFO] Processing sample: $SAMPLE"
    
    TRIM_R1="trimmed_R1_${SAMPLE}.fastq.gz"
    TRIM_R2="trimmed_R2_${SAMPLE}.fastq.gz"
    SAM_OUT="${SAMPLE}.sam"
    BAM_OUT="${SAMPLE}.sorted.bam"
    BAI_OUT="${SAMPLE}.sorted.bam.bai"
    
    # Check if input files exist
    if [[ ! -f "$TRIM_R1" || ! -f "$TRIM_R2" ]]; then
        echo "[WARNING] Input files for $SAMPLE not found, skipping..."
        continue
    fi
    
    # Skip if BAM already exists and is complete
    if [[ -f "$BAM_OUT" && -f "$BAI_OUT" ]]; then
        echo "[INFO] $BAM_OUT already exists, skipping..."
        continue
    fi
    
    # BWA mem alignment
    echo "[INFO] Aligning $SAMPLE with BWA mem..."
    bwa mem -t 20 -M "$REFERENCE" "$TRIM_R1" "$TRIM_R2" > "$SAM_OUT" || {
        echo "[ERROR] BWA alignment failed for $SAMPLE"
        continue
    }
    
    # Convert SAM to BAM, sort, and index
    echo "[INFO] Converting SAM to sorted BAM for $SAMPLE..."
    samtools view -@ 10 -bS "$SAM_OUT" | samtools sort -@ 10 -o "$BAM_OUT" || {
        echo "[ERROR] SAM to BAM conversion failed for $SAMPLE"
        continue
    }
    
    # Index BAM file
    echo "[INFO] Indexing BAM file for $SAMPLE..."
    samtools index "$BAM_OUT" || {
        echo "[ERROR] BAM indexing failed for $SAMPLE"
        continue
    }
    
    # Remove SAM file to save space
    rm "$SAM_OUT"
    echo "[INFO] Completed processing $SAMPLE"
done

# Copy results back to working directory
echo "[INFO] Copying results back to working directory..."
cp *.sorted.bam "$WORKDIR" || echo "[WARNING] Failed to copy some BAM files"
cp *.sorted.bam.bai "$WORKDIR" || echo "[WARNING] Failed to copy some BAI files"

# Generate alignment statistics
echo "[INFO] Generating alignment statistics..."
for SAMPLE in "${SAMPLES[@]}"; do
    BAM_FILE="${SAMPLE}.sorted.bam"
    if [[ -f "$WORKDIR/$BAM_FILE" ]]; then
        samtools flagstat "$WORKDIR/$BAM_FILE" > "$WORKDIR/${SAMPLE}_flagstat.txt"
        samtools stats "$WORKDIR/$BAM_FILE" > "$WORKDIR/${SAMPLE}_stats.txt"
    fi
done

echo "[INFO] Read mapping pipeline completed!"
