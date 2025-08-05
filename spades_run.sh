#!/bin/bash
#PBS -N spades_assembly
#PBS -l select=1:ncpus=20:mem=128gb:scratch_local=400gb
#PBS -l walltime=72:00:00
#PBS -j oe

# Load Conda and activate environment
source /storage/brno12-cerit/home/saikures/miniconda3/etc/profile.d/conda.sh
conda activate spades_assembler_env

# Remote data details
REMOTE_USER="saikures"
REMOTE_HOST="storage-brno12-cerit.metacentrum.cz"
REMOTE_DATA_DIR="~/vp3/trimmomatic_output_minlen_80"
REMOTE_TARGET_DIR="~/vp3/trimmomatic_output_minlen_80/spades_assemblies"

# Ensure cleanup on exit
#trap "clean_scratch" TERM EXIT

# Go to scratch
cd "$SCRATCHDIR"

# Download input reads
scp ${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_DATA_DIR}/EN00006609_2_Vp3_1.paired.fq.gz .
scp ${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_DATA_DIR}/EN00006609_2_Vp3_2.paired.fq.gz .

# Run SPAdes
spades.py \
    -1 EN00006609_2_Vp3_1.paired.fq.gz \
    -2 EN00006609_2_Vp3_2.paired.fq.gz \
    -o EN00006609_spades_out_scp \
    --threads 20 \
    --memory 96

# Check result
if [ -s EN00006609_spades_out_scp/contigs.fasta ]; then
    # Create remote output dir if missing
    ssh ${REMOTE_USER}@${REMOTE_HOST} "mkdir -p ${REMOTE_TARGET_DIR}"

    # Transfer entire output folder (recursive)
    scp -r EN00006609_spades_out_scp ${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_TARGET_DIR}/

    echo "✅ SPAdes output copied successfully."
else
    echo "❌ SPAdes run failed or contigs missing."
    exit 1
fi

echo "🎉 SPAdes assembly + SCP completed."
