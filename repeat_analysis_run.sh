#!/bin/bash
#PBS -N repeat_pipeline
#PBS -l select=1:ncpus=8:mem=64gb
#PBS -l walltime=48:00:00
#PBS -j oe
#PBS -o repeat_pipeline.log

# Load environment or modules
module load blast

source /storage/brno12-cerit/home/saikures/miniconda3/bin/activate variant_calling_snp_env

# Set variables
GENOME_FASTA="/storage/brno12-cerit/home/saikures/vp3/trimmomatic_output_minlen_80_V/spades_assemblies/EN00006609_spades_out_scp/repeatmasking/filtered_scaffolds_200.fasta"
DB_NAME="scaffoldsDB"
THREADS=8

cd "$PBS_O_WORKDIR"

#cp $GENOME_FASTA $SCRATCH_DIR

#echo "[INFO] Building RepeatModeler database..."
#buildDatabase -name $DB_NAME -engine rmblast $GENOME_FASTA

echo "[INFO] Running RepeatModeler..."
RepeatModeler -database $DB_NAME -pa $THREADS

# Wait for RepeatModeler to finish before moving on
wait

# The output library
REPEAT_LIB="${DB_NAME}/consensi.fa.classified"
if [ ! -f "$REPEAT_LIB" ]; then
  echo "[ERROR] RepeatModeler output library not found!"
  exit 1
fi

echo "[INFO] Running RepeatMasker with custom library..."
RepeatMasker -s -pa $THREADS -lib $REPEAT_LIB -xsmall $GENOME_FASTA

#echo "[INFO] Copying results back to working directory..."
#cp -r $SCRATCHDIR/* $PBS_O_WORKDIR/

#echo "[INFO] Cleaning up scratch directory..."
#trap "clean_scratch" TERM EXIT 

echo "[INFO] Repeat masking completed."

