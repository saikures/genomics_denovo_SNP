#!/bin/bash
#PBS -N diamond_blastx_ger
#PBS -l select=1:ncpus=64:mem=64gb
#PBS -l walltime=24:00:00
#PBS -j oe
#PBS -o "$PBS_O_WORKDIR"/diamond_blastx_ger.log

export OMP_NUM_THREADS=64

cd "$PBS_O_WORKDIR"
#mkdir -p DMND_TMP

module load diamond

diamond blastx \
  -d UnirefDMND/uniref50DMND.dmnd \
  -q filtered_scaffolds_200.fasta \
  -o final_filtered_scaffolds_blastx.txt \
  --outfmt 6 qseqid sseqid pident bitscore evalue length qstart qend sstart send \
  --more-sensitive \
  -t DMND_TMP \
  -k 1 \
  -p 64
