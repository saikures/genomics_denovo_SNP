#!/bin/bash

#source /storage/brno12-cerit/home/users/saikures/miniconda3/bin/activate variant_calling_snp_env

# Set paths
REF="filtered_scaffolds_200.fasta"
BAM="mapped.sorted.bam"
VCF_RAW="variants_freebayes.vcf"
VCF_FILTERED="variants_freebayes.filtered.vcf"

# Index reference (only needed once)
samtools faidx "$REF"

# Run FreeBayes
freebayes -f "$REF" "$BAM" > "$VCF_RAW"

# Filter variants by quality (QUAL > 20)
vcffilter -f "QUAL > 20" "$VCF_RAW" > "$VCF_FILTERED"

echo "✅ FreeBayes variant calling and filtering complete."
