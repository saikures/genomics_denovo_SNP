#!/bin/bash
#PBS -N bam_to_vcf
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=8:mem=32gb
#PBS -j oe

# Activate environment if needed
source /storage/brno12-cerit/home/saikures/miniconda3/bin/activate variant_calling_snp_env

WORKDIR=/storage/brno12-cerit/home/saikures/vp3/all_trim

cd "$WORKDIR"

# Reference genome
REF="filtered_scaffolds_200.fasta"

# Index reference if needed
if [[ ! -f "${REF}.fai" ]]; then
    samtools faidx "$REF"
fi

# Sample list
SAMPLES=("Vp3" "AVMp6" "Avp101" "AVp79" "AVSp7")

for SAMPLE in "${SAMPLES[@]}"; do
    BAM="${SAMPLE}.sorted.bam"
    BAI="${SAMPLE}.sorted.bam.bai"
    VCF_RAW="variants_${SAMPLE}.vcf"
    VCF_FILTERED="variants_${SAMPLE}.filtered.vcf"
    FLAGSTAT="${SAMPLE}_flagstat.txt"
    STATS="${SAMPLE}_stats.txt"

    # Skip if final filtered VCF exists
    if [[ -f "$VCF_FILTERED" ]]; then
        echo "✅ Filtered VCF for $SAMPLE already exists, skipping..."
        continue
    fi

    echo "🔄 Processing $SAMPLE..."

    # Check BAM and index
    if [[ ! -f "$BAM" ]]; then
        echo "❌ $BAM not found, skipping."
        continue
    fi

    if [[ ! -f "$BAI" ]]; then
        echo "⚠️ BAM index not found for $SAMPLE, generating..."
        samtools index "$BAM"
    fi

    # Generate flagstat
    samtools flagstat "$BAM" > "$FLAGSTAT"

    # Generate stats
    samtools stats "$BAM" > "$STATS"

    # Run FreeBayes
    freebayes -f "$REF" "$BAM" > "$VCF_RAW"

    # Filter variants by quality
    vcffilter -f "QUAL > 20" "$VCF_RAW" > "$VCF_FILTERED"

    echo "✅ $SAMPLE variant calling and stats done."
    echo "-------------------------------------------"
done

echo "🎉 All samples processed."
