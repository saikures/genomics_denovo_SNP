#!/bin/bash
#PBS -N bcftools_bam_to_vcf
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=8:mem=32gb
#PBS -j oe

# Load conda environment
source /storage/brno12-cerit/home/saikures/miniconda3/bin/activate variant_calling_snp_env

WORKDIR=/storage/brno12-cerit/home/saikures/vp3/all_trim
cd "$WORKDIR"

# Reference genome
REF="filtered_scaffolds_200.fasta"

# Index reference genome if not present
if [[ ! -f "${REF}.fai" ]]; then
    samtools faidx "$REF"
fi

SAMPLES=("Vp3" "AVMp6" "Avp101" "AVp79" "AVSp7")

for SAMPLE in "${SAMPLES[@]}"; do
    BAM="${SAMPLE}.sorted.bam"
    BAI="${SAMPLE}.sorted.bam.bai"
    VCF_RAW="bcf_variants_${SAMPLE}.bcftools.raw.vcf"
    VCF_FILTERED="bcf_variants_${SAMPLE}.bcftools.filtered.vcf"
    FLAGSTAT="$bcf_{SAMPLE}_flagstat.txt"
    STATS="$bcf_{SAMPLE}_stats.txt"

    echo "🔄 Processing $SAMPLE..."

    # Check BAM file
    if [[ ! -f "$BAM" ]]; then
        echo "❌ $BAM not found, skipping..."
        continue
    fi

    # Index BAM if necessary
    if [[ ! -f "$BAI" ]]; then
        echo "⚠️ BAM index not found for $SAMPLE, indexing..."
        samtools index "$BAM"
    fi

    # Generate samtools stats
    samtools flagstat "$BAM" > "$FLAGSTAT"
    samtools stats "$BAM" > "$STATS"

    # Variant calling with bcftools
    bcftools mpileup -Ou -f "$REF" "$BAM" | \
        bcftools call -mv -Ov -o "$VCF_RAW"

    # Filter variants with QUAL > 20
    bcftools filter -i 'QUAL > 20' "$VCF_RAW" -o "$VCF_FILTERED"

    echo "✅ $SAMPLE bcftools variant calling done."
    echo "------------------------------------------"
done

echo "🎉 All samples processed using bcftools."
