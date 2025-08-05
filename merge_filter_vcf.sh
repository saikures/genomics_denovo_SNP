#!/bin/bash

set -e  # Exit on any error
set -u  # Exit on undefined variables

echo "🚀 Starting VCF processing pipeline..."

# Step 1: Fix sample names, sort and index each VCF
echo "📋 Step 1: Fixing sample names, sorting and indexing VCFs..."
for vcf in variants_*.filtered.vcf; do
    if [[ ! -f "$vcf" ]]; then
        echo "⚠️  No matching VCF files found (variants_*.filtered.vcf)"
        exit 1
    fi
    
    # Extract sample name from filename (e.g., variants_AVMp6.filtered.vcf -> AVMp6)
    sample_name=$(basename "$vcf" | sed 's/variants_//; s/\.filtered\.vcf//')
    sorted_vcf="fb_${vcf%.vcf}.sorted.vcf"
    echo "🔄 Processing: $vcf → $sorted_vcf (sample: $sample_name)"
    
    # Fix sample name, sort and compress
    bcftools reheader -s <(echo "$sample_name") "$vcf" | \
    bcftools sort -O z -o "${sorted_vcf}.gz"
    
    # Index the compressed file
    bcftools index --tbi "${sorted_vcf}.gz"
done

# Step 2: Merge all sorted VCFs
echo "📋 Step 2: Merging VCFs..."

# Check if sorted files exist
if ! ls fb_variants_*.sorted.vcf.gz 1> /dev/null 2>&1; then
    echo "❌ No sorted VCF files found for merging"
    exit 1
fi

# Create file list and show what we're merging
ls fb_variants_*.sorted.vcf.gz > vcf_list.txt
echo "🔗 Merging $(cat vcf_list.txt | wc -l) VCF files..."
cat vcf_list.txt | while read file; do
    echo "  - $file"
done

# Use file list approach
bcftools merge -m all --force-samples -O z -o fb_merged_all.vcf.gz \
    -l vcf_list.txt

# Clean up
rm vcf_list.txt
bcftools index --tbi fb_merged_all.vcf.gz

# Step 3: Filter for consensus SNPs
echo "📋 Step 3: Filtering for consensus SNPs..."
echo "🧹 Extracting SNPs, filling tags, and filtering..."

bcftools view -v snps fb_merged_all.vcf.gz | \
  bcftools +fill-tags -- -t AC,AN,AF,MAF | \
  bcftools filter -e 'AC=0 || AC=AN' -O z -o fb_consensus_snps.vcf.gz

# Index the final output
bcftools index --tbi fb_consensus_snps.vcf.gz

# Generate summary statistics
echo "📊 Generating summary statistics..."
echo "=== SUMMARY ==="
echo "Original files processed: $(ls fb_variants_*.sorted.vcf.gz | wc -l)"
echo "Total variants in merged file:"
bcftools view -H fb_merged_all.vcf.gz | wc -l
echo "Consensus SNPs (excluding fixed variants):"
bcftools view -H fb_consensus_snps.vcf.gz | wc -l
echo "Samples in final VCF:"
bcftools query -l fb_consensus_snps.vcf.gz | wc -l

echo "✅ Pipeline complete!"
echo "📁 Final output: fb_consensus_snps.vcf.gz"
echo "📁 Index file: fb_consensus_snps.vcf.gz.tbi"
