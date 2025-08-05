#!/bin/bash

module load gatk
module load samtools

# ====== CONFIGURATION ======
REF_FASTA="filtered_scaffolds_200.fasta"
DICT_FILE="filtered_scaffolds_200.dict"
BAM_INPUT="mapped.sorted.bam"
RG_BAM="vp3.rg.bam"
SORTED_BAM="vp3.sorted.bam"
RGID="vp3"
RGLB="lib1"
RGPL="ILLUMINA"
RGPU="unit1"
RGSM="vp3"

# ====== STEP 1: Create Sequence Dictionary ======
if [[ -f "$DICT_FILE" ]]; then
    echo "[INFO] Dictionary $DICT_FILE already exists. Skipping."
else
    echo "[INFO] Creating sequence dictionary..."
    gatk CreateSequenceDictionary \
        -R "$REF_FASTA"
fi

# ====== STEP 2: Add Read Groups ======
if [[ -f "$RG_BAM" ]]; then
    echo "[INFO] Read group BAM $RG_BAM already exists. Skipping."
else
    echo "[INFO] Adding read groups..."
    gatk AddOrReplaceReadGroups \
        -I "$BAM_INPUT" \
        -O "$RG_BAM" \
        --RGID "$RGID" \
        --RGLB "$RGLB" \
        --RGPL "$RGPL" \
        --RGPU "$RGPU" \
        --RGSM "$RGSM"
fi

# ====== STEP 3: Sort BAM ======
if [[ -f "$SORTED_BAM" ]]; then
    echo "[INFO] Sorted BAM $SORTED_BAM already exists. Skipping."
else
    echo "[INFO] Sorting BAM file..."
    gatk SortSam \
        -I "$RG_BAM" \
        -O "$SORTED_BAM" \
        --SORT_ORDER coordinate
fi

# ====== STEP 4: Index Sorted BAM ======
if [[ -f "$SORTED_BAM.bai" ]]; then
    echo "[INFO] BAM index $SORTED_BAM.bai already exists. Skipping."
else
    echo "[INFO] Indexing sorted BAM..."
    samtools index "$SORTED_BAM"
fi

# ====== STEP 5: Validate BAM ======
echo "[INFO] Validating sorted BAM file..."
gatk ValidateSamFile \
    -I "$SORTED_BAM" \
    --MODE SUMMARY

echo "[INFO] Pipeline completed."
