# Genomics Pipeline Snakemake Workflow
# Author: Generated from pipeline documentation
# Description: Complete genomics workflow from raw reads to gene prediction

import os
from pathlib import Path

# Configuration
configfile: "config.yaml"

# Sample information
SAMPLES = config["samples"]
REFERENCE = config["reference"]

# Output directories
FASTQC_DIR = "results/fastqc"
TRIM_DIR = "results/trimming"
ASSEMBLY_DIR = "results/assembly"
MAPPING_DIR = "results/mapping"
VARIANT_DIR = "results/variants"
REPEAT_DIR = "results/repeats"
ANNOTATION_DIR = "results/annotation"

# Final outputs
rule all:
    input:
        # Quality control
        expand("{fastqc}/pretrim_{sample}_{read}_fastqc.html", 
               fastqc=FASTQC_DIR, sample=SAMPLES, read=["1", "2"]),
        expand("{fastqc}/posttrim_{sample}_{read}_fastqc.html", 
               fastqc=FASTQC_DIR, sample=SAMPLES, read=["1", "2"]),
        
        # Assembly and filtering
        f"{ASSEMBLY_DIR}/scaffolds.fasta",
        f"{ASSEMBLY_DIR}/filtered_scaffolds.fasta",
        f"{ASSEMBLY_DIR}/busco_output/short_summary.specific.euglenozoa_odb10.busco_output.txt",
        
        # Mapping
        expand("{mapping}/{sample}.sorted.bam", mapping=MAPPING_DIR, sample=SAMPLES),
        expand("{mapping}/{sample}.sorted.bam.bai", mapping=MAPPING_DIR, sample=SAMPLES),
        
        # Variant calling
        f"{VARIANT_DIR}/consensus_snps.vcf",
        
        # Repeat masking
        f"{REPEAT_DIR}/filtered_scaffolds_200.fasta.masked",
        
        # Gene prediction
        f"{ANNOTATION_DIR}/predicted_genes.gff",
        f"{ANNOTATION_DIR}/genes_filter_scaffolds_200.txt"

# Quality control - pre-trimming
rule fastqc_pretrim:
    input:
        r1 = "data/{sample}_1.fastq.gz",
        r2 = "data/{sample}_2.fastq.gz"
    output:
        html1 = f"{FASTQC_DIR}/pretrim_{{sample}}_1_fastqc.html",
        html2 = f"{FASTQC_DIR}/pretrim_{{sample}}_2_fastqc.html",
        zip1 = f"{FASTQC_DIR}/pretrim_{{sample}}_1_fastqc.zip",
        zip2 = f"{FASTQC_DIR}/pretrim_{{sample}}_2_fastqc.zip"
    params:
        outdir = FASTQC_DIR
    threads: 16
    resources:
        mem_gb = 32,
        time_min = 120
    shell:
        """
        mkdir -p {params.outdir}
        fastqc {input.r1} {input.r2} -t {threads} -o {params.outdir} --noextract
        # Rename outputs to match expected names
        mv {params.outdir}/$(basename {input.r1} .fastq.gz)_fastqc.html {output.html1}
        mv {params.outdir}/$(basename {input.r1} .fastq.gz)_fastqc.zip {output.zip1}
        mv {params.outdir}/$(basename {input.r2} .fastq.gz)_fastqc.html {output.html2}
        mv {params.outdir}/$(basename {input.r2} .fastq.gz)_fastqc.zip {output.zip2}
        """

# Trimming with FASTP
rule fastp_trimming:
    input:
        r1 = "data/{sample}_1.fastq.gz",
        r2 = "data/{sample}_2.fastq.gz"
    output:
        r1_trimmed = f"{TRIM_DIR}/{{sample}}_1.fastp.fq.gz",
        r2_trimmed = f"{TRIM_DIR}/{{sample}}_2.fastp.fq.gz",
        html = f"{TRIM_DIR}/{{sample}}_fastp.html",
        json = f"{TRIM_DIR}/{{sample}}_fastp.json"
    threads: 8
    resources:
        mem_gb = 16,
        time_min = 180
    shell:
        """
        mkdir -p {TRIM_DIR}
        fastp \
            -i {input.r1} \
            -I {input.r2} \
            -o {output.r1_trimmed} \
            -O {output.r2_trimmed} \
            --detect_adapter_for_pe \
            --trim_poly_g \
            --cut_right \
            --cut_right_window_size 4 \
            --cut_right_mean_quality 20 \
            --thread {threads} \
            --html {output.html} \
            --json {output.json}
        """

# Additional trimming with Trimmomatic
rule trimmomatic_trimming:
    input:
        r1 = f"{TRIM_DIR}/{{sample}}_1.fastp.fq.gz",
        r2 = f"{TRIM_DIR}/{{sample}}_2.fastp.fq.gz"
    output:
        r1_paired = f"{TRIM_DIR}/{{sample}}_1.paired.fq.gz",
        r1_unpaired = f"{TRIM_DIR}/{{sample}}_1.unpaired.fq.gz",
        r2_paired = f"{TRIM_DIR}/{{sample}}_2.paired.fq.gz",
        r2_unpaired = f"{TRIM_DIR}/{{sample}}_2.unpaired.fq.gz"
    params:
        adapters = config["trimmomatic"]["adapters"]
    threads: 8
    resources:
        mem_gb = 16,
        time_min = 180
    shell:
        """
        trimmomatic PE \
            -threads {threads} \
            -phred33 \
            {input.r1} {input.r2} \
            {output.r1_paired} {output.r1_unpaired} \
            {output.r2_paired} {output.r2_unpaired} \
            ILLUMINACLIP:{params.adapters}:2:30:10 \
            SLIDINGWINDOW:4:20 \
            LEADING:20 \
            TRAILING:20 \
            MINLEN:50
        """

# Post-trimming quality control
rule fastqc_posttrim:
    input:
        r1 = f"{TRIM_DIR}/{{sample}}_1.paired.fq.gz",
        r2 = f"{TRIM_DIR}/{{sample}}_2.paired.fq.gz"
    output:
        html1 = f"{FASTQC_DIR}/posttrim_{{sample}}_1_fastqc.html",
        html2 = f"{FASTQC_DIR}/posttrim_{{sample}}_2_fastqc.html",
        zip1 = f"{FASTQC_DIR}/posttrim_{{sample}}_1_fastqc.zip",
        zip2 = f"{FASTQC_DIR}/posttrim_{{sample}}_2_fastqc.zip"
    params:
        outdir = FASTQC_DIR
    threads: 16
    resources:
        mem_gb = 32,
        time_min = 120
    shell:
        """
        fastqc {input.r1} {input.r2} -t {threads} -o {params.outdir} --noextract
        # Rename outputs to match expected names
        mv {params.outdir}/$(basename {input.r1} .paired.fq.gz)_fastqc.html {output.html1}
        mv {params.outdir}/$(basename {input.r1} .paired.fq.gz)_fastqc.zip {output.zip1}
        mv {params.outdir}/$(basename {input.r2} .paired.fq.gz)_fastqc.html {output.html2}
        mv {params.outdir}/$(basename {input.r2} .paired.fq.gz)_fastqc.zip {output.zip2}
        """

# De novo assembly with SPAdes
rule spades_assembly:
    input:
        r1 = f"{TRIM_DIR}/{SAMPLES[0]}_1.paired.fq.gz",  # Using first sample for assembly
        r2 = f"{TRIM_DIR}/{SAMPLES[0]}_2.paired.fq.gz"
    output:
        scaffolds = f"{ASSEMBLY_DIR}/scaffolds.fasta",
        contigs = f"{ASSEMBLY_DIR}/contigs.fasta"
    params:
        outdir = f"{ASSEMBLY_DIR}/spades_out"
    threads: 20
    resources:
        mem_gb = 256,
        time_min = 1440  # 24 hours
    shell:
        """
        spades.py \
            -1 {input.r1} \
            -2 {input.r2} \
            -o {params.outdir} \
            --threads {threads}
        
        cp {params.outdir}/scaffolds.fasta {output.scaffolds}
        cp {params.outdir}/contigs.fasta {output.contigs}
        """

# Filter scaffolds by coverage
rule filter_scaffolds:
    input:
        scaffolds = f"{ASSEMBLY_DIR}/scaffolds.fasta"
    output:
        filtered = f"{ASSEMBLY_DIR}/filtered_scaffolds.fasta",
        coverage_tsv = f"{ASSEMBLY_DIR}/scaffold_coverage.tsv"
    params:
        min_coverage = config["filtering"]["min_coverage"]
    shell:
        """
        # Extract coverage information
        echo -e "Scaffold_ID\tCoverage" > {output.coverage_tsv}
        grep "^>" {input.scaffolds} | while read -r line; do
            scaffold=$(echo "$line" | cut -d' ' -f1 | sed 's/>//')
            coverage=$(echo "$line" | grep -oP "cov_\K[0-9]+(\.[0-9]+)?")
            echo -e "${{scaffold}}\t${{coverage}}" >> {output.coverage_tsv}
        done
        
        # Sort by coverage and filter
        sort -k2,2nr {output.coverage_tsv} > {output.coverage_tsv}.sorted
        awk -v min_cov={params.min_coverage} '$2 >= min_cov {{print $1}}' {output.coverage_tsv}.sorted > keep_scaffolds.txt
        
        # Extract filtered scaffolds
        seqtk subseq {input.scaffolds} keep_scaffolds.txt > {output.filtered}
        """

# BUSCO analysis
rule busco_analysis:
    input:
        assembly = f"{ASSEMBLY_DIR}/contigs.fasta"
    output:
        summary = f"{ASSEMBLY_DIR}/busco_output/short_summary.specific.euglenozoa_odb10.busco_output.txt"
    params:
        outdir = f"{ASSEMBLY_DIR}/busco_output",
        lineage = config["busco"]["lineage"]
    threads: 16
    resources:
        mem_gb = 64,
        time_min = 480
    shell:
        """
        busco -i {input.assembly} \
              -l {params.lineage} \
              -o busco_output \
              -m genome \
              --cpu {threads} \
              --out_path {ASSEMBLY_DIR}
        """

# BWA indexing
rule bwa_index:
    input:
        reference = f"{ASSEMBLY_DIR}/filtered_scaffolds.fasta"
    output:
        index = f"{ASSEMBLY_DIR}/filtered_scaffolds.fasta.bwt"
    resources:
        mem_gb = 32,
        time_min = 120
    shell:
        """
        bwa index {input.reference}
        """

# Read mapping
rule bwa_mapping:
    input:
        reference = f"{ASSEMBLY_DIR}/filtered_scaffolds.fasta",
        index = f"{ASSEMBLY_DIR}/filtered_scaffolds.fasta.bwt",
        r1 = f"{TRIM_DIR}/{{sample}}_1.paired.fq.gz",
        r2 = f"{TRIM_DIR}/{{sample}}_2.paired.fq.gz"
    output:
        bam = f"{MAPPING_DIR}/{{sample}}.sorted.bam",
        bai = f"{MAPPING_DIR}/{{sample}}.sorted.bam.bai",
        flagstat = f"{MAPPING_DIR}/{{sample}}_flagstat.txt",
        stats = f"{MAPPING_DIR}/{{sample}}_stats.txt"
    threads: 20
    resources:
        mem_gb = 64,
        time_min = 480
    shell:
        """
        mkdir -p {MAPPING_DIR}
        
        # BWA alignment
        bwa mem -t {threads} -M {input.reference} {input.r1} {input.r2} | \
            samtools view -@ 10 -bS - | \
            samtools sort -@ 10 -o {output.bam}
        
        # Index BAM
        samtools index {output.bam}
        
        # Generate statistics
        samtools flagstat {output.bam} > {output.flagstat}
        samtools stats {output.bam} > {output.stats}
        """

# Variant calling with FreeBayes
rule freebayes_variants:
    input:
        reference = f"{ASSEMBLY_DIR}/filtered_scaffolds.fasta",
        bam = f"{MAPPING_DIR}/{{sample}}.sorted.bam",
        bai = f"{MAPPING_DIR}/{{sample}}.sorted.bam.bai"
    output:
        vcf_raw = f"{VARIANT_DIR}/variants_{{sample}}.freebayes.raw.vcf",
        vcf_filtered = f"{VARIANT_DIR}/variants_{{sample}}.freebayes.filtered.vcf"
    resources:
        mem_gb = 32,
        time_min = 360
    shell:
        """
        mkdir -p {VARIANT_DIR}
        
        # Index reference if needed
        if [[ ! -f "{input.reference}.fai" ]]; then
            samtools faidx {input.reference}
        fi
        
        # FreeBayes variant calling
        freebayes -f {input.reference} {input.bam} > {output.vcf_raw}
        
        # Filter variants
        vcffilter -f "QUAL > 20" {output.vcf_raw} > {output.vcf_filtered}
        """

# Variant calling with BCFtools
rule bcftools_variants:
    input:
        reference = f"{ASSEMBLY_DIR}/filtered_scaffolds.fasta",
        bam = f"{MAPPING_DIR}/{{sample}}.sorted.bam",
        bai = f"{MAPPING_DIR}/{{sample}}.sorted.bam.bai"
    output:
        vcf_raw = f"{VARIANT_DIR}/variants_{{sample}}.bcftools.raw.vcf",
        vcf_filtered = f"{VARIANT_DIR}/variants_{{sample}}.bcftools.filtered.vcf"
    resources:
        mem_gb = 32,
        time_min = 360
    shell:
        """
        mkdir -p {VARIANT_DIR}
        
        # Index reference if needed
        if [[ ! -f "{input.reference}.fai" ]]; then
            samtools faidx {input.reference}
        fi
        
        # BCFtools variant calling
        bcftools mpileup -Ou -f {input.reference} {input.bam} | \
            bcftools call -mv -Ov -o {output.vcf_raw}
        
        # Filter variants
        bcftools filter -i 'QUAL > 20' {output.vcf_raw} -o {output.vcf_filtered}
        """

# Merge and filter consensus SNPs
rule merge_variants:
    input:
        vcfs = expand(f"{VARIANT_DIR}/variants_{{sample}}.bcftools.filtered.vcf", sample=SAMPLES)
    output:
        merged = f"{VARIANT_DIR}/merged_all.vcf",
        consensus = f"{VARIANT_DIR}/consensus_snps.vcf"
    resources:
        mem_gb = 16,
        time_min = 120
    shell:
        """
        # Sort and index VCFs
        for vcf in {input.vcfs}; do
            bcftools sort "$vcf" -o "${{vcf%.vcf}}.sorted.vcf"
            bcftools index --tbi "${{vcf%.vcf}}.sorted.vcf"
        done
        
        # Merge all sorted VCFs
        bcftools merge -m all -o {output.merged} {VARIANT_DIR}/variants_*.sorted.vcf
        bcftools index --tbi {output.merged}
        
        # Filter for consensus SNPs
        bcftools view -v snps {output.merged} | \
            bcftools +fill-tags | \
            bcftools filter -e 'AC!=0 && AC!=AN' -o {output.consensus}
        """

# Repeat masking
rule repeat_modeler:
    input:
        genome = f"{ASSEMBLY_DIR}/filtered_scaffolds.fasta"
    output:
        library = f"{REPEAT_DIR}/scaffoldsDB-families.fa"
    params:
        db_name = "scaffoldsDB",
        workdir = REPEAT_DIR
    threads: 8
    resources:
        mem_gb = 128,
        time_min = 1440  # 24 hours
    shell:
        """
        mkdir -p {params.workdir}
        cd {params.workdir}
        
        # Copy genome to repeat directory
        cp {input.genome} ./filtered_scaffolds_200.fasta
        
        # Build database
        BuildDatabase -name {params.db_name} -engine rmblast filtered_scaffolds_200.fasta
        
        # Run RepeatModeler
        RepeatModeler -database {params.db_name} -pa {threads}
        
        # Copy results
        cp {params.db_name}-families.fa {output.library}
        """

rule repeat_masker:
    input:
        genome = f"{ASSEMBLY_DIR}/filtered_scaffolds.fasta",
        library = f"{REPEAT_DIR}/scaffoldsDB-families.fa"
    output:
        masked = f"{REPEAT_DIR}/filtered_scaffolds_200.fasta.masked"
    params:
        workdir = REPEAT_DIR
    threads: 8
    resources:
        mem_gb = 64,
        time_min = 480
    shell:
        """
        cd {params.workdir}
        cp {input.genome} ./filtered_scaffolds_200.fasta
        
        RepeatMasker -s -pa {threads} -lib {input.library} -xsmall filtered_scaffolds_200.fasta
        
        cp filtered_scaffolds_200.fasta.masked {output.masked}
        """

# Gene prediction with AUGUSTUS
rule augustus_prediction:
    input:
        genome = f"{REPEAT_DIR}/filtered_scaffolds_200.fasta.masked"
    output:
        gff = f"{ANNOTATION_DIR}/predicted_genes.gff"
    params:
        species = config["augustus"]["species"]
    resources:
        mem_gb = 32,
        time_min = 240
    shell:
        """
        mkdir -p {ANNOTATION_DIR}
        augustus --species={params.species} {input.genome} > {output.gff}
        """

# DIAMOND BLASTx for gene annotation
rule diamond_blastx:
    input:
        genome = f"{REPEAT_DIR}/filtered_scaffolds_200.fasta.masked"
    output:
        blast_results = f"{ANNOTATION_DIR}/genes_filter_scaffolds_200.txt"
    params:
        diamond_db = config["diamond"]["database"]
    threads: 16
    resources:
        mem_gb = 64,
        time_min = 480
    shell:
        """
        mkdir -p {ANNOTATION_DIR}
        
        diamond blastx \
            -d {params.diamond_db} \
            -q {input.genome} \
            -o {output.blast_results} \
            --outfmt 6 qseqid sseqid pident bitscore evalue length qstart qend sstart send \
            --more-sensitive \
            -t DMND_TMP \
            -k 1 \
            --threads {threads}
        """
