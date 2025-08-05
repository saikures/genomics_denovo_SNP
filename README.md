# Genomics Pipeline - Snakemake Workflow

A comprehensive Snakemake pipeline for genomics analysis including quality control, trimming, de novo assembly, variant calling, repeat masking, and gene prediction.

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Configuration](#configuration)
- [Input Data](#input-data)
- [Pipeline Steps](#pipeline-steps)
- [Output Structure](#output-structure)
- [Running the Pipeline](#running-the-pipeline)
- [Monitoring and Troubleshooting](#monitoring-and-troubleshooting)
- [Advanced Usage](#advanced-usage)

## Overview

This pipeline processes paired-end Illumina sequencing data through a complete genomics workflow:

1. **Quality Control** - FastQC analysis before and after trimming
2. **Read Trimming** - FASTP and Trimmomatic for adapter removal and quality trimming
3. **De Novo Assembly** - SPAdes genome assembly
4. **Assembly Quality Assessment** - BUSCO completeness analysis
5. **Scaffold Filtering** - Coverage-based filtering of assembled scaffolds
6. **Read Mapping** - BWA alignment of reads to filtered assembly
7. **Variant Calling** - FreeBayes and BCFtools for SNP/indel detection
8. **Consensus Variant Analysis** - Merging and filtering variants across samples
9. **Repeat Masking** - RepeatModeler and RepeatMasker for repeat identification
10. **Gene Prediction** - AUGUSTUS for ab initio gene prediction
11. **Functional Annotation** - DIAMOND BLASTx against UniRef database

## Features

- **Modular Design**: Each step is independent and can be run separately
- **PBS Integration**: Native support for PBS/Torque job scheduling systems
- **Automatic Resource Management**: Optimized memory and CPU allocation per step
- **Scratch Directory Handling**: Efficient use of local scratch space on HPC clusters
- **Comprehensive Quality Control**: Multiple QC checkpoints throughout the pipeline
- **Flexible Configuration**: Easy customization via YAML configuration files
- **Resume Capability**: Pipeline can resume from interruption points
- **Detailed Logging**: Comprehensive logging for debugging and monitoring

## Requirements

### Software Dependencies

- **Snakemake** (≥6.0)
- **FastQC** (v0.12.1)
- **FASTP** (v0.25.0)
- **Trimmomatic** (v0.39)
- **SPAdes** (v4.2.0)
- **BUSCO** (v6.0.0)
- **BWA** (v0.7.19)
- **SAMtools** (v1.22)
- **BCFtools** (v1.22)
- **FreeBayes** (v1.3.10)
- **RepeatModeler**
- **RepeatMasker**
- **AUGUSTUS**
- **DIAMOND**
- **seqtk**
- **vcffilter**

### System Requirements

- **HPC Cluster**: PBS/Torque job scheduler
- **Memory**: Minimum 512GB available across cluster nodes
- **Storage**: Fast scratch storage recommended for temporary files
- **CPU**: Multi-core nodes (minimum 8 cores recommended)

## Installation

### 1. Clone the Repository

```bash
git clone <repository_url>
cd genomics-pipeline
```

### 2. Create Conda Environment

```bash
# Create base Snakemake environment
conda create -n snakemake_env snakemake pandas pyyaml

# Create variant calling environment
conda create -n variant_calling_snp_env \
    samtools bcftools freebayes vcflib diamond

# Create additional environments as needed
conda create -n samtools_env samtools bwa seqtk
```

### 3. Install Additional Software

Some software may need manual installation:

- **RepeatModeler/RepeatMasker**: Follow installation instructions from official sources
- **AUGUSTUS**: Install with species-specific training data
- **BUSCO databases**: Download required lineage databases

## Quick Start

### 1. Prepare Input Data

```bash
# Create data directory
mkdir -p data

# Place your paired-end FASTQ files in data/
# Files should be named: {sample}_1.fastq.gz and {sample}_2.fastq.gz
```

### 2. Configure Pipeline

```bash
# Edit config.yaml with your sample names and paths
vim config.yaml

# Edit cluster.yaml for your PBS system if needed
vim cluster.yaml
```

### 3. Test Pipeline (Dry Run)

```bash
# Test the pipeline without executing
DRY_RUN=true bash submit_pipeline.sh
```

### 4. Run Pipeline

```bash
# Submit to PBS scheduler
qsub submit_pipeline.sh

# Or run with custom parameters
MAX_JOBS=30 qsub submit_pipeline.sh
```

## Configuration

### config.yaml

Main configuration file containing:

- **Sample names**: List of sample identifiers
- **Software parameters**: Tool-specific settings
- **File paths**: Database and reference file locations
- **Resource allocation**: Default memory and time limits

Key sections to modify:

```yaml
# Update with your sample names
samples:
  - "your_sample_1"
  - "your_sample_2"

# Update database paths
diamond:
  database: "/path/to/your/diamond/database.dmnd"

# Update species for gene prediction
augustus:
  species: "your_species"
```

### cluster.yaml

PBS-specific configuration:

- **Resource requests**: Memory, CPU, walltime per rule
- **Scratch space**: Local storage requirements
- **Queue parameters**: PBS-specific options

## Input Data

### Required Files

1. **Paired-end FASTQ files**:
   - Format: `{sample}_1.fastq.gz` and `{sample}_2.fastq.gz`
   - Location: `data/` directory
   - Quality: Illumina paired-end sequencing data

2. **Reference databases**:
   - **DIAMOND database**: UniRef or custom protein database
   - **BUSCO lineage**: Appropriate lineage database for your organism
   - **Adapter sequences**: Trimmomatic adapter file (e.g., TruSeq3-PE.fa)

### Directory Structure

```
project/
├── data/
│   ├── sample1_1.fastq.gz
│   ├── sample1_2.fastq.gz
│   ├── sample2_1.fastq.gz
│   └── sample2_2.fastq.gz
├── config.yaml
├── cluster.yaml
├── Snakefile
└── submit_pipeline.sh
```

## Pipeline Steps

### 1. Quality Control (FastQC)
- **Input**: Raw FASTQ files
- **Output**: HTML quality reports
- **Purpose**: Assess read quality, adapter content, sequence composition

### 2. Read Trimming
- **FASTP**: Remove polyG tails, detect adapters, quality trimming
- **Trimmomatic**: Aggressive quality trimming, adapter removal
- **Output**: High-quality paired reads

### 3. De Novo Assembly (SPAdes)
- **Input**: Trimmed paired reads from first sample
- **Output**: Scaffolds and contigs
- **Parameters**: Optimized for genome assembly

### 4. Assembly Quality Assessment
- **BUSCO**: Assess completeness using single-copy orthologs
- **Coverage filtering**: Remove low-coverage scaffolds
- **Output**: Filtered, high-quality scaffolds

### 5. Read Mapping (BWA)
- **Input**: Trimmed reads from all samples
- **Reference**: Filtered scaffolds
- **Output**: Sorted, indexed BAM files

### 6. Variant Calling
- **FreeBayes**: Population-aware variant caller
- **BCFtools**: Alternative variant calling approach
- **Output**: High-quality SNPs and indels

### 7. Variant Analysis
- **Merging**: Combine variants across samples
- **Filtering**: Consensus SNPs, quality filtering
- **Output**: Population-level variant set

### 8. Repeat Masking
- **RepeatModeler**: De novo repeat identification
- **RepeatMasker**: Mask repeats in genome
- **Output**: Repeat-masked genome sequence

### 9. Gene Prediction (AUGUSTUS)
- **Input**: Repeat-masked genome
- **Species model**: Configurable species parameters
- **Output**: Predicted gene coordinates (GFF format)

### 10. Functional Annotation (DIAMOND)
- **Input**: Repeat-masked genome
- **Database**: UniRef50 protein database
- **Output**: Functional annotations for predicted sequences

## Output Structure

```
results/
├── fastqc/                    # Quality control reports
├── trimming/                  # Trimmed sequences
├── assembly/                  # Assembly files and statistics
│   ├── scaffolds.fasta
│   ├── filtered_scaffolds.fasta
│   └── busco_output/
├── mapping/                   # Alignment files
│   ├── {sample}.sorted.bam
│   └── {sample}_flagstat.txt
├── variants/                  # Variant calling results
│   ├── variants_{sample}.*.vcf
│   └── consensus_snps.vcf
├── repeats/                   # Repeat analysis
│   └── filtered_scaffolds_200.fasta.masked
└── annotation/                # Gene prediction and annotation
    ├── predicted_genes.gff
    └── genes_filter_scaffolds_200.txt
```

## Running the Pipeline

### Basic Execution

```bash
# Submit to PBS queue
qsub submit_pipeline.sh
```

### Custom Parameters

```bash
# Run with specific number of jobs
MAX_JOBS=20 qsub submit_pipeline.sh

# Dry run to test configuration
DRY_RUN=true bash submit_pipeline.sh

# Force rerun of specific steps
FORCE_RERUN=true qsub submit_pipeline.sh
```

### Interactive Execution

```bash
# For testing or small datasets
conda activate snakemake_env
snakemake --configfile config.yaml --cores 8
```

## Monitoring and Troubleshooting

### Job Monitoring

```bash
# Check job status
qstat -u $USER

# Monitor pipeline progress
tail -f logs/pipeline_master.log

# Check individual job logs
ls logs/
tail -f logs/{rule}_{sample}.log
```

### Common Issues

1. **Job Failures**:
   - Check individual job logs in `logs/` directory
   - Verify resource allocation in `cluster.yaml`
   - Check software availability and versions

2. **File Not Found Errors**:
   - Verify input files exist in `data/` directory
   - Check file naming convention matches config
   - Ensure database paths are correct

3. **Memory Issues**:
   - Increase memory allocation in `cluster.yaml`
   - Check scratch space availability
   - Monitor node resource usage

4. **Network/Storage Issues**:
   - Verify scratch directory accessibility
   - Check file system quotas
   - Ensure network connectivity to databases

### Resuming Failed Runs

```bash
# Snakemake automatically resumes from last successful step
qsub submit_pipeline.sh

# Force rerun incomplete jobs
snakemake --rerun-incomplete --configfile config.yaml
```

## Advanced Usage

### Running Specific Steps

```bash
# Run only quality control
snakemake --configfile config.yaml results/fastqc/pretrim_sample1_1_fastqc.html

# Run up to assembly
snakemake --configfile config.yaml results/assembly/filtered_scaffolds.fasta

# Run variant calling only
snakemake --configfile config.yaml results/variants/consensus_snps.vcf
```

### Custom Resource Allocation

```bash
# Override default resources
snakemake --configfile config.yaml --resources mem_gb=256 --cores 32
```

### Pipeline Customization

1. **Adding New Rules**: Edit `Snakefile` to add new analysis steps
2. **Modifying Parameters**: Update `config.yaml` for tool-specific options
3. **Resource Optimization**: Adjust `cluster.yaml` for your HPC environment

### Performance Optimization

- **Parallel Execution**: Increase `MAX_JOBS` for more concurrent jobs
- **Scratch Usage**: Utilize fast local storage for I/O intensive steps
- **Resource Tuning**: Match resource requests to actual requirements
- **Job Dependency**: Leverage Snakemake's automatic dependency resolution

## Support and Citation

For issues, questions, or contributions, please:

1. Check the troubleshooting section
2. Review log files for error messages
3. Consult tool-specific documentation
4. Open an issue in the repository

When using this pipeline, please cite the relevant tools and databases used in your analysis.

## License

This pipeline is distributed under [LICENSE]. Individual tools retain their respective licenses.
