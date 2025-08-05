#!/bin/bash
#PBS -N genomics_pipeline_master
#PBS -l select=1:ncpus=2:mem=8gb
#PBS -l walltime=72:00:00
#PBS -j oe
#PBS -o logs/pipeline_master.log

# Master submission script for Snakemake genomics pipeline
# This script sets up the environment and runs Snakemake with PBS cluster execution

# Set working directory
cd "$PBS_O_WORKDIR" || exit 1
WORKDIR="$PBS_O_WORKDIR"

echo "[INFO] Starting genomics pipeline at $(date)"
echo "[INFO] Working directory: $WORKDIR"

# Create necessary directories
mkdir -p logs
mkdir -p results/{fastqc,trimming,assembly,mapping,variants,repeats,annotation}

# Load conda environment for Snakemake
echo "[INFO] Loading Snakemake environment..."
source /storage/brno12-cerit/home/saikures/miniconda3/bin/activate snakemake_env

# Check if Snakemake is available
if ! command -v snakemake &> /dev/null; then
    echo "[ERROR] Snakemake not found. Please install or check conda environment."
    exit 1
fi

# Validate configuration files
echo "[INFO] Validating configuration..."
if [[ ! -f "config.yaml" ]]; then
    echo "[ERROR] config.yaml not found!"
    exit 1
fi

if [[ ! -f "cluster.yaml" ]]; then
    echo "[ERROR] cluster.yaml not found!"
    exit 1
fi

# Check input data directory
if [[ ! -d "data" ]]; then
    echo "[ERROR] Data directory not found. Please create 'data' directory with input files."
    exit 1
fi

# Print pipeline configuration
echo "[INFO] Pipeline configuration:"
echo "  - Snakefile: Snakefile"
echo "  - Config: config.yaml"
echo "  - Cluster config: cluster.yaml"
echo "  - Jobs: ${MAX_JOBS:-50}"
echo "  - Dry run: ${DRY_RUN:-false}"

# Set Snakemake parameters
SNAKEMAKE_ARGS=(
    --snakefile Snakefile
    --configfile config.yaml
    --cluster-config cluster.yaml
    --cluster "qsub {cluster}"
    --jobs ${MAX_JOBS:-50}
    --rerun-incomplete
    --keep-going
    --printshellcmds
    --reason
    --stats snakemake_stats.json
)

# Add dry-run if specified
if [[ "${DRY_RUN:-false}" == "true" ]]; then
    SNAKEMAKE_ARGS+=(--dry-run)
    echo "[INFO] Running in DRY-RUN mode"
fi

# Add force rerun if specified
if [[ "${FORCE_RERUN:-false}" == "true" ]]; then
    SNAKEMAKE_ARGS+=(--forceall)
    echo "[INFO] Forcing rerun of all rules"
fi

# Create cluster status script for better job monitoring
cat > cluster_status.py << 'EOF'
#!/usr/bin/env python3
import sys
import subprocess
import re

jobid = sys.argv[1]

try:
    result = subprocess.run(['qstat', '-f', jobid], capture_output=True, text=True)
    if result.returncode == 0:
        if 'job_state = R' in result.stdout:
            print("running")
        elif 'job_state = Q' in result.stdout:
            print("running")  # Treat queued as running for Snakemake
        elif 'job_state = C' in result.stdout:
            print("success")
        else:
            print("failed")
    else:
        print("failed")
except:
    print("failed")
EOF

chmod +x cluster_status.py

# Add cluster status checking
SNAKEMAKE_ARGS+=(--cluster-status ./cluster_status.py)

echo "[INFO] Starting Snakemake execution..."
echo "[INFO] Command: snakemake ${SNAKEMAKE_ARGS[*]}"

# Run Snakemake
snakemake "${SNAKEMAKE_ARGS[@]}"

SNAKEMAKE_EXIT_CODE=$?

echo "[INFO] Snakemake completed with exit code: $SNAKEMAKE_EXIT_CODE"

# Generate summary report
if [[ -f "snakemake_stats.json" ]]; then
    echo "[INFO] Generating pipeline summary..."
    python3 << 'EOF'
import json
import sys

try:
    with open('snakemake_stats.json', 'r') as f:
        stats = json.load(f)
    
    print("\n" + "="*50)
    print("PIPELINE EXECUTION SUMMARY")
    print("="*50)
    
    if 'rules' in stats:
        total_jobs = sum(rule_stats.get('output-files', 0) for rule_stats in stats['rules'].values())
        print(f"Total jobs completed: {len(stats['rules'])}")
        print(f"Total output files: {total_jobs}")
        
        print("\nRule execution summary:")
        for rule, rule_stats in stats['rules'].items():
            runtime = rule_stats.get('mean-runtime', 0)
            print(f"  {rule}: {runtime:.2f}s avg runtime")
    
    print("="*50)
    
except Exception as e:
    print(f"Could not parse statistics: {e}")
EOF
fi

# Clean up temporary files
rm -f cluster_status.py

# Final status
if [[ $SNAKEMAKE_EXIT_CODE -eq 0 ]]; then
    echo "[SUCCESS] Pipeline completed successfully at $(date)"
    echo "[INFO] Results available in: results/"
    echo "[INFO] Logs available in: logs/"
else
    echo "[ERROR] Pipeline failed with exit code: $SNAKEMAKE_EXIT_CODE"
    echo "[INFO] Check logs for details: logs/"
    exit $SNAKEMAKE_EXIT_CODE
fi
