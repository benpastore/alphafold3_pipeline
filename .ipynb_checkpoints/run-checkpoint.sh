#!/bin/bash

set -euo pipefail

#######################################
# User-provided variables
####################################### 
FASTA1="/fs/ess/PCON0160/ben/pipelines/alphafold3_pipeline/af3_inputs/E3_NEDD_SUMO_SCREEN/c_elegans_E3.fa"
FASTA2="/fs/ess/PCON0160/ben/pipelines/alphafold3_pipeline/af3_inputs/E3_NEDD_SUMO_SCREEN/uniprot_reformatted.fa"
OUTDIR="$PWD/af3_output/E3_NEDD_SUMO_SCREEN"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MAINNF="$SCRIPT_DIR/main.nf"
DATE_TAG=$(date +%Y%m%d)
FINAL_OUTPUT="$OUTDIR/merged_pae_scores.tsv"
WRAPPER_LOG="$SCRIPT_DIR/pipeline_launcher.log"

#######################################
# Step 1: Check for pipeline completion
#######################################
echo "[$(date)] Checking for output file: $FINAL_OUTPUT" | tee -a "$WRAPPER_LOG"

if [[ -f "$FINAL_OUTPUT" ]]; then
  echo "[$(date)] ✅ Pipeline complete. No resubmission needed." | tee -a "$WRAPPER_LOG"
  exit 0
fi

#######################################
# Step 2: Submit the pipeline job via Slurm
#######################################
echo "[$(date)] 🚀 Checking for existing job..." | tee -a "$WRAPPER_LOG"

existing_job=$(squeue --me --name=af3_run --format=%A --noheader | head -n 1)

if [[ -n "$existing_job" ]]; then
    echo "[$(date)] 🔁 Found existing job (ID: $existing_job). Cancelling before resubmission." | tee -a "$WRAPPER_LOG"
    scancel "$existing_job"
    sleep 10
fi

echo "[$(date)] 📤 Submitting new Nextflow job with -resume" | tee -a "$WRAPPER_LOG"

echo """#!/bin/bash
#SBATCH --job-name=af3_run
#SBATCH --output=$SCRIPT_DIR/af3_run_%j.out
#SBATCH --error=$SCRIPT_DIR/af3_run_%j.err
#SBATCH --account=PCON0160
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=64:00:00
#SBATCH --partition=nextgen
#SBATCH --cluster=ascend

set -euo pipefail

nextflow $MAINNF --fasta1 $FASTA1 --fasta2 $FASTA2 --outdir $OUTDIR -resume

""" > run.sbatch
sbatch run.sbatch

#######################################
# Step 3: Schedule the script to re-run in 23 hours if not already scheduled
#######################################
echo "[$(date)] 📆 Checking if script is already scheduled..." | tee -a "$WRAPPER_LOG"

if atq | grep -F "$(whoami)" | grep -F "$SCRIPT_DIR/run.sh" > /dev/null; then
    echo "[$(date)] ⏱ Script already scheduled. Skipping re-schedule." | tee -a "$WRAPPER_LOG"
else
    echo "sh $SCRIPT_DIR/run.sh" | at now + 63 hours
    echo "[$(date)] 📆 Scheduled next check in 63 hours." | tee -a "$WRAPPER_LOG"
fi
