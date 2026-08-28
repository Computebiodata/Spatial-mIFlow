#!/usr/bin/env bash
#SBATCH --job-name=Spatial-mIFlow
#SBATCH --partition=long
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=240:00:00
#SBATCH --output=logs/Spatial-mIFlow_%j.out
#SBATCH --error=logs/Spatial-mIFlow_%j.err
set -euo pipefail

# ==============================================================================
# SLURM launcher for the COMET multiplex IF pipeline
# ==============================================================================
#
# Usage (from the project root, so the relative logs/ path above resolves):
#   sbatch run_pipeline_slurm.sh
#   sbatch run_pipeline_slurm.sh -resume
#   sbatch run_pipeline_slurm.sh -w /custom/work/dir   # overrides the default below
#
#
# ==============================================================================

PROJECT_DIR="$(pwd)"
PIPELINE_DIR="${PROJECT_DIR}/sopa-plus"
PARAMS_FILE="${PROJECT_DIR}/params.yml"
MAIN_NF="${PIPELINE_DIR}/main.nf"

# Checks if nextflow available
if ! command -v nextflow >/dev/null 2>&1; then
    echo "ERROR: nextflow command not found. Load Nextflow before submitting this script."
    exit 127
fi

# Checks if program and params paths exist
for f in "${MAIN_NF}" "${PARAMS_FILE}"; do
    if [[ ! -s "${f}" ]]; then
        echo "ERROR: required file not found or empty: ${f}"
        exit 1
    fi
done

# Sets workDir if not passed via extra args
WORKDIR_ARGS=()
HAS_WORKDIR_FLAG=false
for arg in "$@"; do
    if [[ "${arg}" == "-w" || "${arg}" == "-work-dir" ]]; then
        HAS_WORKDIR_FLAG=true
        break
    fi
done
if [[ "${HAS_WORKDIR_FLAG}" == false ]]; then
    WORKDIR_ARGS=(-w "/scratch/${USER}/comet_pipeline/work")
fi

# Fixes for conda environment issues with local clash and cache
export PYTHONNOUSERSITE=1
export PIP_NO_CACHE_DIR=1

echo "============================================================"
echo "Spatial-mIFlow workflow launcher"
echo "============================================================"
echo "Pipeline dir:  ${PIPELINE_DIR}"
echo "Params file:   ${PARAMS_FILE}"
echo "SLURM job ID:  ${SLURM_JOB_ID:-none (not running under sbatch)}"
echo "Extra args:    $*"
echo "============================================================"

nextflow -version
echo "============================================================"

nextflow run "${MAIN_NF}" \
    -profile slurm \
    -params-file "${PARAMS_FILE}" \
    "${WORKDIR_ARGS[@]}" \
    "$@"
