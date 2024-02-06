#!/usr/bin/env bash
set -eo pipefail

# Authors:
#   - Arun Manoharan (arun@primediscoveries.com)
# Created August 2019 by Arun Manoharan
# Last Updated August 2019 by Arun Manoharan
#
# This wrapper script processes a directory of samples in parallel. Only the
# upstream "single-sample" portion of the 16S marker gene pipeline is executed
# for each sample (merging and downstream analyses are not performed; see
# process_downstream.sh to do that).
#
# This script is provided mainly as a convenience for R&D and
# manual/exploratory analysis. run_marker_gene_pipeline.sh is intended to be
# used directly in production, where each sample is processed independently of
# other samples.
#
# IMPORTANT: This pipeline is configured for Illumina 16S sequence data ONLY.
#
# Requirements:
#     - Use install/create_conda_envs.sh to create necessary conda
#       environments.
#
#     - This script requires the GNU `parallel` tool to be installed.
#
# Script options:
#     --samples_dir SAMPLES_DIR: directory of samples to process. Each sample
#       directory must contain a 'raw_seqs' subdirectory of single-end or
#       paired-end fastq files. The name of each sample directory will be the
#       sample ID used in all downstream analyses [REQUIRED]
#
#     --deblur_trim_len DEBLUR_TRIM_LEN: trim length for deblur [REQUIRED]
#
#     --target_subfragment TARGET_SUBFRAGMENT: target subfragment of the 16S
#       gene each sample's data was generated from [default: V4]
#
#     --trim_515f: if provided, trim the V4 515F primer sequence from each
#       forward read (i.e. the first 19 bases). A maximum of one base mismatch
#       is allowed. Reads that don't start with the primer sequence, or have
#       more than the allowed number of mismatches, will be discarded. CAUTION:
#       only use this option if reads begin with the V4 515F primer sequence
#       (this applies to Akesogen data); if primer sequences have already been
#       removed, do not use this option [default: no bases are trimmed]
#
#     --bloom_seqs BLOOM_SEQS: path to fasta file of bloom sequences to filter
#       from ASV table and representative sequences. A representative sequence
#       is identified as bloom if it is an exact match to a bloom sequence, or
#       if one is a prefix of the other (i.e. if representative sequence and
#       bloom sequence lengths differ) [default: no bloom filtering]
#
#     --jobs JOBS: number of samples to process in parallel [default: 8]
#
# Example usage:
#     bash process_upstream.sh
#         --samples_dir samples
#         --deblur_trim_len 150

SAMPLES_DIR=""
DEBLUR_TRIM_LEN=""
TARGET_SUBFRAGMENT="V4"
TRIM_515F="NO"
BLOOM_SEQS=""
JOBS="8"

while [[ $# -gt 0 ]]; do
    key="${1}"

    case "${key}" in
        --samples_dir)
            SAMPLES_DIR="${2}"
            shift 2
            ;;
        --deblur_trim_len)
            DEBLUR_TRIM_LEN="${2}"
            shift 2
            ;;
        --target_subfragment)
            TARGET_SUBFRAGMENT="${2}"
            shift 2
            ;;
        --trim_515f)
            TRIM_515F="YES"
            shift 1
            ;;
        --bloom_seqs)
            BLOOM_SEQS="${2}"
            shift 2
            ;;
        --jobs)
            JOBS="${2}"
            shift 2
            ;;
        *)
            echo "Error: Unknown option '${key}'. Please see script header for usage instructions."
            exit 1
            ;;
    esac
done

if [ -z "${SAMPLES_DIR}" ]; then
    echo "Error: --samples_dir is a required option."
    exit 1
fi

if [ -z "${DEBLUR_TRIM_LEN}" ]; then
    echo "Error: --deblur_trim_len is a required option."
    exit 1
fi

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)"

function process_sample {
    SAMPLE_DIR="$1"
    DEBLUR_TRIM_LEN="$2"
    TARGET_SUBFRAGMENT="$3"
    TRIM_515F="$4"
    BLOOM_SEQS="$5"
    PIPELINE_DIR="$6"

    if [ -f "${SAMPLE_DIR}/pipeline_completed.log" ]; then
        echo "Skipping '${SAMPLE_DIR}' because it has already been successfully processed."
    else
        echo "Processing '${SAMPLE_DIR}'..."

        declare -a SCRIPT_OPTS
        SCRIPT_OPTS=(--sample_dir "${SAMPLE_DIR}" --deblur_trim_len "${DEBLUR_TRIM_LEN}" --target_subfragment "${TARGET_SUBFRAGMENT}" --threads 1)

        if [[ "${TRIM_515F}" == "YES" ]]; then
            SCRIPT_OPTS+=(--trim_515f)
        fi

        if [ ! -z "${BLOOM_SEQS}" ]; then
            SCRIPT_OPTS+=(--bloom_seqs "${BLOOM_SEQS}")
        fi

        bash "${PIPELINE_DIR}/run_marker_gene_pipeline.sh" "${SCRIPT_OPTS[@]}" >> "${SAMPLE_DIR}.log" 2>&1

        if [ -f "${SAMPLE_DIR}/pipeline_completed.log" ]; then
            echo "Successfully processed '${SAMPLE_DIR}'."
        else
            echo "An error occurred while processing '${SAMPLE_DIR}'. Please see '${SAMPLE_DIR}.log' for details."
        fi
    fi
}
export -f process_sample

# Use --line-buffer to have status messages printed as they occur instead of at
# the end of each job.
find "${SAMPLES_DIR}" -mindepth 1 -maxdepth 1 -type d | parallel --max-args 1 --no-notice -j "${JOBS}" --line-buffer -q process_sample {} "${DEBLUR_TRIM_LEN}" "${TARGET_SUBFRAGMENT}" "${TRIM_515F}" "${BLOOM_SEQS}" "${PIPELINE_DIR}"
