#!/usr/bin/env bash
set -eo pipefail

# Authors:
#   - Arun Manoharan (arun@primediscoveries.com)
# Created August 2019 by Arun Manoharan
# Last Updated August 2019 by Arun Manoharan
#
# This wrapper script processes samples in parallel, performing downstream
# "merged samples" and "report data" analyses (see process_upstream.sh for
# upstream analyses). The samples have their data merged, and taxonomy
# classification and picrust predictions are generated. If "report data"
# analyses are enabled, report data is generated for each sample, treating all
# other samples as the "baseline". There is also an option to export data in
# common file formats (e.g. TSV, BIOM, FASTA) for compatibility with external
# tools.
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
#     --samples SAMPLES: file containing list of sample directories to process.
#       Each sample directory must contain "single-sample" results (e.g. from
#       running process_upstream.sh or run_marker_gene_pipeline.sh) [REQUIRED]
#
#     --ref_dir REF_DIR: directory containing reference files required by the
#       pipeline. See run_marker_gene_pipeline.sh for details [REQUIRED]
#
#     --deblur_trim_len DEBLUR_TRIM_LEN: deblur trim length used in upstream
#       processing [REQUIRED]
#
#     --output_dir OUTPUT_DIR: output directory in which to write results
#       [REQUIRED]
#
#     --target_subfragment TARGET_SUBFRAGMENT: target subfragment of the 16S
#       gene each sample's data was generated from [default: V4]
#
#     --report_data REPORT_DATA: if provided, generate report data for each
#       sample, treating all other samples as the "baseline" [default: report
#       data is not generated]
#
#     --export_data EXPORT_DATA: if provided, data are exported from qiime2
#       .qza files into common file formats (e.g. TSV, BIOM, FASTA) for
#       compatibility with external tools [default: data are not exported]
#
#     --jobs JOBS: number of jobs to run in parallel [default: 8]
#
# Example usage:
#     bash process_downstream.sh
#         --samples samples.txt
#         --ref_dir refs
#         --deblur_trim_len 150
#         --output_dir merged_samples

SAMPLES=""
REF_DIR=""
DEBLUR_TRIM_LEN=""
OUTPUT_DIR=""
TARGET_SUBFRAGMENT="V4"
REPORT_DATA="NO"
EXPORT_DATA="NO"
JOBS="8"

while [[ $# -gt 0 ]]; do
    key="${1}"

    case "${key}" in
        --samples)
            SAMPLES="${2}"
            shift 2
            ;;
        --ref_dir)
            REF_DIR="${2}"
            shift 2
            ;;
        --deblur_trim_len)
            DEBLUR_TRIM_LEN="${2}"
            shift 2
            ;;
        --output_dir)
            OUTPUT_DIR="${2}"
            shift 2
            ;;
        --target_subfragment)
            TARGET_SUBFRAGMENT="${2}"
            shift 2
            ;;
        --report_data)
            REPORT_DATA="YES"
            shift 1
            ;;
        --export_data)
            EXPORT_DATA="YES"
            shift 1
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

if [ -z "${SAMPLES}" ]; then
    echo "Error: --samples is a required option."
    exit 1
fi

if [ -z "${REF_DIR}" ]; then
    echo "Error: --ref_dir is a required option."
    exit 1
fi

if [ -z "${DEBLUR_TRIM_LEN}" ]; then
    echo "Error: --deblur_trim_len is a required option."
    exit 1
fi

if [ -z "${OUTPUT_DIR}" ]; then
    echo "Error: --output_dir is a required option."
    exit 1
fi

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)"

source activate qiime2-2019.4

python ${PIPELINE_DIR}/merged_samples/merge_samples.py \
    --sample_list "${SAMPLES}" \
    --output_dir "${OUTPUT_DIR}" \
    --denoiser deblur \
    --trim_len "${DEBLUR_TRIM_LEN}" \
    --KN 6:5000

conda deactivate

if [[ "${TARGET_SUBFRAGMENT}" == "V4" ]]; then
    CLASSIFIER="${REF_DIR}/gg-13-8-99-515-806-nb-classifier.qza"
    TAXONOMY_DIR="${OUTPUT_DIR}/taxonomy/gg-13-8-99-515-806-nb-classifier"
else
    CLASSIFIER="${REF_DIR}/gg-13-8-99-nb-classifier.qza"
    TAXONOMY_DIR="${OUTPUT_DIR}/taxonomy/gg-13-8-99-nb-classifier"
fi

bash ${PIPELINE_DIR}/merged_samples/classify_taxonomy.sh \
    --merged_dir "${OUTPUT_DIR}" \
    --classifier "${CLASSIFIER}" \
    --threads "${JOBS}"

bash ${PIPELINE_DIR}/merged_samples/picrust_predictions.sh \
    --merged_dir "${OUTPUT_DIR}" \
    --ref_seqs "${REF_DIR}/gg-13-8-97-rep-set.qza" \
    --perc_identity 0.97 \
    --threads "${JOBS}"

source activate qiime2-2019.4

# Functional abundances are currently computed by compute_functional_scores.py.
# We'll use the first sample ID in the file as a placeholder; the functional
# abundances themselves are computed independent of which sample ID is used for
# scoring. To minimize potential confusion, remove the score files computed by
# the script; we only want the functional abundance tables at this stage.
#
# In the future, it'd be better to compute functional abundance tables in a
# separate step to avoid this hack.
SAMPLE_ID="$(basename $(head -n1 "${SAMPLES}"))"
python ${PIPELINE_DIR}/report_data/compute_functional_scores.py \
    --picrust_dir "${OUTPUT_DIR}/picrust_predictions" \
    --sample_id "${SAMPLE_ID}" \
    --output_dir "${OUTPUT_DIR}/picrust_predictions/functional_abundances" \
    --save_abundances

rm "${OUTPUT_DIR}"/picrust_predictions/functional_abundances/gaba/score.tsv
rm "${OUTPUT_DIR}"/picrust_predictions/functional_abundances/vitamins/scores.tsv
rm "${OUTPUT_DIR}"/picrust_predictions/functional_abundances/scfa/scores.tsv

conda deactivate

if [[ "${EXPORT_DATA}" == "YES" ]]; then
    source activate qiime2-2019.4

    qiime tools export \
        --input-path "${OUTPUT_DIR}/asvs/feature_table.qza" \
        --output-path "${OUTPUT_DIR}/asvs/feature_table.biom" \
        --output-format BIOMV210Format

    biom convert \
        -i "${OUTPUT_DIR}/asvs/feature_table.biom" \
        -o "${OUTPUT_DIR}/asvs/feature_table.tsv" \
        --to-tsv

    qiime tools export \
        --input-path "${OUTPUT_DIR}/asvs/rep_seqs.qza" \
        --output-path "${OUTPUT_DIR}/asvs/rep_seqs.fasta" \
        --output-format DNAFASTAFormat

    qiime tools export \
        --input-path "${TAXONOMY_DIR}/taxonomy_classifications.qza" \
        --output-path "${TAXONOMY_DIR}/taxonomy_classifications.tsv" \
        --output-format TSVTaxonomyFormat

    # TODO some of the exporting could be parallelized (it's not a bottleneck
    # though, low priority)
    LEVELS=(2 3 4 5 6)
    LEVEL_NAMES=(phylum class order family genus)
    for ((i=0;i<${#LEVELS[@]};++i)); do
        LEVEL="${LEVELS[i]}"
        LEVEL_NAME="${LEVEL_NAMES[i]}"
        LEVEL_DIR="${TAXONOMY_DIR}/level_${LEVEL}_${LEVEL_NAME}"

        qiime tools export \
            --input-path "${LEVEL_DIR}/${LEVEL_NAME}_level_table.qza" \
            --output-path "${LEVEL_DIR}/${LEVEL_NAME}_level_table.biom" \
            --output-format BIOMV210Format

        biom convert \
            -i "${LEVEL_DIR}/${LEVEL_NAME}_level_table.biom" \
            -o "${LEVEL_DIR}/${LEVEL_NAME}_level_table.tsv" \
            --to-tsv
    done

    FUNCTIONS=(gaba vitamins scfa)
    for ((i=0;i<${#FUNCTIONS[@]};++i)); do
        FUNCTION="${FUNCTIONS[i]}"
        FUNCTION_DIR="${OUTPUT_DIR}/picrust_predictions/functional_abundances/${FUNCTION}"

        qiime tools export \
            --input-path "${FUNCTION_DIR}/functional_abundances.qza" \
            --output-path "${FUNCTION_DIR}/functional_abundances.biom" \
            --output-format BIOMV210Format

        biom convert \
            -i "${FUNCTION_DIR}/functional_abundances.biom" \
            -o "${FUNCTION_DIR}/functional_abundances.tsv" \
            --to-tsv
    done

    conda deactivate
fi

function generate_report_data {
    set -eo pipefail

    SAMPLE_DIR="$1"
    OUTPUT_DIR="$2"
    TAXONOMY_DIR="$3"
    PIPELINE_DIR="$4"
    DEBLUR_TRIM_LEN="$5"

    SAMPLE_ID="$(basename ${SAMPLE_DIR})"
    REPORT_DATA_DIR="${OUTPUT_DIR}/report_data/${SAMPLE_ID}/report_data"
    REPORT_DATA_LOGS_DIR="${OUTPUT_DIR}/report_data_logs/${SAMPLE_ID}"
    mkdir -p "${REPORT_DATA_DIR}"
    mkdir -p "${REPORT_DATA_LOGS_DIR}"

    source activate qiime2-2019.4

    python ${PIPELINE_DIR}/report_data/compute_diversity_score.py \
        --entropy_file "${OUTPUT_DIR}/kmers/entropy/6-mers_5000_entropy.tsv" \
        --sample_id "${SAMPLE_ID}" \
        --output_dir "${REPORT_DATA_DIR}/diversity"

    # TODO use --underflow_threshold once Eric decides on an appropriate value
    python ${PIPELINE_DIR}/report_data/compute_functional_scores.py \
        --picrust_dir "${OUTPUT_DIR}/picrust_predictions" \
        --sample_id "${SAMPLE_ID}" \
        --output_dir "${REPORT_DATA_DIR}/functional"

    # TODO use --underflow_threshold once Eric decides on an appropriate value
    python ${PIPELINE_DIR}/report_data/compute_probiotic_scores.py \
        --taxonomy_dir "${TAXONOMY_DIR}" \
        --sample_id "${SAMPLE_ID}" \
        --probiotic_taxa "${PIPELINE_DIR}/report_data/data/probiotic_taxa.txt" \
        --output_dir "${REPORT_DATA_DIR}/probiotic"

    python ${PIPELINE_DIR}/report_data/format_taxonomy_report_data.py \
        --taxonomy_dir "${TAXONOMY_DIR}" \
        --sample_id "${SAMPLE_ID}" \
        --output_dir "${REPORT_DATA_DIR}/taxonomy"

    # Copying the intermediate files for batch qc report
    ROOT_DIR=$(dirname "${PIPELINE_DIR}")
    source "${ROOT_DIR}/lib/util/batch_qc_files.sh"
    prepare_batch_qc_files "${SAMPLE_ID}" "${SAMPLE_DIR}" \
        "${REPORT_DATA_DIR}" "${TAXONOMY_DIR}" "${DEBLUR_TRIM_LEN}"

    mkdir "${REPORT_DATA_DIR}/qc"
    cp "${SAMPLE_DIR}/qc_stats.tsv" "${REPORT_DATA_DIR}/qc"
    cp -r "${SAMPLE_DIR}/quality_control/single_end_qc/fastqc" "${REPORT_DATA_DIR}/qc"

    mkdir "${REPORT_DATA_DIR}/results"
    mkdir "${REPORT_DATA_LOGS_DIR}/diversity"
    mkdir "${REPORT_DATA_LOGS_DIR}/functional"
    mkdir "${REPORT_DATA_LOGS_DIR}/probiotic"
    mkdir "${REPORT_DATA_LOGS_DIR}/taxonomy"

    mv "${REPORT_DATA_DIR}"/diversity/*.log "${REPORT_DATA_LOGS_DIR}/diversity"
    mv "${REPORT_DATA_DIR}/diversity" "${REPORT_DATA_DIR}/results"

    mv "${REPORT_DATA_DIR}"/functional/*.log "${REPORT_DATA_LOGS_DIR}/functional"
    mv "${REPORT_DATA_DIR}"/functional/* "${REPORT_DATA_DIR}/results"
    rmdir "${REPORT_DATA_DIR}/functional"

    mv "${REPORT_DATA_DIR}"/probiotic/*.log "${REPORT_DATA_LOGS_DIR}/probiotic"
    mv "${REPORT_DATA_DIR}/probiotic" "${REPORT_DATA_DIR}/results"

    mv "${REPORT_DATA_DIR}"/taxonomy/*.log "${REPORT_DATA_LOGS_DIR}/taxonomy"
    mv "${REPORT_DATA_DIR}"/taxonomy/* "${REPORT_DATA_DIR}/results"
    rmdir "${REPORT_DATA_DIR}/taxonomy"

    conda deactivate
}
export -f generate_report_data

if [[ "${REPORT_DATA}" == "YES" ]]; then
    mkdir "${OUTPUT_DIR}/report_data"
    mkdir "${OUTPUT_DIR}/report_data_logs"

    parallel --max-args 1 --no-notice -j "${JOBS}" generate_report_data {} "${OUTPUT_DIR}" "${TAXONOMY_DIR}" "${PIPELINE_DIR}" "${DEBLUR_TRIM_LEN}" < "${SAMPLES}"
fi
