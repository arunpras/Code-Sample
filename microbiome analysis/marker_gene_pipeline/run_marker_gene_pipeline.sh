#!/usr/bin/env bash
set -eo pipefail

# Authors:
#   - Galen Xing (galen@primediscoveries.com)
#   - Arun Manoharan (arun@primediscoveries.com)
# Created May 2019 by Galen Xing
# Last Updated August 2019 by Arun Manoharan
#
# This wrapper script executes an end-to-end 16S marker gene analysis of a
# client's sample, starting with raw sequence data and ending with report data.
# The client's sample is compared to the provided set of baseline samples to
# produce various scores and comparative results presented in reports.
#
# IMPORTANT: This pipeline is configured for Illumina 16S sequence data ONLY.
#
# Requirements:
#     - Use install/create_conda_envs.sh to create necessary conda
#       environments.
#
#     - Download the following reference files into a single directory to use
#       with --ref_dir option (do not rename the files):
#
#       If samples are 16S V4 data:
#
#       https://data.qiime2.org/2019.4/common/gg-13-8-99-515-806-nb-classifier.qza
#
#       If samples are non-V4 16S data:
#
#       https://data.qiime2.org/2019.4/common/gg-13-8-99-nb-classifier.qza
#
#     - Download Greengenes 13_8 and extract the archive:
#
#       ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz
#
#       Import the 97% OTUs reference sequences into qiime2 .qza format, and
#       store the file in the --ref_dir directory:
#
#       qiime tools import --input-path gg_13_8_otus/rep_set/97_otus.fasta --output-path gg-13-8-97-rep-set.qza --type FeatureData[Sequence]
#
# Script options:
#     --sample_dir SAMPLE_DIR: path to client sample directory, must contain
#       'raw_seqs' subdirectory of single-end or paired-end fastq files. The
#       name of SAMPLE_DIR will be the sample ID used in all downstream
#       analyses [REQUIRED]
#
#     --deblur_trim_len DEBLUR_TRIM_LEN: trim length for deblur [REQUIRED]
#
#     --baseline BASELINE: file containing list of baseline sample directories
#       to compare to client sample. If not provided, only single-sample
#       results will be generated for the client sample, and the sample will
#       not be merged and compared to baseline samples [default: no baseline]
#
#     --ref_dir REF_DIR: directory containing reference files required by the
#       pipeline. This option is only required if providing baseline samples
#       [default: no reference directory]
#
#     --target_subfragment TARGET_SUBFRAGMENT: target subfragment of the 16S
#       gene this sample's data was generated from [default: V4]
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
#     --threads THREADS: number of threads to use for parallel processes
#       [default: 8]
#
# Example usage:
#     bash run_marker_gene_pipeline.sh
#         --sample_dir ERR1328338
#         --deblur_trim_len 150
#         --baseline baseline_samples.txt
#         --ref_dir refs

SAMPLE_DIR=""
DEBLUR_TRIM_LEN=""
BASELINE=""
REF_DIR=""
TARGET_SUBFRAGMENT="V4"
TRIM_515F="NO"
BLOOM_SEQS=""
THREADS="8"

# https://stackoverflow.com/a/14203146/3776794
while [[ $# -gt 0 ]]; do
    key="${1}"

    case "${key}" in
        --sample_dir)
            SAMPLE_DIR="${2}"
            shift 2
            ;;
        --deblur_trim_len)
            DEBLUR_TRIM_LEN="${2}"
            shift 2
            ;;
        --baseline)
            BASELINE="${2}"
            shift 2
            ;;
        --ref_dir)
            REF_DIR="${2}"
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
        --threads)
            THREADS="${2}"
            shift 2
            ;;
        *)
            echo "Error: Unknown option '${key}'. Please see script header for usage instructions."
            exit 1
            ;;
    esac
done

if [ -z "${SAMPLE_DIR}" ]; then
    echo "Error: --sample_dir is a required option."
    exit 1
fi

if [ -z "${DEBLUR_TRIM_LEN}" ]; then
    echo "Error: --deblur_trim_len is a required option."
    exit 1
fi

if [[ -z "${REF_DIR}" && ! -z "${BASELINE}" ]]; then
    echo "Error: --ref_dir is a required option when --baseline is provided."
    exit 1
fi

SAMPLE_DIR_CONTENTS="$(ls -1 "${SAMPLE_DIR}")"

if [[ "${SAMPLE_DIR_CONTENTS}" != "raw_seqs" ]]; then
    echo "Error: ${SAMPLE_DIR} directory must only contain a "raw_seqs" subdirectory. The subdirectory either does not exist, or there are additional files/directories present (e.g. from previous processing of this sample)."
    exit 1
fi

# Writes message with timestamp to log file and also displays on stdout. Will
# create log file if it doesn't already exist.
function log_message {
    set -eo pipefail

    MSG="$1"
    SCRIPT_LOG="${SAMPLE_DIR}/run_marker_gene_pipeline.log"

    if [ ! -f "${SCRIPT_LOG}" ]; then
        echo -e "run_marker_gene_pipeline.sh log file\n" > "${SCRIPT_LOG}"
    fi

    echo "$(date +%Y-%m-%d:%H:%M:%S) ${MSG}" | tee -a "${SCRIPT_LOG}"
}

PIPELINE_STARTED="$(date +%Y-%m-%d:%H:%M:%S)"
log_message "Pipeline execution started."

log_message "Writing file of input parameters started."
PARAMS_LOG="${SAMPLE_DIR}/run_marker_gene_pipeline_params.log"

echo "SAMPLE_DIR: ${SAMPLE_DIR}" > "${PARAMS_LOG}"
echo "DEBLUR_TRIM_LEN: ${DEBLUR_TRIM_LEN}" >> "${PARAMS_LOG}"
echo "BASELINE: ${BASELINE}" >> "${PARAMS_LOG}"
echo "REF_DIR: ${REF_DIR}" >> "${PARAMS_LOG}"
echo "TARGET_SUBFRAGMENT: ${TARGET_SUBFRAGMENT}" >> "${PARAMS_LOG}"
echo "TRIM_515F: ${TRIM_515F}" >> "${PARAMS_LOG}"
echo "BLOOM_SEQS: ${BLOOM_SEQS}" >> "${PARAMS_LOG}"
echo "THREADS: ${THREADS}" >> "${PARAMS_LOG}"

log_message "Writing file of input parameters completed."

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)"

log_message "single_sample/quality_control.sh started."
if [[ "${TRIM_515F}" == "YES" ]]; then
    bash ${PIPELINE_DIR}/single_sample/quality_control.sh \
        --sample_dir ${SAMPLE_DIR} \
        --target_gene 16S \
        --target_subfragment ${TARGET_SUBFRAGMENT} \
        --trim_515f \
        --threads ${THREADS}
else
    bash ${PIPELINE_DIR}/single_sample/quality_control.sh \
        --sample_dir ${SAMPLE_DIR} \
        --target_gene 16S \
        --target_subfragment ${TARGET_SUBFRAGMENT} \
        --threads ${THREADS}
fi
log_message "single_sample/quality_control.sh completed."

log_message "single_sample/generate_asvs.sh started."
if [ -z "${BLOOM_SEQS}" ]; then
    bash ${PIPELINE_DIR}/single_sample/generate_asvs.sh \
        --sample_dir ${SAMPLE_DIR} \
        --denoiser deblur \
        --trim_len ${DEBLUR_TRIM_LEN}
else
    bash ${PIPELINE_DIR}/single_sample/generate_asvs.sh \
        --sample_dir ${SAMPLE_DIR} \
        --denoiser deblur \
        --trim_len ${DEBLUR_TRIM_LEN} \
        --bloom_seqs ${BLOOM_SEQS}
fi
log_message "single_sample/generate_asvs.sh completed."

log_message "single_sample/generate_kmers.sh started."
bash ${PIPELINE_DIR}/single_sample/generate_kmers.sh \
    --asv_dir ${SAMPLE_DIR}/asvs/deblur_single_end_trim_len_${DEBLUR_TRIM_LEN} \
    --KN 6:5000
log_message "single_sample/generate_kmers.sh completed."

source activate qiime2-2019.4

log_message "single_sample/compute_kmer_entropy.sh started."
python ${PIPELINE_DIR}/single_sample/compute_kmer_entropy.py \
    --kmer_dir ${SAMPLE_DIR}/asvs/deblur_single_end_trim_len_${DEBLUR_TRIM_LEN}/kmers
log_message "single_sample/compute_kmer_entropy.sh completed."

log_message "lib/util/sample_qc_stats.py started."
python "$(dirname "${PIPELINE_DIR}")/lib/util/sample_qc_stats.py" \
    --qc_dir ${SAMPLE_DIR}/quality_control \
    --asv_dir ${SAMPLE_DIR}/asvs/deblur_single_end_trim_len_${DEBLUR_TRIM_LEN} \
    --output_path ${SAMPLE_DIR}/qc_stats.tsv
log_message "lib/util/sample_qc_stats.py completed."

if [ -z "${BASELINE}" ]; then
    log_message "Since --baseline was not provided, only single-sample results have been generated. The sample will not be merged and compared to baseline samples."
else
    log_message "merged_samples/merge_samples.py started."

    cp ${BASELINE} ${SAMPLE_DIR}/samples_merged.txt
    echo -e "\n${SAMPLE_DIR}" >> ${SAMPLE_DIR}/samples_merged.txt

    MERGE_DIR="${SAMPLE_DIR}/merged_with_baseline"

    python ${PIPELINE_DIR}/merged_samples/merge_samples.py --sample_list ${SAMPLE_DIR}/samples_merged.txt --output_dir ${MERGE_DIR} --denoiser deblur --trim_len ${DEBLUR_TRIM_LEN} --KN 6:5000

    mv ${SAMPLE_DIR}/samples_merged.txt ${MERGE_DIR}/

    log_message "merged_samples/merge_samples.py completed."

    log_message "merged_samples/classify_taxonomy.sh started."

    if [[ "${TARGET_SUBFRAGMENT}" == "V4" ]]; then
        CLASSIFIER="${REF_DIR}/gg-13-8-99-515-806-nb-classifier.qza"
        TAXONOMY_DIR="${MERGE_DIR}/taxonomy/gg-13-8-99-515-806-nb-classifier"
    else
        CLASSIFIER="${REF_DIR}/gg-13-8-99-nb-classifier.qza"
        TAXONOMY_DIR="${MERGE_DIR}/taxonomy/gg-13-8-99-nb-classifier"
    fi

    bash ${PIPELINE_DIR}/merged_samples/classify_taxonomy.sh \
        --merged_dir ${MERGE_DIR} \
        --classifier ${CLASSIFIER} \
        --threads ${THREADS}
    log_message "merged_samples/classify_taxonomy.sh completed."

    log_message "merged_samples/picrust_predictions.sh started."
    bash ${PIPELINE_DIR}/merged_samples/picrust_predictions.sh \
        --merged_dir ${MERGE_DIR} \
        --ref_seqs ${REF_DIR}/gg-13-8-97-rep-set.qza \
        --perc_identity 0.97 \
        --threads ${THREADS}
    log_message "merged_samples/picrust_predictions.sh completed."

    SAMPLE_ID="$(basename ${SAMPLE_DIR})"
    REPORT_DATA_DIR="${SAMPLE_DIR}/report_data"
    mkdir "${REPORT_DATA_DIR}"

    log_message "report_data/compute_diversity_score.py started."
    python ${PIPELINE_DIR}/report_data/compute_diversity_score.py \
        --entropy_file ${MERGE_DIR}/kmers/entropy/6-mers_5000_entropy.tsv \
        --sample_id "${SAMPLE_ID}" \
        --output_dir ${REPORT_DATA_DIR}/diversity
    log_message "report_data/compute_diversity_score.py completed."

    log_message "report_data/compute_functional_scores.py started."
    # TODO use --underflow_threshold once Eric decides on an appropriate value
    python ${PIPELINE_DIR}/report_data/compute_functional_scores.py \
        --picrust_dir ${MERGE_DIR}/picrust_predictions \
        --sample_id "${SAMPLE_ID}" \
        --output_dir ${REPORT_DATA_DIR}/functional
    log_message "report_data/compute_functional_scores.py completed."

    log_message "report_data/compute_probiotic_scores.py started."
    # TODO use --underflow_threshold once Eric decides on an appropriate value
    python ${PIPELINE_DIR}/report_data/compute_probiotic_scores.py \
        --taxonomy_dir ${TAXONOMY_DIR} \
        --sample_id "${SAMPLE_ID}" \
        --probiotic_taxa ${PIPELINE_DIR}/report_data/data/probiotic_taxa.txt \
        --output_dir ${REPORT_DATA_DIR}/probiotic
    log_message "report_data/compute_probiotic_scores.py completed."

    log_message "report_data/format_taxonomy_report_data.py started."
    python ${PIPELINE_DIR}/report_data/format_taxonomy_report_data.py \
        --taxonomy_dir ${TAXONOMY_DIR} \
        --sample_id "${SAMPLE_ID}" \
        --output_dir ${REPORT_DATA_DIR}/taxonomy
    log_message "report_data/format_taxonomy_report_data.py completed."

    # Reorganize report data output based on Michal's request (separate log
    # files from results, and store results in a flattened hierarchy). This
    # output format is required by Michal's report app/visualizers.
    log_message "Reorganizing report data output started."

    # Copying the intermediate files for batch qc report
    ROOT_DIR=$(dirname "${PIPELINE_DIR}")
    source "${ROOT_DIR}/lib/util/batch_qc_files.sh"
    prepare_batch_qc_files "${SAMPLE_ID}" "${SAMPLE_DIR}" \
        "${REPORT_DATA_DIR}" "${TAXONOMY_DIR}" "${DEBLUR_TRIM_LEN}"

    # Engineering will provide sample QC stats to report reviewers.
    mkdir "${REPORT_DATA_DIR}/qc"
    mv "${SAMPLE_DIR}/qc_stats.tsv" "${REPORT_DATA_DIR}/qc"
    mv "${SAMPLE_DIR}/quality_control/single_end_qc/fastqc" "${REPORT_DATA_DIR}/qc"

    mkdir "${REPORT_DATA_DIR}/results"

    # TODO consider providing *all* pipeline logs (organized into some sort of
    # hierarchy) with the report data in the future (could be useful to display
    # in report review UI). For now, we're just moving the report data logs
    # into a separate directory because it was confusing to display only those
    # terminal logs to report reviewers (it gives an incomplete picture of the
    # execution).
    mkdir "${SAMPLE_DIR}/report_data_logs"
    mkdir "${SAMPLE_DIR}/report_data_logs/diversity"
    mkdir "${SAMPLE_DIR}/report_data_logs/functional"
    mkdir "${SAMPLE_DIR}/report_data_logs/probiotic"
    mkdir "${SAMPLE_DIR}/report_data_logs/taxonomy"

    mv "${REPORT_DATA_DIR}"/diversity/*.log "${SAMPLE_DIR}/report_data_logs/diversity"
    mv "${REPORT_DATA_DIR}/diversity" "${REPORT_DATA_DIR}/results"

    mv "${REPORT_DATA_DIR}"/functional/*.log "${SAMPLE_DIR}/report_data_logs/functional"
    mv "${REPORT_DATA_DIR}"/functional/* "${REPORT_DATA_DIR}/results"
    rmdir "${REPORT_DATA_DIR}/functional"

    mv "${REPORT_DATA_DIR}"/probiotic/*.log "${SAMPLE_DIR}/report_data_logs/probiotic"
    mv "${REPORT_DATA_DIR}/probiotic" "${REPORT_DATA_DIR}/results"

    mv "${REPORT_DATA_DIR}"/taxonomy/*.log "${SAMPLE_DIR}/report_data_logs/taxonomy"
    mv "${REPORT_DATA_DIR}"/taxonomy/* "${REPORT_DATA_DIR}/results"
    rmdir "${REPORT_DATA_DIR}/taxonomy"

    log_message "Reorganizing report data output completed."
fi

conda deactivate

PIPELINE_COMPLETED="$(date +%Y-%m-%d:%H:%M:%S)"

echo "${PIPELINE_STARTED} Pipeline execution started." > "${SAMPLE_DIR}/pipeline_completed.log"
echo "${PIPELINE_COMPLETED} Pipeline execution completed." >> "${SAMPLE_DIR}/pipeline_completed.log"
log_message "Pipeline execution completed."
