#!/usr/bin/env bash
set -eo pipefail

# Authors:
#   - Arun Manoharan (arun@primediscoveries.com)
# Created February 2019 by Arun Manoharan

# This script computes a number of alpha and beta diversity metrics from a
# feature table and (optionally) a phylogenetic tree:
#
# Alpha diversity:
#
# - Observed OTUs
# - Shannon’s diversity index
# - Pielou’s Evenness index
# - Simpson index
# - Faith's Phylogenetic Diversity (only computed if a phylogenetic tree is provided)
#
# Beta diversity:
#
# - Jaccard distance
# - Bray-Curtis distance
# - Unweighted UniFrac distance (only computed if a phylogenetic tree is provided)
# - Weighted UniFrac distance (only computed if a phylogenetic tree is provided)
# - Generalized UniFrac (alpha=0.5) distance (only computed if a phylogenetic tree is provided)
#
# Additionally, PCoA ordination results are computed for each of the beta
# diversity distance matrices.
#
# Output files are stored in qiime2 .qza format.

# TODO
#
# - Support dynamic sampling depths in addition to static? For example, using a
#   sampling depth that includes all samples, or median sampling depth.

# Requirements:
#
# - FEATURE_TABLE can be any type of feature table, including ASV, picrust,
#   taxonomy, etc.
#
# - If supplying a phylogenetic tree, it's recommended to use the phylogenetic
#   tree and filtered ASV feature table produced by generate_phylogeny.sh.

# Program requirements:
#     qiime2-2019.4 conda environment (use install/create_conda_envs.sh to create one)

# Script options:
#     --feature_table FEATURE_TABLE: path to feature table in qiime2 .qza format [REQUIRED]
#     --sampling_depth SAMPLING_DEPTH: number of sequences to rarefy each sample to (even sampling depth) [REQUIRED]
#     --output_dir OUTPUT_DIR: output directory in which to write diversity metric results [REQUIRED]
#     --phylogeny PHYLOGENY: path to phylogenetic tree in qiime2 .qza format [OPTIONAL]
#     --threads THREADS: number of threads to use for parallel processes [default: 8]
#     --qiime2_env QIIME2_ENV: qiime2 conda environment name [default: 'qiime2-2019.4']

# Outputs:
#     ${OUTPUT_DIR}/
#         diversity_metrics.log: log of commands started and completed during running of this script
#         rarefied_table.qza: feature table rarefied to specified SAMPLING_DEPTH
#         alpha_diversity: directory containing alpha diversity metric results
#         beta_diversity: directory containing beta diversity metric results

# Example usage:
#     bash diversity_metrics.sh
#         --feature_table merged_samples/phylogeny/sepp_default_gg_13_8_99_tree/2_filtered_asvs/feature_table.qza
#         --sampling_depth 5000
#         --output_dir diversity_metrics_5000_depth
#         --phylogeny merged_samples/phylogeny/sepp_default_gg_13_8_99_tree/1_sepp_placement/insertion_tree.qza
#         --threads 8

FEATURE_TABLE=""
SAMPLING_DEPTH=""
OUTPUT_DIR=""
PHYLOGENY=""
THREADS="8"
QIIME2_ENV="qiime2-2019.4"

# https://stackoverflow.com/a/14203146/3776794
while [[ $# -gt 0 ]]; do
    key="${1}"

    case "${key}" in
        --feature_table)
            FEATURE_TABLE="${2}"
            shift 2
            ;;
        --sampling_depth)
            SAMPLING_DEPTH="${2}"
            shift 2
            ;;
        --output_dir)
            OUTPUT_DIR="${2}"
            shift 2
            ;;
        --phylogeny)
            PHYLOGENY="${2}"
            shift 2
            ;;
        --threads)
            THREADS="${2}"
            shift 2
            ;;
        --qiime2_env)
            QIIME2_ENV="${2}"
            shift 2
            ;;
        *)
            echo "Error: Unknown option '${key}'. Please see script header for usage instructions."
            exit 1
            ;;
    esac
done

if [ -z "${FEATURE_TABLE}" ]; then
    echo "Error: --feature_table is a required option."
    exit 1
fi

if [ -z "${SAMPLING_DEPTH}" ]; then
    echo "Error: --sampling_depth is a required option."
    exit 1
fi

if [ -z "${OUTPUT_DIR}" ]; then
    echo "Error: --output_dir is a required option."
    exit 1
fi

if [ -d "${OUTPUT_DIR}" ]; then
  echo "Error: ${OUTPUT_DIR} directory already exists. Please either remove the directory and rerun this script or specify a different --output_dir."
  exit 1
fi

mkdir -p ${OUTPUT_DIR}

# Writes message with timestamp to log file and also displays on stdout. Will
# create log file if it doesn't already exist.
function log_message {
    set -eo pipefail

    MSG="${1}"
    SCRIPT_LOG="${OUTPUT_DIR}/diversity_metrics.log"

    if [ ! -f "${SCRIPT_LOG}" ]; then
        echo -e "diversity_metrics.sh log file\n" > "${SCRIPT_LOG}"
    fi

    echo "$(date +%Y-%m-%d:%H:%M:%S) ${MSG}" | tee -a "${SCRIPT_LOG}"
}

# Save pipeline start timestamp for writing to pipeline_completed.log file at
# end of script.
PIPELINE_STARTED="$(date +%Y-%m-%d:%H:%M:%S)"
log_message "Pipeline execution started."

# Create directories for output
log_message "Creating output directories started."

#mkdir -p ${OUTPUT_DIR}/alpha_diversity
mkdir -p ${OUTPUT_DIR}/beta_diversity

log_message "Creating output directories completed."

# export a copy of input parameters
log_message "Writing file of input parameters started."
PARAMS_LOG="${OUTPUT_DIR}/diversity_metrics_params.log"

echo "FEATURE_TABLE: $FEATURE_TABLE" > "${PARAMS_LOG}"
echo "SAMPLING_DEPTH: $SAMPLING_DEPTH" >> "${PARAMS_LOG}"
echo "OUTPUT_DIR: $OUTPUT_DIR" >> "${PARAMS_LOG}"
echo "PHYLOGENY: $PHYLOGENY" >> "${PARAMS_LOG}"
echo "THREADS: $THREADS" >> "${PARAMS_LOG}"
echo "QIIME2_ENV: $QIIME2_ENV" >> "${PARAMS_LOG}"

log_message "Writing file of input parameters completed."

source activate ${QIIME2_ENV}

# Note: We're not using `qiime diversity core-metrics` or
# `qiime diversity core-metrics-phylogenetic` because those commands require
# sample metadata in order to produce Emperor ordination plots. We also want to
# compute additional diversity metrics that aren't included in the
# `core-metrics*` commands.

# Rarefy feature table to even sampling depth.
log_message "Rarefying feature table to even sampling depth (${SAMPLING_DEPTH} seqs/sample) started."

RAREFIED_TABLE="${OUTPUT_DIR}/rarefied_table.qza"

qiime feature-table rarefy \
    --i-table "${FEATURE_TABLE}" \
    --o-rarefied-table "${RAREFIED_TABLE}" \
    --p-sampling-depth "${SAMPLING_DEPTH}" \
    --verbose

log_message "Rarefying feature table to even sampling depth (${SAMPLING_DEPTH} seqs/sample) completed."

# Compute non-phylogenetic alpha diversity metrics.
#log_message "Computing non-phylogenetic alpha diversity metrics started."
#
#for METRIC in "observed_otus" "shannon" "pielou_e" "simpson"; do
#    qiime diversity alpha \
#        --i-table "${RAREFIED_TABLE}" \
#        --o-alpha-diversity "${OUTPUT_DIR}/alpha_diversity/${METRIC}_vector.qza" \
#        --p-metric "${METRIC}" \
#        --verbose
#done
#
#log_message "Computing non-phylogenetic alpha diversity metrics completed."
#
#if [ ! -z "${PHYLOGENY}" ]; then
#    # Compute phylogenetic alpha diversity metrics.
#    log_message "Computing phylogenetic alpha diversity metrics started."
#
#    for METRIC in "faith_pd"; do
#        qiime diversity alpha-phylogenetic \
#            --i-table "${RAREFIED_TABLE}" \
#            --i-phylogeny "${PHYLOGENY}" \
#            --o-alpha-diversity "${OUTPUT_DIR}/alpha_diversity/${METRIC}_vector.qza" \
#            --p-metric "${METRIC}" \
#            --verbose
#    done
#
#    log_message "Computing phylogenetic alpha diversity metrics completed."
#fi

# Compute non-phylogenetic beta diversity metrics.
#log_message "Computing non-phylogenetic beta diversity metrics started."
#
#for METRIC in "jaccard" "braycurtis"; do
#    qiime diversity beta \
#        --i-table "${RAREFIED_TABLE}" \
#        --o-distance-matrix "${OUTPUT_DIR}/beta_diversity/${METRIC}_distance_matrix.qza" \
#        --p-metric "${METRIC}" \
#        --p-n-jobs "${THREADS}" \
#        --verbose
#done
#
#log_message "Computing non-phylogenetic beta diversity metrics completed."

if [ ! -z "${PHYLOGENY}" ]; then
    # Compute phylogenetic beta diversity metrics.
    log_message "Computing phylogenetic beta diversity metrics started."

#    for METRIC in "unweighted_unifrac" "weighted_unifrac" "generalized_unifrac"; do
    for METRIC in "unweighted_unifrac" "weighted_unifrac"; do
        if [[ "${METRIC}" == "generalized_unifrac" ]]; then
            # Note: Can only use --p-alpha with the generalized_unifrac metric.
            qiime diversity beta-phylogenetic \
                --i-table "${RAREFIED_TABLE}" \
                --i-phylogeny "${PHYLOGENY}" \
                --o-distance-matrix "${OUTPUT_DIR}/beta_diversity/${METRIC}_distance_matrix.qza" \
                --p-metric "${METRIC}" \
                --p-alpha 0.5 \
                --p-n-jobs "${THREADS}" \
                --verbose
        else
            # Note: Can only use --p-alpha with the generalized_unifrac metric.
            qiime diversity beta-phylogenetic \
                --i-table "${RAREFIED_TABLE}" \
                --i-phylogeny "${PHYLOGENY}" \
                --o-distance-matrix "${OUTPUT_DIR}/beta_diversity/${METRIC}_distance_matrix.qza" \
                --p-metric "${METRIC}" \
                --p-n-jobs "${THREADS}" \
                --verbose
        fi
    done

    log_message "Computing phylogenetic beta diversity metrics completed."
fi

# Compute PCoA ordination on distance matrices.
log_message "Computing PCoA ordination on distance matrices started."

for DISTANCE_MATRIX in ${OUTPUT_DIR}/beta_diversity/*_distance_matrix.qza; do
    METRIC="$(basename "${DISTANCE_MATRIX}" _distance_matrix.qza)"

    qiime diversity pcoa \
        --i-distance-matrix "${DISTANCE_MATRIX}" \
        --o-pcoa "${OUTPUT_DIR}/beta_diversity/${METRIC}_pcoa_results.qza" \
        --verbose
done

log_message "Computing PCoA ordination on distance matrices completed."

# Deactivate qiime2 conda env.
conda deactivate

PIPELINE_COMPLETED="$(date +%Y-%m-%d:%H:%M:%S)"

echo "${PIPELINE_STARTED} Pipeline execution started." > "${OUTPUT_DIR}/pipeline_completed.log"
echo "${PIPELINE_COMPLETED} Pipeline execution completed." >> "${OUTPUT_DIR}/pipeline_completed.log"
log_message "Pipeline execution completed."
