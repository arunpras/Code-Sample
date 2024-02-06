# Authors:
#   - Arun Manoharan(arun@primediscoveries.com)
# Created January 2020 by Arun Manoharan
# Last Updated January 2020 by Arun Manoharan

function prepare_batch_qc_files {
    set -eo pipefail

    SAMPLE_ID="$1"
    SAMPLE_DIR="$2"
    REPORT_DATA_DIR="$3"
    TAXONOMY_DIR="$4"
    DEBLUR_TRIM_LEN="$5"

    mkdir -p "${REPORT_DATA_DIR}/intermediates/kmers"
    mkdir -p "${REPORT_DATA_DIR}/intermediates/taxonomy"

    SAMPLE_TMP_DIR="${SAMPLE_DIR}/tmp"
    SAMPLE_ID_TEMP_FILE="${SAMPLE_TMP_DIR}/sampleid.tsv"
    mkdir -p "${SAMPLE_TMP_DIR}"
    echo -e "sampleid\n${SAMPLE_ID}" > "${SAMPLE_ID_TEMP_FILE}"

    # copy the kmer data
    find "${SAMPLE_DIR}/asvs/deblur_single_end_trim_len_${DEBLUR_TRIM_LEN}/kmers" -wholename "*.npz" | xargs -I{} \
        cp -u {} "${REPORT_DATA_DIR}/intermediates/kmers/"

    # copy the taxonomy data
    find "${TAXONOMY_DIR}" -wholename "*_level_table.qza" | while read -r TABLE; do
        BASE_QZA=$(basename "${TABLE}")
        BASE_TSV=${BASE_QZA/.qza/.tsv}

        # filter out baseline sample data and leave only the processed sample
        qiime feature-table filter-samples \
            --i-table "${TABLE}" \
            --o-filtered-table "${SAMPLE_TMP_DIR}/${BASE_QZA}" \
            --m-metadata-file "${SAMPLE_ID_TEMP_FILE}"

        # convert the resulting file to tsv
        qiime tools export \
            --input-path "${SAMPLE_TMP_DIR}/${BASE_QZA}" \
            --output-path "${SAMPLE_TMP_DIR}/feature_table.biom" \
            --output-format BIOMV210Format
        biom convert \
            -i "${SAMPLE_TMP_DIR}/feature_table.biom" \
            -o "${REPORT_DATA_DIR}/intermediates/taxonomy/${BASE_TSV}" \
            --to-tsv
    done

    rm -r "${SAMPLE_TMP_DIR}"
}
