#!/usr/bin/env bash
set -eo pipefail

# Authors:
#   - Sarah Nadeau (sarah@primediscoveries.com)
#   - Arun Manoharan (arun@primediscoveries.com)
# Created July 2018 by Sarah Nadeau
# Last Updated May 2019 by Arun Manoharan

# This script performs quality control (QC) on a single demultiplexed sample by:
#
# 1. Optionally trimming primers from beginning of reads
# 2. Trimming adapter contaminants from read ends
# 3. Quality trimming/filtering reads based on Phred score
# 4. Producing fastqc reports on raw and QC'd reads
#
# The script QC's single-end reads (forward reads only when analyzing paired-end data). Read merging is not yet supported.
#
# Note: PhiX reads and chimeric sequences are not filtered by this script. The denoisers (deblur and dada2) handle this filtering in `generate_asvs.sh`.
#
# Many of the QC steps performed in this script are based on qiita/deblur recommendations for performing meta-analyses (i.e. combining studies).
# There isn't a comprehensive or authoritative list of recommendations, and meta-analysis "best practices" are always evolving.
# The following resources were used to determine the current best practices for meta-analysis:
#
# - https://qiita.ucsd.edu/static/doc/html/processingdata/processing-recommendations.html
# - https://qiita.ucsd.edu/static/doc/html/processingdata/deblur_quality.html
# - https://cmi-workshop.readthedocs.io/en/latest/qiita-16S-analysis.html#analysis-of-deblur-processed-data
# - https://forum.qiime2.org/t/transferring-qiita-artifacts-to-qiime2/4790
# - https://forum.qiime2.org/t/running-deblur-on-paired-end-sequences/3387
# - https://forum.qiime2.org/t/combining-datasets-with-2-sets-of-primers/3073/12
# - https://genome.cshlp.org/content/23/10/1704.long

# To-DOs
#
# - Support optionally merging paired-end reads prior to quality filtering. Currently only the forward reads are used because this is what qiita supports,
#   and since the recommended read trim lengths are 90, 100, 150, and possibly 250, merged reads would end up being truncated to shorter length
#   anyways during denoising. Merging reads is desirable but we'll want a strategy to determine when too many reads are being discarded during merging,
#   because downstream pipelines rely on minimum sampling depths (e.g. kmer pipeline requires 5k seqs/sample, and sample comparison analyses may need
#   more seqs/sample to detect a signal depending on effect size of variable being investigated).

# Requirements:
#
# - The input sample directory name is used as the sample ID (e.g. sample name) in analyses.
#   The names of the sample's fastq files are *not* used to determine the sample ID.
#
# - For paired end reads, forward and reverse reads must be distinguished by _1 and _2 or _R1 and _R2 in filename suffixes.
#
# - Reads must have file extension ".fastq" or ".fq", and may optionally be gzipped (.gz).
#
# - Phred quality scores are assumed to be Phred33.
#
# - The QC has been optimized for Illumina sequence data. Exercise caution if
#   using this script with data produced by other sequencing platforms (e.g. 454, Sanger).

# Program requirements:
#     bbmap, qiime2-2019.4, and fastqc conda environments (use install/create_conda_envs.sh to create them)

# Script options:
#     --sample_dir SAMPLE_DIR: path to sample directory, must contain 'raw_seqs' subdirectory of single-end or paired-end fastq files. The name of SAMPLE_DIR should be the sample ID used in all downstream analyses [REQUIRED]
#     --target_gene TARGET_GENE: target gene for this sample's amplicon sequence data. Currently the only supported option is '16S' [REQUIRED]
#     --target_subfragment TARGET_SUBFRAGMENT: target subfragment of the TARGET_GENE (e.g. amplified region of the 16S gene). Supported options are: V1, V1V2, V1V3, V1V4, V1V5, V1V9, V2, V2V6, V3, V3V4, V3V5, V4, V4V5, V4V6, V5V6, V6, V6V7, V6V8, V8, V9 [REQUIRED]
#     --trim_515f: if provided, trim the V4 515F primer sequence from each
#       forward read (i.e. the first 19 bases). A maximum of one base mismatch
#       is allowed. Reads that don't start with the primer sequence, or have
#       more than the allowed number of mismatches, will be discarded. CAUTION:
#       only use this option if reads begin with the V4 515F primer sequence
#       (this applies to Akesogen data); if primer sequences have already been
#       removed, do not use this option [default: no bases are trimmed]
#     --threads THREADS: number of threads to use during QC steps that can be parallelized [default: 8]
#     --bbmap_env BBMAP_ENV: bbmap conda environment name [default: 'bbmap']
#     --qiime2_env QIIME2_ENV: qiime2 conda environment name [default: 'qiime2-2019.4']
#     --fastqc_env FASTQC_ENV: fastqc conda environment name [default: 'fastqc']

# Outputs:
#     ${SAMPLE_DIR}/quality_control/
#         qc.log: log of commands started and completed during running of this script
#         single_end_qc/
#             qc_seqs.qza: single-end QC'd sequences in qiime2 demultiplexed format
#             fastqc: directory containing fastqc reports for raw and QC'd sequences
#             read_counts.tsv: read counts remaining after each QC step

# Example usage:
#     bash quality_control.sh
#         --sample_dir /path/to/sample/directory
#         --target_gene 16S
#         --target_subfragment V4

SAMPLE_DIR=""
TARGET_GENE=""
TARGET_SUBFRAGMENT=""
TRIM_515F="NO"
THREADS="8"
BBMAP_ENV="bbmap"
QIIME2_ENV="qiime2-2019.4"
FASTQC_ENV="fastqc"

# https://stackoverflow.com/a/14203146/3776794
while [[ $# -gt 0 ]]; do
    key="${1}"

    case "${key}" in
        --sample_dir)
            SAMPLE_DIR="${2}"
            shift 2
            ;;
        --target_gene)
            TARGET_GENE="${2}"
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
        --threads)
            THREADS="${2}"
            shift 2
            ;;
        --bbmap_env)
            BBMAP_ENV="${2}"
            shift 2
            ;;
        --qiime2_env)
            QIIME2_ENV="${2}"
            shift 2
            ;;
        --fastqc_env)
            FASTQC_ENV="${2}"
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

if [ -z "${TARGET_GENE}" ]; then
    echo "Error: --target_gene is a required option."
    exit 1
elif [[ "${TARGET_GENE}" != "16S" ]]; then
    echo "Error: '${TARGET_GENE}' is not a supported target gene. The only supported target gene at this time is '16S'."
    exit 1
fi

if [ -z "${TARGET_SUBFRAGMENT}" ]; then
    echo "Error: --target_subfragment is a required option."
    exit 1
fi

case "${TARGET_SUBFRAGMENT}" in
    V1|V1V2|V1V3|V1V4|V1V5|V1V9|V2|V2V6|V3|V3V4|V3V5|V4|V4V5|V4V6|V5V6|V6|V6V7|V6V8|V8|V9)
        ;;
    *)
        echo "Error: '${TARGET_SUBFRAGMENT}' is not a supported target subfragment. Supported target subfragments are: V1, V1V2, V1V3, V1V4, V1V5, V1V9, V2, V2V6, V3, V3V4, V3V5, V4, V4V5, V4V6, V5V6, V6, V6V7, V6V8, V8, V9"
        exit 1
        ;;
esac

if [[ "${TRIM_515F}" == "YES" && "${TARGET_SUBFRAGMENT}" != "V4" ]]; then
    echo "Error: --trim_515f may only be used with --target_subfragment V4"
    exit 1
fi

QC_DIR="${SAMPLE_DIR}/quality_control"

if [ -d "${QC_DIR}" ]; then
    echo "Error: ${QC_DIR} directory already exists. Please either remove the directory and rerun this script, or specify a different --sample_dir."
    exit 1
fi

# loading common functions
# spellcheck source=../../lib/util/helpers.sh
. "$(dirname ${BASH_SOURCE[0]})/../../lib/util/helpers.sh"

mkdir -p ${QC_DIR}

# Writes message with timestamp to log file and also displays on stdout. Will
# create log file if it doesn't already exist.
function log_message {
    set -eo pipefail

    MSG="${1}"
    SCRIPT_LOG="${QC_DIR}/qc.log"

    if [ ! -f "${SCRIPT_LOG}" ]; then
        echo -e "quality_control.sh log file\n" > "${SCRIPT_LOG}"
    fi

    echo "$(date +%Y-%m-%d:%H:%M:%S) ${MSG}" | tee -a "${SCRIPT_LOG}"
}

# Save pipeline start timestamp for writing to pipeline_completed.log file at
# end of script.
PIPELINE_STARTED="$(date +%Y-%m-%d:%H:%M:%S)"
log_message "Pipeline execution started."

# Create directories for output
log_message "Creating output directories started."

mkdir -p ${QC_DIR}/single_end_qc

# Note: each QC step is performed within `tmp`. Final results are moved from
# `tmp` at the end of the script, and `tmp` is removed. This is done to
# decrease disk space usage, because each QC step has the potential to create a
# copy of the sample's sequences (and we're really only interested in the final
# QC'd output for downstream analyses). Stats are collected about the reads at
# each step before intermediate files are cleaned up.
mkdir -p ${QC_DIR}/single_end_qc/tmp

log_message "Creating output directories completed."

# export a copy of input parameters
log_message "Writing file of input parameters started."
PARAMS_LOG="${QC_DIR}/qc_params.log"

echo "SAMPLE_DIR: $SAMPLE_DIR" > "${PARAMS_LOG}"
echo "TARGET_GENE: $TARGET_GENE" >> "${PARAMS_LOG}"
echo "TARGET_SUBFRAGMENT: $TARGET_SUBFRAGMENT" >> "${PARAMS_LOG}"
echo "TRIM_515F: $TRIM_515F" >> "${PARAMS_LOG}"
echo "THREADS: $THREADS" >> "${PARAMS_LOG}"
echo "BBMAP_ENV: $BBMAP_ENV" >> "${PARAMS_LOG}"
echo "QIIME2_ENV: $QIIME2_ENV" >> "${PARAMS_LOG}"
echo "FASTQC_ENV: $FASTQC_ENV" >> "${PARAMS_LOG}"

log_message "Writing file of input parameters completed."

# Get the ID (name) of the sample we're processing, which is assumed to be the
# name of the sample directory.
log_message "Determining sample ID based on SAMPLE_DIR name started."

SAMPLE_ID="$(basename ${SAMPLE_DIR})"
log_message "Sample ID is: '${SAMPLE_ID}'"

log_message "Determining sample ID based on SAMPLE_DIR name completed."

# Perform basic validation of raw sequences directory to determine if sample is
# single or paired end, and if the necessary files are there.
log_message "Validating raw sequences directory started."

RAW_SEQS_DIR="${SAMPLE_DIR}/raw_seqs"

if [ ! -d "${RAW_SEQS_DIR}" ]; then
    echo "Error: ${RAW_SEQS_DIR} directory does not exist. The raw sequences directory must have exactly one or two fastq files containing this sample's single- or paired-end sequence data, respectively."
    exit 1
fi

FILE_COUNT="$(ls -1 "${RAW_SEQS_DIR}" | wc -l | sed 's/^[ \t]*//;s/[ \t]*$//')"

if [ "${FILE_COUNT}" -eq 1 ]; then
    READ_TYPE="SingleEnd"
    FILENAME="$(ls -1 "${RAW_SEQS_DIR}" | head -1)"

    case "${FILENAME}" in
        *.fastq|*.fastq.gz|*.fq|*.fq.gz)
            FWD_READS="${RAW_SEQS_DIR}/${FILENAME}"
            ;;
        *)
            echo "Error: Files in ${RAW_SEQS_DIR} must be .fastq, .fastq.gz, .fq, or .fq.gz files. Found this path: ${FILENAME}"
            exit 1
            ;;
    esac
elif [ "${FILE_COUNT}" -eq 2 ]; then
    READ_TYPE="PairedEnd"
    FWD_READS=""
    REV_READS=""
    for FILE in "${RAW_SEQS_DIR}"/*; do
        FILENAME="$(basename ${FILE})"

        case "${FILENAME}" in
            *_1.fastq|*_1.fastq.gz|*_1.fq|*_1.fq.gz|*_R1.fastq|*_R1.fastq.gz|*_R1.fq|*_R1.fq.gz)
                FWD_READS="${RAW_SEQS_DIR}/${FILENAME}"
                ;;
            *_2.fastq|*_2.fastq.gz|*_2.fq|*_2.fq.gz|*_R2.fastq|*_R2.fastq.gz|*_R2.fq|*_R2.fq.gz)
                REV_READS="${RAW_SEQS_DIR}/${FILENAME}"
                ;;
            *)
                echo "Error: Paired-end files in ${RAW_SEQS_DIR} must be .fastq, .fastq.gz, .fq, or .fq.gz files with read direction denoted by _1, _R1, _2, or _R2 suffix (e.g. seqs_R1.fastq). Found this path: ${FILENAME}"
                exit 1
                ;;
        esac
    done

    if [ -z "${FWD_READS}" ]; then
        echo "Error: Unable to find forward reads file in ${RAW_SEQS_DIR}. Forward reads must have _1 or _R1 suffix in the filename prior to the file extension (e.g. seqs_R1.fastq)."
        exit 1
    fi

    if [ -z "${REV_READS}" ]; then
        echo "Error: Unable to find reverse reads file in ${RAW_SEQS_DIR}. Reverse reads must have _2 or _R2 suffix in the filename prior to the file extension (e.g. seqs_R2.fastq)."
        exit 1
    fi
else
    echo "Error: ${RAW_SEQS_DIR} contains ${FILE_COUNT} files or directories. The raw sequences directory must have exactly one or two fastq files containing this sample's single- or paired-end sequence data, respectively."
    exit 1
fi

log_message "Sample read type: ${READ_TYPE}"

if [[ $READ_TYPE == "SingleEnd" ]]; then
    log_message "Single-end reads: ${FWD_READS}"
else
    log_message "Forward reads: ${FWD_READS}"
    log_message "Reverse reads: ${REV_READS}"
fi

log_message "Validating raw sequences directory completed."

# Function to count reads in a gzipped fastq file.
function count_reads {
    set -eo pipefail

    FASTQ="${1}"

    # Divide file line count by 4 to get number of fastq reads.
    echo "$(( $(gunzip -c ${FASTQ} | wc -l | cut -f1) / 4 ))"
}

log_message "Copying, renaming, and gzipping raw single-end/forward reads into single_end_qc directory started."

SINGLE_END_FILENAME="${SAMPLE_ID}.fastq.gz"

mkdir "${QC_DIR}/single_end_qc/tmp/raw"

case "${FWD_READS}" in
    *.gz)
        # Check that gzip file ends as expected (have had trouble with per-base
        # quality score filtering failing due to incomplete files).
        gzip -t "${FWD_READS}" || { echo "Error: ${FWD_READS} is not a valid gzip file. The file may be corrupted; please try re-downloading the file."; exit 1; }
        cp "${FWD_READS}" "${QC_DIR}/single_end_qc/tmp/raw/${SINGLE_END_FILENAME}"
        ;;
    *)
        cp "${FWD_READS}" "${QC_DIR}/single_end_qc/tmp/raw/${SAMPLE_ID}.fastq"
        gzip "${QC_DIR}/single_end_qc/tmp/raw/${SAMPLE_ID}.fastq"
        ;;
esac

RAW_READ_COUNT="$(count_reads "${QC_DIR}/single_end_qc/tmp/raw/${SINGLE_END_FILENAME}")"
log_message "Raw read count: ${RAW_READ_COUNT}"

log_message "Copying, renaming, and gzipping raw single-end/forward reads into single_end_qc directory completed."

source activate ${FASTQC_ENV}

# Run fastqc to generate quality reports for the raw sequences. These reports
# can be used to troubleshoot QC issues.
function runfastqc() {
    set -eo pipefail

    INPUT_FILE="${1}"
    OUTPUT_DIR="${2}"

    # Not using --threads because fastqc uses a separate thread per file, but
    # we're only processing a single file here.
    #
    # Extract the FastQC report data for further parsing (used to flag
    # samples).
    fastqc ${INPUT_FILE} \
        --outdir ${OUTPUT_DIR} \
        --extract
}

log_message "Running fastqc on raw sequences (single-end) started."

mkdir -p "${QC_DIR}/single_end_qc/tmp/fastqc/raw_seqs"

runfastqc "${QC_DIR}/single_end_qc/tmp/raw/${SINGLE_END_FILENAME}" "${QC_DIR}/single_end_qc/tmp/fastqc/raw_seqs"

log_message "Running fastqc on raw sequences (single-end) completed."

# Deactivate fastqc env
conda deactivate

function trim_primers {
    set -eo pipefail

    INPUT_FILE="${1}"
    OUTPUT_FILE="${2}"

    # The following command trims the V4 515F primer from the start of each
    # forward read.
    #
    # `--error-rate 0.053` was chosen to allow for a single base mismatch. The
    # primer is 19bp, and a single base mismatch has an error rate of
    # 1/19 = 0.0526... so we choose a value slightly higher than that as the
    # maximum acceptable error rate.
    #
    # `--no-indels` allows only mismatches in the alignment.
    #
    # `--times 1` specifies that the primer should only be found a single time.
    #
    # `--overlap 19` specifies that the primer must have 19bp of overlap with
    # the read.
    #
    # `--action trim` and `--discard-untrimmed` specifies that the primer is
    # trimmed from the reads, and reads without a matching primer are
    # discarded.
    #
    # Note: we use these criteria instead of simply trimming off the first 19
    # bases of each read because data generated by Akesogen (and potentially
    # other datasets) may contain a small proportion of reads that don't start
    # with the primer (they start with a frame-shifted or partial primer). In
    # this case we want to discard those reads.
    cutadapt \
        --front '^GTGYCAGCMGCCGCGGTAA' \
        --error-rate 0.053 \
        --no-indels \
        --times 1 \
        --overlap 19 \
        --action trim \
        --discard-untrimmed \
        --output "${OUTPUT_FILE}" \
        --cores "${THREADS}" \
        "${INPUT_FILE}" \
        2>&1 | tee -a "${QC_DIR}/cutadapt_primer_trimming.log"
}

mkdir "${QC_DIR}/single_end_qc/tmp/primers_trimmed"

if [[ "${TRIM_515F}" == "YES" ]]; then
    # cutadapt is installed in the qiime2 conda environment.
    source activate ${QIIME2_ENV}

    log_message "Trimming V4 515F primers from forward reads with cutadapt started."

    trim_primers "${QC_DIR}/single_end_qc/tmp/raw/${SINGLE_END_FILENAME}" "${QC_DIR}/single_end_qc/tmp/primers_trimmed/${SINGLE_END_FILENAME}"

    PRIMERS_TRIMMED_READ_COUNT="$(count_reads "${QC_DIR}/single_end_qc/tmp/primers_trimmed/${SINGLE_END_FILENAME}")"
    log_message "Read count after trimming primers: ${PRIMERS_TRIMMED_READ_COUNT}"

    log_message "Trimming V4 515F primers from forward reads with cutadapt completed."

    # Deactivate qiime2 env
    conda deactivate
else
    # No primer trimming requested; copy raw seqs into `primers_trimmed`
    # output directory for the next QC step.
    cp "${QC_DIR}/single_end_qc/tmp/raw/${SINGLE_END_FILENAME}" "${QC_DIR}/single_end_qc/tmp/primers_trimmed/${SINGLE_END_FILENAME}"
    PRIMERS_TRIMMED_READ_COUNT="N/A"
fi

source activate ${BBMAP_ENV}

# Locate bbmap resources directory in conda environment. The resources
# directory contains adapters to use with bbduk.
BBMAP_RESOURCES="$(dirname "$(which bbduk.sh)")/../opt/bbmap*/resources"
BBMAP_RESOURCES=( ${BBMAP_RESOURCES} )
BBMAP_RESOURCES="${BBMAP_RESOURCES[0]}"

# Trim adapter contaminants using bbduk.
#
# The parameters used here are based on the bbduk author's guides:
#
# - https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/data-preprocessing/
# - https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/
#
# The ddbuk author recommends quality filtering/trimming after adapter
# trimming. bbduk/bbduk2 can do this, but we choose to use qiime2's
# q2-quality-filter plugin in order to match qiita/deblur meta-analysis
# recommendations.
function trim_adapters {
    set -eo pipefail

    INPUT_FILE="${1}"
    OUTPUT_FILE="${2}"

    # bbduk should be able to auto-detect if quality scores are Phred33 or Phred64
    # ref: reference fasta sequences which, if a read kmer matches them, that kmer and anything the right should be trimmed (i.e. adapter contamination)
    # k: kmer length to match contaminants by - contaminants shorter than k will not be removed
    # tossbrokenreads will discard reads where sequence and quality score lengths differ
    #
    # Specify 1GB max memory for the JVM. We're overriding bbduk's default
    # automatic memory detection because it sometimes fails if this script
    # is running in multiple parallel processes (e.g. using the GNU `parallel`
    # command). 1GB should be more than enough memory for how we're using
    # bbduk. If we get out-of-memory errors in the future, consider making
    # `-Xmx` a script option that users can set based on their system's
    # RAM and number of parallel jobs they're executing.
    #
    # TODO investigate using `tpe` and `tbo` parameters when processing
    # paired-end reads.
    bbduk.sh \
    in=${INPUT_FILE} \
    out=${OUTPUT_FILE} \
    ref=${BBMAP_RESOURCES}/adapters.fa \
    ktrim=r \
    k=23 \
    mink=11 \
    hdist=1 \
    threads=${THREADS} \
    tossbrokenreads \
    -Xmx1g \
    2>&1 | tee -a ${QC_DIR}/bbduk_adapter_trimming.log
}

log_message "Adapter trimming single-end reads with bbduk started."

mkdir "${QC_DIR}/single_end_qc/tmp/adapters_trimmed"

trim_adapters "${QC_DIR}/single_end_qc/tmp/primers_trimmed/${SINGLE_END_FILENAME}" "${QC_DIR}/single_end_qc/tmp/adapters_trimmed/${SINGLE_END_FILENAME}"

ADAPTERS_TRIMMED_READ_COUNT="$(count_reads "${QC_DIR}/single_end_qc/tmp/adapters_trimmed/${SINGLE_END_FILENAME}")"
log_message "Read count after trimming adapters: ${ADAPTERS_TRIMMED_READ_COUNT}"

log_message "Adapter trimming single-end reads with bbduk completed."

# Deactivate bbmap env
conda deactivate

source activate ${QIIME2_ENV}

function import_qiime2_manifest {
    set -eo pipefail

    INPUT_FILE="${1}"
    OUTPUT_MANIFEST="${2}"
    OUTPUT_QZA="${3}"

    echo "sample-id,absolute-filepath,direction" > ${OUTPUT_MANIFEST}

    if [[ "${INPUT_FILE}" = /* ]]; then
      ABSOLUTE_FILEPATH="${INPUT_FILE}"
    else
      ABSOLUTE_FILEPATH="$(pwd)/${INPUT_FILE}"
    fi

    echo "${SAMPLE_ID},${ABSOLUTE_FILEPATH},forward" >> ${OUTPUT_MANIFEST}

    qiime tools import \
        --type SampleData[SequencesWithQuality] \
        --input-path ${OUTPUT_MANIFEST} \
        --output-path ${OUTPUT_QZA} \
        --input-format SingleEndFastqManifestPhred33
}

log_message "Making sample manifest and importing adapter-trimmed single-end sequences to qiime2 started."

import_qiime2_manifest "${QC_DIR}/single_end_qc/tmp/adapters_trimmed/${SINGLE_END_FILENAME}" "${QC_DIR}/single_end_qc/tmp/adapters_trimmed/MANIFEST.csv" "${QC_DIR}/single_end_qc/tmp/adapters_trimmed/adapter_trimmed_seqs.qza"

log_message "Making sample manifest and importing adapter-trimmed single-end sequences to qiime2 completed."

function quality_filter {
    set -eo pipefail

    INPUT_FILE="${1}"
    OUTPUT_DIR="${2}"

    if ! qiime quality-filter q-score \
        --i-demux ${INPUT_FILE} \
        --p-min-quality 4 \
        --o-filtered-sequences ${OUTPUT_DIR}/quality_filtered_seqs.qza \
        --o-filter-stats ${OUTPUT_DIR}/quality_filter_stats.qza \
        --verbose \
        2>&1 | tee -a ${QC_DIR}/q2_quality_filter.log; then
            notify_value_error_from "${QC_DIR}/q2_quality_filter.log"
            exit 1
    fi
}

log_message "Quality filtering single-end reads with qiime2 started."

mkdir "${QC_DIR}/single_end_qc/tmp/quality_filtered"

quality_filter "${QC_DIR}/single_end_qc/tmp/adapters_trimmed/adapter_trimmed_seqs.qza" "${QC_DIR}/single_end_qc/tmp/quality_filtered"

log_message "Quality filtering single-end reads with qiime2 completed."

function export_quality_filtered_seqs {
    set -eo pipefail

    INPUT_FILE="${1}"
    OUTPUT_DIR="${2}"

    qiime tools export \
        --input-path ${INPUT_FILE} \
        --output-path ${OUTPUT_DIR} \
        --output-format SingleLanePerSampleSingleEndFastqDirFmt

    # Exported per-sample fastq files follow the format "[SAMPLEID]_[DIGITS]_L001_R1_001.fastq.gz".
    # Rename the file to be "[SAMPLEID].fastq.gz" for clarity.
    #
    # Note: the exported file extension will always be .fastq, even if the original file extension was .fq
    mv "${OUTPUT_DIR}"/*.fastq.gz "${OUTPUT_DIR}/${SINGLE_END_FILENAME}"

    # Clean up other exported files that aren't needed; we just want the per-sample fastq file.
    rm ${OUTPUT_DIR}/MANIFEST
    rm ${OUTPUT_DIR}/metadata.yml
}

log_message "Exporting quality filtered single-end sequences from qiime2 started."

export_quality_filtered_seqs "${QC_DIR}/single_end_qc/tmp/quality_filtered/quality_filtered_seqs.qza" "${QC_DIR}/single_end_qc/tmp/exported_quality_filtered_seqs"

QUALITY_FILTERED_READ_COUNT="$(count_reads "${QC_DIR}/single_end_qc/tmp/exported_quality_filtered_seqs/${SINGLE_END_FILENAME}")"
log_message "Read count after quality filtering reads: ${QUALITY_FILTERED_READ_COUNT}"

log_message "Exporting quality filtered single-end sequences from qiime2 completed."

# Deactivate qiime2 env
conda deactivate

source activate ${FASTQC_ENV}

log_message "Running fastqc on QC'd single-end sample started."

mkdir "${QC_DIR}/single_end_qc/tmp/fastqc/qc_seqs"

runfastqc "${QC_DIR}/single_end_qc/tmp/exported_quality_filtered_seqs/${SINGLE_END_FILENAME}" "${QC_DIR}/single_end_qc/tmp/fastqc/qc_seqs"

log_message "Running fastqc on QC'd single-end sample completed."

# Deactivate fastqc env
conda deactivate

# Clean up intermediate files and move final output files out of `tmp`
# directory.
mv "${QC_DIR}/single_end_qc/tmp/quality_filtered/quality_filtered_seqs.qza" "${QC_DIR}/single_end_qc/qc_seqs.qza"
mv "${QC_DIR}/single_end_qc/tmp/fastqc" "${QC_DIR}/single_end_qc"
rm -r "${QC_DIR}/single_end_qc/tmp"

# Write out read counts remaining after each QC step. Since only a single
# sample is being processed, the file will always contain a single row of data.
# This file format makes it easy to concatenate multiple read count logs (e.g.
# when aggregating results across batches of samples).
log_message "Logging single-end read counts for each QC step started."

echo -e "sample_id\traw_read_count\tprimers_trimmed_read_count\tadapters_trimmed_read_count\tquality_filtered_read_count" > "${QC_DIR}/single_end_qc/read_counts.tsv"
echo -e "${SAMPLE_ID}\t${RAW_READ_COUNT}\t${PRIMERS_TRIMMED_READ_COUNT}\t${ADAPTERS_TRIMMED_READ_COUNT}\t${QUALITY_FILTERED_READ_COUNT}" >> "${QC_DIR}/single_end_qc/read_counts.tsv"

log_message "Logging single-end read counts for each QC step completed."

PIPELINE_COMPLETED="$(date +%Y-%m-%d:%H:%M:%S)"

echo "${PIPELINE_STARTED} Pipeline execution started." > "${QC_DIR}/pipeline_completed.log"
echo "${PIPELINE_COMPLETED} Pipeline execution completed." >> "${QC_DIR}/pipeline_completed.log"
log_message "Pipeline execution completed."
