#!/usr/bin/env bash
set -eo pipefail

# Authors:
#   - Arun Manoharan (arun@primediscoveries.com)
# Created November 2019 by Arun Manoharan
# Last Updated November 2019 by Arun Manoharan
#
# This script executes an end-to-end shotgun metagenomics analysis on a single
# sample, starting with raw sequence data and ending with taxonomic,
# functional, and kmer profiles.
#
# The following steps are performed:
#
# - Download raw data from NCBI SRA.
# - Perform initial QC with Sunbeam (adapter trimming,
#   quality trimming/filtering, complexity filtering, FastQC reports)
# - Perform decontamination with Sunbeam. Human and mouse host DNA is removed,
#   as well as PhiX contamination.
# - Generate kmer profiles from decontaminated reads.
# - Generate taxonomy and functional profiles from decontaminated reads using
#   humann2/metaphlan2.
#
# Requirements:
#
# - See the `bio-snap` repository for a Dockerfile that configures a suitable
#   environment for running this pipeline.
#
# Script options:
#
# --output_dir: output directory in which to write results [REQUIRED]
#
# --ref_dir: directory containing reference files required by the pipeline
#   [REQUIRED]
#
# --data_source: source from which to obtain raw sample data. Supported
#   options are 'sra' and 'local' [default: sra]
#
# --sample_id: sample accession to download and process from --data_source. If
#   --data_source is `sra`, this should be an NCBI sample accession (SAMN,
#   SAMEA, SAMD, SRS, ERS, or DRS). Required for --data_source `sra`.
#
# --sample_dir: Directory containing sample files to process. Required for
#   --data_source `local`.
#
# --kmers: type of kmers to generate from QC'd/decontaminated reads. See kmer
#   pipeline's `generate_kmers.py --KN` script option for details
#   [default: 6:-1,12:-1]
#
# --cores: maximum number of cores to use [default: 8]
#
# --keep_intermediates: keep all intermediate output files (useful for advanced
#   usage or debugging) [default: remove non-essential intermediate output]
#
# Example usage:
#
# bash run_shotgun_pipeline.sh
#     --sample_id SAMN10718830
#     --output_dir /data/SAMN10718830
#     --ref_dir /refs

function log_message {
    set -eo pipefail

    MSG="${1}"
    echo "$(date +%Y-%m-%d:%H:%M:%S) ${MSG}"
}

log_message "Pipeline execution started"

# Parse script options.
SAMPLE_ID=""
SAMPLE_DIR=""
OUTPUT_DIR=""
REF_DIR=""
DATA_SOURCE="sra"
KMERS="6:-1,12:-1"
CORES="8"
KEEP_INTERMEDIATES="NO"

# https://stackoverflow.com/a/14203146/3776794
while [[ $# -gt 0 ]]; do
    key="${1}"

    case "${key}" in
        --output_dir)
            OUTPUT_DIR="${2}"
            shift 2
            ;;
        --ref_dir)
            REF_DIR="${2}"
            shift 2
            ;;
        --data_source)
            DATA_SOURCE="${2}"
            shift 2
            ;;
        --sample_id)
            SAMPLE_ID="${2}"
            shift 2
            ;;
        --sample_dir)
            SAMPLE_DIR="${2}"
            shift 2
            ;;
        --kmers)
            KMERS="${2}"
            shift 2
            ;;
        --cores)
            CORES="${2}"
            shift 2
            ;;
        --keep_intermediates)
            KEEP_INTERMEDIATES="YES"
            shift 1
            ;;
        *)
            echo "Error: Unknown option '${key}'. Please see script header for usage instructions."
            exit 1
            ;;
    esac
done

# Make sure the required options are present.
if [[ -z "${OUTPUT_DIR}" ]]; then
    echo "Error: --output_dir is a required option."
    exit 1
fi

if [[ -z "${REF_DIR}" ]]; then
    echo "Error: --ref_dir is a required option."
    exit 1
fi

# Validate script options for samples depending on data source.
case "${DATA_SOURCE}" in
    sra)
        if [[ -z "${SAMPLE_ID}" ]]; then
            echo "Error: --sample_id is a required option for --data_source sra."
            exit 1
        fi

        case "${SAMPLE_ID}" in
            SAMN*|SAMEA*|SAMD*|SRS*|ERS*|DRS*|SRR*)
                ;;
            *)
                echo "Error: '${SAMPLE_ID}' is not a supported sample accession. Sample accessions must start with 'SAMN', 'SAMEA', 'SAMD', 'SRS', 'ERS', 'DRS', or 'SRR' followed by numbers. Please specify a different sample accession with --sample_id."
                exit 1
                ;;
        esac
        ;;
    local)
        if [[ -z "${SAMPLE_DIR}" ]]; then
            echo "Error: --sample_dir is a required option for --data_source local."
            exit 1
        fi
        FIRST_SAMPLE=$(find "${SAMPLE_DIR}" -type f | head -n 1)
        SAMPLE_PATH=$(echo -n "${FIRST_SAMPLE}" | sed -E 's/_R?[1-2].*$//' | sed -e 's/.fastq.gz//')
        OTHER_SAMPLES=$(find "${SAMPLE_DIR}" -type f -not -wholename "${SAMPLE_PATH}*" | wc -l)
        if [[ "${OTHER_SAMPLES}" -ne 0 ]]; then
            echo "Error: more samples found in directory '${SAMPLE_DIR}'. Only single sample runs are supported."
            exit 1
        fi
        SAMPLE_ID=$(basename "${SAMPLE_PATH}")
        log_message "Running in single sample local mode for sample '${SAMPLE_ID}'"
        ;;
    *)
        echo "Error: '${DATA_SOURCE}' is not a supported data source. Supported data sources are 'sra' and 'local'."
        exit 1
        ;;
esac


if [ -d "${OUTPUT_DIR}" ]; then
    echo "Error: '${OUTPUT_DIR}' directory already exists. Please either remove the directory and rerun this script or specify a different --output_dir."
    exit 1
fi

if [ ! -d "${REF_DIR}" ]; then
    echo "Error: '${REF_DIR}' directory does not exist or is not accessible. Please specify a valid reference files directory with --ref_dir."
    exit 1
fi

log_message "Creating output directory"
mkdir -p "${OUTPUT_DIR}"

log_message "Saving pipeline configuration"

PROPS_FILE="${OUTPUT_DIR}/pipeline.properties"
REPO_ROOT_DIR="$(dirname "$(cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)")"
PIPELINE_COMMIT="$(cd "${REPO_ROOT_DIR}" && git rev-parse --verify HEAD)"
REF_DIR_VERSION="$(cat "${REF_DIR}/VERSION")"

echo "# Versioning" > "${PROPS_FILE}"
echo "pipeline_commit = ${PIPELINE_COMMIT}" >> "${PROPS_FILE}"
echo "ref_dir_version = ${REF_DIR_VERSION}" >> "${PROPS_FILE}"

echo "# Script options" >> "${PROPS_FILE}"
echo "sample_id = ${SAMPLE_ID}" >> "${PROPS_FILE}"
echo "output_dir = ${OUTPUT_DIR}" >> "${PROPS_FILE}"
echo "ref_dir = ${REF_DIR}" >> "${PROPS_FILE}"
echo "data_source = ${DATA_SOURCE}" >> "${PROPS_FILE}"
echo "kmers = ${KMERS}" >> "${PROPS_FILE}"
echo "cores = ${CORES}" >> "${PROPS_FILE}"
echo "keep_intermediates = ${KEEP_INTERMEDIATES}" >> "${PROPS_FILE}"

source activate sunbeam

log_message "Initializing sunbeam project"
case "${DATA_SOURCE}" in
    sra)
        sunbeam init "${OUTPUT_DIR}" --data_acc "${SAMPLE_ID}"
        ;;
    local)
        # sunbeam sometimes has trouble detecting the proper format, trying several approaches
        sunbeam init --data_fp "${SAMPLE_DIR}" "${OUTPUT_DIR}" || \
            sunbeam init --data_fp "${SAMPLE_DIR}" --format "{sample}_{r,R?}{rp,[1-2]}{suf,_?\d*}.fastq.gz" --force "${OUTPUT_DIR}" || \
            sunbeam init --data_fp "${SAMPLE_DIR}" --format "{sample}.fastq.gz" --force "${OUTPUT_DIR}"
        ;;
esac


log_message "Configuring sunbeam project"
SUNBEAM_CONFIG="${OUTPUT_DIR}/sunbeam_config.yml"

# Sunbeam uses trimmomatic's Nextera paired-end adapters by default. Update to
# use a custom fasta file containing trimmomatic's Nextera, TruSeq2, and
# TruSeq3 single and paired-end adapters.
sunbeam config modify -i -s "qc: {adapter_fp: ${REF_DIR}/sunbeam/qc/illumina_adapters.fasta}" "${SUNBEAM_CONFIG}"
sunbeam config modify -i -s "qc: {host_fp: ${REF_DIR}/sunbeam/decontam}" "${SUNBEAM_CONFIG}"
sunbeam config modify -i -s "qc: {threads: ${CORES}}" "${SUNBEAM_CONFIG}"
# Default output directory is "sunbeam_output". Change to "pipeline_output"
# because we'll be generating output in addition to sunbeam's output.
sunbeam config modify -i -s "all: {output_fp: pipeline_output}" "${SUNBEAM_CONFIG}"

# Sunbeam treats each run accession as a separate sample, which isn't always
# the case (sometimes samples are split across multiple experiments and runs in
# NCBI).
SUNBEAM_SAMPLES="${OUTPUT_DIR}/samples.csv"
NUM_RUNS="$(wc -l "${SUNBEAM_SAMPLES}" | cut -f1 -d' ')"
if [[ "${NUM_RUNS}" != "1" ]]; then
    echo "Error: Sample accession '${SAMPLE_ID}' is composed of ${NUM_RUNS} runs across all experiments. Currently only sample accessions composed of a single run are supported. Please specify a different --sample_id."
    exit 1
fi

# Use the sample accession instead of the run accession as the sample ID in
# sunbeam analyses. We do this by replacing the first column's contents with
# the sample accession.
DOWNLOAD_PATHS="$(cut -d, -f2- "${SUNBEAM_SAMPLES}")"
echo "${SAMPLE_ID},${DOWNLOAD_PATHS}" > "${SUNBEAM_SAMPLES}"

log_message "Detecting if sample is single or paired-end data"
# Using xargs to strip leading and trailing whitespace from grepped line.
PAIRED_END="$(grep 'paired_end:' "${SUNBEAM_CONFIG}" | xargs | cut -f2 -d' ')"
if [[ "${PAIRED_END}" == "true" ]]; then
    PAIRED_END="YES"
    log_message "Sample is paired-end data"
elif [[ "${PAIRED_END}" == "false" ]]; then
    PAIRED_END="NO"
    log_message "Sample is single-end data"
else
    # This case shouldn't happen in practice. It's here in case Sunbeam changes
    # its config file in the future.
    echo "Error: Could not detect whether sample is single or paired-end data. Found this value for 'paired_end:' key in ${SUNBEAM_CONFIG}: '${PAIRED_END}'"
    exit 1
fi

log_message "Running sunbeam QC"
sunbeam run -- --configfile "${SUNBEAM_CONFIG}" --cores "${CORES}" all_qc

# Perform some basic sanity-checks on the FastQC reports (these reports were
# generated on the QC'd data, prior to decontamination). These checks can
# likely be added to or made more sophisticated in the future.
QC_DIR="${OUTPUT_DIR}/pipeline_output/qc"
FASTQC_DIR="${QC_DIR}/reports"

declare -a FASTQC_REPORTS=("${FASTQC_DIR}/${SAMPLE_ID}_1_fastqc/fastqc_data.txt")
if [[ "${PAIRED_END}" == "YES" ]]; then
    FASTQC_REPORTS+=("${FASTQC_DIR}/${SAMPLE_ID}_2_fastqc/fastqc_data.txt")
fi

for REPORT in "${FASTQC_REPORTS[@]}"; do
    log_message "Checking FastQC report for any red flags: '${REPORT}'"

    ENCODING="$(grep 'Encoding' "${REPORT}" | cut -f2)"
    if [[ "${ENCODING}" != "Sanger / Illumina 1.9" ]]; then
        echo "Error: fastq quality score encoding appears to be phred64. The pipeline currently supports phred33 encoding only. FastQC detected the following encoding: '${ENCODING}'"
        exit 1
    fi

    PER_BASE_QUAL="$(grep '>>Per base sequence quality' "${REPORT}" | cut -f2)"
    if [[ "${PER_BASE_QUAL}" != "pass" && "${PER_BASE_QUAL}" != "warn" ]]; then
        echo "Error: FastQC check for per-base sequence quality did not pass. Status: '${PER_BASE_QUAL}'"
        exit 1
    fi

    PER_SEQ_QUAL="$(grep '>>Per sequence quality scores' "${REPORT}" | cut -f2)"
    if [[ "${PER_SEQ_QUAL}" != "pass" && "${PER_SEQ_QUAL}" != "warn" ]]; then
        echo "Error: FastQC check for per-sequence quality scores did not pass. Status: '${PER_SEQ_QUAL}'"
        exit 1
    fi

    N_CONTENT="$(grep '>>Per base N content' "${REPORT}" | cut -f2)"
    if [[ "${N_CONTENT}" != "pass" ]]; then
        echo "Error: FastQC check for per-base N content did not pass. Status: '${N_CONTENT}'"
        exit 1
    fi

    ADAPTER_CONTENT="$(grep '>>Adapter Content' "${REPORT}" | cut -f2)"
    if [[ "${ADAPTER_CONTENT}" != "pass" ]]; then
        echo "Error: FastQC check for adapter content did not pass. Status: '${ADAPTER_CONTENT}'"
        exit 1
    fi
done

log_message "Running sunbeam decontam"
sunbeam run -- --configfile "${SUNBEAM_CONFIG}" --cores "${CORES}" all_decontam

conda deactivate

DECONTAM_DIR="${QC_DIR}/decontam"
SAMPLE_R1="${DECONTAM_DIR}/${SAMPLE_ID}_1.fastq.gz"

if [[ "${PAIRED_END}" == "YES" ]]; then
    SAMPLE_R2="${DECONTAM_DIR}/${SAMPLE_ID}_2.fastq.gz"
    SAMPLE_CONCAT="${DECONTAM_DIR}/${SAMPLE_ID}_R1R2_concat.fastq.gz"

    log_message "Concatenating R1 and R2 reads into single fastq file and rewriting sequence headers for humann2/metaphlan2"

    # Activate an environment that has python3 (that's all the script
    # requires).
    source activate sunbeam

    # See script docs for why this step is necessary.
    python "${REPO_ROOT_DIR}/shotgun_pipeline/util/concat_paired_end.py" \
        "${SAMPLE_R1}" "${SAMPLE_R2}" "${SAMPLE_CONCAT}"

    conda deactivate
fi

KMERS_DIR="${OUTPUT_DIR}/pipeline_output/kmers"
KMERS_SCRIPT="${REPO_ROOT_DIR}/lib/kmer_pipeline/generate_kmers.sh"

log_message "Running kmer pipeline"
# any conda environment with java would be just as good
source activate sunbeam

bash "${KMERS_SCRIPT}" \
    --inaddr "${DECONTAM_DIR}" \
    --out "${KMERS_DIR}" \
    --name "${SAMPLE_ID}" \
    --KN "${KMERS}" \
    --filetype fastq.gz \
    --cores "${CORES}"

conda deactivate

# Note: The following humann2/metaphlan2 analyses are inspired by the
# `curatedMetagenomicData_pipeline.sh` script in
# https://github.com/waldronlab/curatedMetagenomicDataHighLoad
source activate humann2

if [[ "${PAIRED_END}" == "YES" ]]; then
    HUMANN2_INPUT_FASTQ="${SAMPLE_CONCAT}"
else
    HUMANN2_INPUT_FASTQ="${SAMPLE_R1}"
fi

# TODO we may want this check to happen earlier (prior to running the kmer
# pipeline) if we determine a minimum number of reads necessary for generating
# kmers.
log_message "Checking number of reads are sufficient for humann2"
# Divide file line count by 4 to get number of fastq reads.
READ_COUNT="$(( $(gunzip -c "${HUMANN2_INPUT_FASTQ}" | wc -l | cut -f1) / 4 ))"

log_message "Sample has ${READ_COUNT} reads after QC, decontamination, and concatenating R1 and R2 reads into a single fastq file (if data are paired-end)."

# TODO make this value configurable with a script option.
# TODO this value is just a rough approximation of a reasonable cutoff for
# humann2 (humann2 needs more than 1x genome coverage, which is ~5 million
# bases). It'd be better to count the number of bases in the file (since reads
# are of different lengths at this point). 250,000 100bp reads yields ~5x
# coverage, which should usually be sufficient.
MIN_READ_COUNT="250000"
if [[ "${READ_COUNT}" -lt "${MIN_READ_COUNT}" ]]; then
    echo "Error: Sample has ${READ_COUNT} reads, which is less than the required minimum read count of ${MIN_READ_COUNT}. This sample may not have enough reads for humann2 analyses."
    exit 1
fi

HUMANN2_DIR="${OUTPUT_DIR}/pipeline_output/humann2"
NUC_DB="${REF_DIR}/humann2/chocophlan"
PROT_DB="${REF_DIR}/humann2/uniref90_diamond"

log_message "Running humann2/metaphlan2"
# TODO try out `--memory-use maximum`. This will speed up the program at the
# cost of using more memory. I can't find info about how much memory is
# required, so we may need to do some trial-and-error to see if our instances
# can handle it.
humann2 \
    --input "${HUMANN2_INPUT_FASTQ}" \
    --output "${HUMANN2_DIR}" \
    --nucleotide-database "${NUC_DB}" \
    --protein-database "${PROT_DB}" \
    --metaphlan-options "-t rel_ab --bowtie2db ${REF_DIR}/metaphlan2/mpa_v20_m200 --index v20_m200 --sample_id ${SAMPLE_ID}" \
    --output-basename "${SAMPLE_ID}" \
    --threads "${CORES}" \
    --input-format fastq.gz \
    --verbose

log_message "Normalizing humann2 tables to relative abundance"
humann2_renorm_table \
    --input "${HUMANN2_DIR}/${SAMPLE_ID}_genefamilies.tsv" \
    --output "${HUMANN2_DIR}/${SAMPLE_ID}_genefamilies_relab.tsv" \
    --units relab \
    --update-snames

humann2_renorm_table \
    --input "${HUMANN2_DIR}/${SAMPLE_ID}_pathabundance.tsv" \
    --output "${HUMANN2_DIR}/${SAMPLE_ID}_pathabundance_relab.tsv" \
    --units relab \
    --update-snames

HUMANN2_TEMP_DIR="${HUMANN2_DIR}/${SAMPLE_ID}_humann2_temp"
METAPHLAN2_DIR="${OUTPUT_DIR}/pipeline_output/metaphlan2"

mkdir "${METAPHLAN2_DIR}"

log_message "Copying and renaming important humann2 temp files"
cp \
    "${HUMANN2_TEMP_DIR}/${SAMPLE_ID}.log" \
    "${HUMANN2_DIR}/"
cp \
    "${HUMANN2_TEMP_DIR}/${SAMPLE_ID}_metaphlan_bugs_list.tsv" \
    "${METAPHLAN2_DIR}/${SAMPLE_ID}_bugs_list_relab.tsv"
cp \
    "${HUMANN2_TEMP_DIR}/${SAMPLE_ID}_metaphlan_bowtie2.txt" \
    "${METAPHLAN2_DIR}/${SAMPLE_ID}_bowtie2.txt"

# Perform some basic sanity-checks on humann2/metaphlan2 output. These checks
# can likely be added to or made more sophisticated in the future.
log_message "Checking humann2/metaphlan2 output for any red flags"

if [ ! -s "${METAPHLAN2_DIR}/${SAMPLE_ID}_bowtie2.txt" ]; then
    echo "Error: metaphlan2 bowtie2 output is empty; no species were selected during the metaphlan2 prescreen step. This may happen for a number of reasons, including but not limited to the following: 1) the sample may be non-shotgun data (e.g. marker gene data); or 2) none of the sample's reads are long enough for metaphlan2's minimum read length filter (default is 70bp)."
    exit 1
fi

# This is gross: parse out the percent unaligned sequences from humann2 log
# file and convert to an integer for bash integer comparison. Raise an error if
# >90% of the reads were unaligned.
PERCENT_UNALIGNED="$(printf "%.0f" "$(grep 'Unaligned reads after translated alignment:' "${HUMANN2_DIR}/${SAMPLE_ID}.log" | rev | cut -f2 -d' ' | rev)")"
MAX_UNALIGNED="90"
if [[ "${PERCENT_UNALIGNED}" -gt "${MAX_UNALIGNED}" ]]; then
    echo "Error: ${PERCENT_UNALIGNED}% of input reads were unaligned after translated alignment."
    exit 1
fi

# Note: The following commands are inspired by the `run_markers2.py` script in
# https://github.com/waldronlab/curatedMetagenomicDataHighLoad
log_message "Running metaphlan2 to produce marker presence and abundance tables"
metaphlan2.py \
    -t marker_pres_table \
    --bowtie2db "${REF_DIR}/metaphlan2/mpa_v20_m200" \
    --index v20_m200 \
    --sample_id "${SAMPLE_ID}" \
    --nproc "${CORES}" \
    --input_type bowtie2out \
    "${METAPHLAN2_DIR}/${SAMPLE_ID}_bowtie2.txt" \
    "${METAPHLAN2_DIR}/${SAMPLE_ID}_marker_presence.tsv"

metaphlan2.py \
    -t marker_ab_table \
    --bowtie2db "${REF_DIR}/metaphlan2/mpa_v20_m200" \
    --index v20_m200 \
    --sample_id "${SAMPLE_ID}" \
    --nproc "${CORES}" \
    --input_type bowtie2out \
    "${METAPHLAN2_DIR}/${SAMPLE_ID}_bowtie2.txt" \
    "${METAPHLAN2_DIR}/${SAMPLE_ID}_marker_abundance.tsv"

conda deactivate

if [[ "${KEEP_INTERMEDIATES}" == "NO" ]]; then
    log_message "Removing non-essential intermediate output"

    rm -rf "${OUTPUT_DIR}/pipeline_output/.snakemake"
    rm -rf "${OUTPUT_DIR}/pipeline_output/download"

    rm -r "${QC_DIR}/00_samples"
    rm -r "${QC_DIR}/01_cutadapt"
    rm -r "${QC_DIR}/02_trimmomatic"
    rm -r "${QC_DIR}/03_komplexity"
    rm -r "${QC_DIR}/cleaned"

    rm -r "${HUMANN2_TEMP_DIR}"
    # removing the decontaminated fastq files
    rm -rf "${DECONTAM_DIR}"
fi

log_message "Pipeline execution completed"
