INPUT_DIR=""
OUTPUT_DIR=""
SAMPLE_NAME=""
KMERS=""
FILE_TYPE="fastq.gz"
CORES=8
KANALYZE_PATH="/opt/kanalyze-2.0.0/kanalyze.jar"

# https://stackoverflow.com/a/14203146/3776794
while [[ $# -gt 0 ]]; do
    key="${1}"

    case "${key}" in
        --inaddr)
            INPUT_DIR="${2}"
            shift 2
            ;;
        --out)
            OUTPUT_DIR="${2}"
            shift 2
            ;;
        --name)
            SAMPLE_NAME="${2}"
            shift 2
            ;;
        --KN)
            KMERS="${2}"
            shift 2
            ;;
        --filetype)
            FILE_TYPE="${2}"
            shift 2
            ;;
        --cores)
            CORES="${2}"
            shift 2
            ;;
        *)
            echo "Error: Unknown option '${key}'. Please see script header for usage instructions."
            exit 1
            ;;
    esac
done

if [[ -z "${INPUT_DIR}" ]]; then
    echo "Error: --inaddr is a required option."
    exit 1
fi
if [[ -z "${OUTPUT_DIR}" ]]; then
    echo "Error: --out is a required option."
    exit 1
fi
if [[ -z "${SAMPLE_NAME}" ]]; then
    echo "Error: --name is a required option."
    exit 1
fi
if [[ -z "${KMERS}" ]]; then
    echo "Error: --KN is a required option."
    exit 1
fi

declare -a kmer_vals=(${KMERS//,/ })
for kmer in "${kmer_vals[@]}"; do
    kmer_len=$(echo -n "$kmer" | cut -d: -f1)
    sampling=$(echo -n "$kmer" | cut -d: -f2)

    if [[ $sampling -ne "-1" ]]; then
        echo "warn: subsampling is not supported right now"
    fi

    while read -r fastq; do
        fastq_name=$(basename "${fastq}")
        output_file="${OUTPUT_DIR}/${fastq_name%%.*}_${kmer_len}-mers.tsv"
        echo "Generating ${kmer_len}-mers for ${output_file}"
        java -jar "${KANALYZE_PATH}" count -rcanonical -k $kmer_len --outfmt dec \
            -f fastqgz -o "${output_file}" -d $CORES -l $CORES -t $CORES "${fastq}"
        gzip "${output_file}"
    done < <(find "${INPUT_DIR}" -type f -name "*.${FILE_TYPE}" -print)
done
