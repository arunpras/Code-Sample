#!/usr/bin/env bash
set -euo pipefail

# Written by Arun Manoharan (arun@primediscoveries.com)
# This script determines which ASVs are biomarkers for a user-supplied metadata class using lefse. 

# Requirements:
# A metadata file and a pyloseq object produced by otu_pipeline_script.sh
# format_lefse_input.R script must be in the same directory this script is run from

# Program Requirements:
#	lefse
#   phyloseq for loading ASV table
#	R for formatting ASV table correctly and appending metadata class information to it

# Inputs are:
    # PROJECT_PATH: absolute path to project directory *required
    # METADATA_PATH: relative path from project directory to metadata file, including file name and extension *required
    # LEFSE_SUBJECT: column name to use as subject (i.e. person or field a sample came from); must correspond to sample names in ASV table *required
    # LEFSE_CLASS: metadata column name representing factor to determine differential abundance based on *required
    # LEFSE_SUBCLASS: [Optional] stratefying factor under LEFSE_CLASS; Options are [<factor column name>, NA]
    # LEFSE_PATH: absolute path to directory containing lefse python scripts lefse-format_input.py and run_lefse.py

# Outputs are:
#	batch_visualization_data/biomarker_lefse

# Example usage:
    # bash biomarker_lefse_script.sh
    #     /home/arun/workdir/test_datasets/ncbi_paired_fastq_test
    #     /metadata/list.txt         
    #     sample_name         
    #     class_col         
    #     strat_col
    #     /home/ubuntu/miniconda2/bin/

# parse positional command line input
MIN_INPUTS=4
if [ $# -lt ${MIN_INPUTS} ]; then
    echo -e "$# command line arguements found but at least ${MIN_INPUTS} required:
    PROJECT_PATH: absolute path to project directory *required
    METADATA_PATH: relative path from project directory to metadata file, including file name and extension *required
    LEFSE_SUBJECT: column name to use as subject (i.e. person or field a sample came from); must correspond to sample names in ASV table *required
    LEFSE_CLASS: metadata column name representing factor to determine differential abundance based on *required
    LEFSE_SUBCLASS: [Optional] stratefying factor under LEFSE_CLASS; Options are [<factor column name>, NA]
    LEFSE_PATH: absolute path to directory containing lefse python scripts lefse-format_input.py and run_lefse.py

    This script requires format_lefse_input.R to be located in the same directory the script is called from.

    Example Usage:
        bash biomarker_lefse_script.sh
        /home/arun/workdir/test_datasets/ncbi_paired_fastq_test
        /metadata/list.txt         
        sample_name         
        class_col         
        strat_col
        /home/ubuntu/miniconda2/bin/"
    exit
fi

PROJECT_PATH=$1
METADATA_PATH=$2
LEFSE_SUBJECT=$3
LEFSE_CLASS=$4
LEFSE_SUBCLASS=${5:-"NA"}
LEFSE_PATH=${6:-"/home/ubuntu/miniconda2/bin/"}


####################################
# Make Directories
####################################

# make output folders if they don't exist
echo -e "\nmaking directories \n"
BATCH_DATA_DIR=${PROJECT_PATH}/batch_visualization_data
LEFSE_OUTPUT_PATH=${BATCH_DATA_DIR}/biomarker_lefse
mkdir -p ${BATCH_DATA_DIR}
mkdir -p ${LEFSE_OUTPUT_PATH}

# make log file
echo "$(date +%Y-%m-%d:%H:%M:%S) Creating directories for output completed." > ${BATCH_DATA_DIR}/biomarker_lefse_script_log.txt

# export a copy of input parameters
echo "Exporting parameters started."
echo "PROJECT_PATH: $PROJECT_PATH" > ${BATCH_DATA_DIR}/biomarker_lefse_script_parameters.txt
echo "METADATA_PATH: $METADATA_PATH" >> ${BATCH_DATA_DIR}/biomarker_lefse_script_parameters.txt
echo "LEFSE_SUBJECT: $LEFSE_SUBJECT" >> ${BATCH_DATA_DIR}/biomarker_lefse_script_parameters.txt
echo "LEFSE_CLASS: $LEFSE_CLASS" >> ${BATCH_DATA_DIR}/biomarker_lefse_script_parameters.txt
echo "LEFSE_SUBCLASS: $LEFSE_SUBCLASS" >> ${BATCH_DATA_DIR}/biomarker_lefse_script_parameters.txt
echo "Exporting parameters completed."

####################################
# Differential Abudnance with Lefse 
####################################

echo "$(date +%Y-%m-%d:%H:%M:%S) Appending metadata to feature table started." >> ${BATCH_DATA_DIR}/biomarker_lefse_script_log.txt
# Append 1-2 columns of metadata to feature table for lefse analysis
# command line arguments are
# 1) column name to use as subject (i.e. person or field a sample came from), must correspond to otu table sample names
# 2) class (i.e. what we would like differential abundance results for)
# 3) subclass (a strategying column with subgroups of class) [optional]
Rscript format_lefse_input.R ${PROJECT_PATH} ${METADATA_PATH} ${LEFSE_SUBJECT} ${LEFSE_CLASS} ${LEFSE_SUBCLASS}
echo "$(date +%Y-%m-%d:%H:%M:%S) Appending metadata to feature table completed." >> ${BATCH_DATA_DIR}/biomarker_lefse_script_log.txt

# format the biom table as lefse input
echo "$(date +%Y-%m-%d:%H:%M:%S) Calling lefse script 'lefse-format_input.py' started." >> ${BATCH_DATA_DIR}/biomarker_lefse_script_log.txt
# -f specifies features, i.e. samples and metadata info, in rows (r) or columns (c)
# -c specifies feature to use as class (i.e. metadata row/column indicating healthy or diseased)
# -s specifies feature(s) to use as subclass (stratifying subgroup of class)
# -u specifies feature to use as subject (i.e. person the sample came from)
if [ LEFSE_SUBCLASS != "NA" ]; then
    ${LEFSE_PATH}/lefse-format_input.py ${LEFSE_OUTPUT_PATH}/lefse_input_table.txt ${LEFSE_OUTPUT_PATH}/lefse_data.txt -f r -c 2 -u 1
else
    ${LEFSE_PATH}/lefse-format_input.py ${LEFSE_OUTPUT_PATH}/lefse_input_table.txt ${LEFSE_OUTPUT_PATH}/lefse_data.txt -f r -c 3 -s 2 -u 1
fi
echo "$(date +%Y-%m-%d:%H:%M:%S) Calling lefse script 'lefse-format_input.py' completed." >> ${BATCH_DATA_DIR}/biomarker_lefse_script_log.txt

# # run lefse analysis and generate results using default settings
echo "$(date +%Y-%m-%d:%H:%M:%S) Calling lefse script 'run_lefse.py' started." >> ${BATCH_DATA_DIR}/biomarker_lefse_script_log.txt
${LEFSE_PATH}/run_lefse.py ${LEFSE_OUTPUT_PATH}/lefse_data.txt ${LEFSE_OUTPUT_PATH}/lefse_results.txt
echo "$(date +%Y-%m-%d:%H:%M:%S) Calling lefse script 'run_lefse.py' completed." >> ${BATCH_DATA_DIR}/biomarker_lefse_script_log.txt
