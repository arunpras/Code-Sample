#!/usr/bin/env bash
set -eo pipefail

# Authors:
#   - Arun Manoharan (arun@primediscoveries.com)
# Created January 2019 by Arun Manoharan
# Last Updated June 2019 by Arun Manoharan

# This script creates conda environments for tools/dependencies required by the pipelines.
# Each tool is pinned to the version that has been tested with the pipelines.
#
# Creates these conda environments, each with the following names:
#     - bbmap
#     - qiime2-2019.4
#     - picrust1
#     - fastqc
#     - kmer-pipeline
#
# Any existing environments matching the names above will be skipped and will not be modified.

# Program requirements:
#     Miniconda3, Miniconda2, or Anaconda installation: https://conda.io/miniconda.html

# Example usage:
#     bash create_conda_envs.sh

if (( $# != 0 )); then
    echo "This script does not accept any arguments. Please see script header for usage instructions."
    exit 1
fi

# bbmap 38.22
# https://bioconda.github.io/recipes/bbmap/README.html
if (( $(conda env list | grep -E "^bbmap " | wc -l) == 1)); then
    echo "conda environment 'bbmap' already exists, skipping environment creation."
else
    echo "Creating conda environment 'bbmap':"
    conda create --yes --no-default-packages --override-channels -c conda-forge -c bioconda -c defaults -n bbmap bbmap=38.22
fi

# qiime2-2019.4
# https://docs.qiime2.org/2019.4/install/native/
if (( $(conda env list | grep -E "^qiime2-2019\.4 " | wc -l) == 1)); then
    echo "conda environment 'qiime2-2019.4' already exists, skipping environment creation."
else
    # https://stackoverflow.com/a/17072017/3776794
    if [ "$(uname)" == "Darwin" ]; then
        echo "Detected system platform: macOS"
        URL="https://data.qiime2.org/distro/core/qiime2-2019.4-py36-osx-conda.yml"
        FILENAME="qiime2-2019.4-py36-osx-conda.yml"
    elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
        echo "Detected system platform: Linux"
        URL="https://data.qiime2.org/distro/core/qiime2-2019.4-py36-linux-conda.yml"
        FILENAME="qiime2-2019.4-py36-linux-conda.yml"
    else
        echo "Couldn't detect whether system platform is macOS or Linux while attempting to create qiime2-2019.4 environment."
        exit 1
    fi

    echo "Creating conda environment 'qiime2-2019.4':"
    curl -sL "${URL}" > "${FILENAME}"
    # `conda env create` doesn't have a --yes option currently, so use --force instead.
    conda env create --force -n qiime2-2019.4 --file "${FILENAME}"
    rm "${FILENAME}"
fi

# picrust 1.1.3
# http://picrust.github.io/picrust/install.html
if (( $(conda env list | grep -E "^picrust1 " | wc -l) == 1)); then
    echo "conda environment 'picrust1' already exists, skipping environment creation."
else
    echo "Creating conda environment 'picrust1':"
    conda create --yes --no-default-packages --override-channels -c conda-forge -c bioconda -c defaults -n picrust1 picrust=1.1.3
    source activate picrust1
    echo "Downloading picrust1 reference data:"
    # Downloads reference data and stores it within the picrust1 conda
    # environment.
    download_picrust_files.py -t ko -g 13_5
    conda deactivate
fi

# fastqc 0.11.8
# https://bioconda.github.io/recipes/fastqc/README.html
if (( $(conda env list | grep -E "^fastqc " | wc -l) == 1)); then
    echo "conda environment 'fastqc' already exists, skipping environment creation."
else
    echo "Creating conda environment 'fastqc':"
    conda create --yes --no-default-packages --override-channels -c conda-forge -c bioconda -c defaults -n fastqc fastqc=0.11.8
fi

# kmer pipeline (code lives within this repository)
if (( $(conda env list | grep -E "^kmer-pipeline " | wc -l) == 1)); then
    echo "conda environment 'kmer-pipeline' already exists, skipping environment creation."
else
    echo "Creating conda environment 'kmer-pipeline':"

    # Install the kmer pipeline's dependencies via conda for better
    # cross-platform compatibility and to avoid having to reference an external
    # requirements file. The pinned versions here match the pinned versions in
    # lib/kmer_pipeline/requirements.txt. If the pinned versions change in the
    # future, please be sure to update this file *and*
    # lib/kmer_pipeline/requirements.txt.
    #
    # Note: we're only installing packages required to generate kmers
    # (`generate_kmers.py`) to avoid installing unnecessary dependencies. See
    # `lib/kmer_pipeline/requirements.txt` for the full list of dependencies
    # required by other kmer pipeline scripts.
    conda create --yes --no-default-packages --override-channels -c conda-forge -c bioconda -c defaults -n kmer-pipeline \
        python=3.6 \
        numpy=1.14.3 \
        biom-format=2.1.6 \
        scipy=1.0.0 \
        scikit-learn=0.19.1 \
        tqdm=4.19.5 \
        biopython=1.71
fi
