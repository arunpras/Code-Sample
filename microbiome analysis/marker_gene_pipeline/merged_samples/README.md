# `merged_samples` directory

This directory contains scripts for performing comparative analyses on a set of samples.

## Overview

- `merge_samples.py`: Script for combining single-sample results for further comparative analyses. Use this script before executing the other scripts below.

- `classify_taxonomy.sh`: Performs taxonomic classification of ASV representative sequences using qiime2's naive bayes feature classifier.

- `picrust_predictions.sh`: Generates picrust1 metagenome function prediction tables from an ASV feature table and representative sequences.
