# `single_sample` directory

This directory contains scripts for processing a single sample at a time, independent of any other samples. Begin processing your sample's raw sequence data here.

## Overview

- `quality_control.sh`: Performs QC on raw single-end or paired-end fastq files for a single sample.

- `generate_asvs.sh`: Generates an ASV feature table and representative sequences for a single sample.

- `generate_kmers.sh`: Generates kmers from a single sample's ASV sequences.

- `compute_kmer_entropy.py`: Computes Shannon entropy for a single sample's kmers.
