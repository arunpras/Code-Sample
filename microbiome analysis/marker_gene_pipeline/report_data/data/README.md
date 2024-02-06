# `data` directory

This directory contains data files used by the report data scripts.

## Overview

- `probiotic_taxa.txt`: List of probiotic taxa to use with `compute_probiotic_scores.py`.

The following files are used by `compute_functional_scores.py`. KO lists were compiled by Alec Reed using KEGG Modules and [this paper](https://www.nature.com/articles/s41564-018-0337-x#MOESM3) (Supplementary Dataset 1) as references. Note that some KOs are omitted because they are not predicted by picrust1 (e.g. KOs above `K15039` and some others are missing from picrust1's Greengenes->IMG->KO mapping).

- `vitamin_KOs.tsv`: KOs associated with vitamin synthesis.

- `scfa_KOs.tsv`: KOs associated with SCFA synthesis.

- `gaba_KOs.tsv`: KOs associated with GABA synthesis.
