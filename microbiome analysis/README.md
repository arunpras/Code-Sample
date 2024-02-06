# Bioinformatics pipelines

This repository contains bioinformatics pipelines for processing raw sequence data into end-user reports.

See the [marker gene pipeline README](marker_gene_pipeline/README.md) and the
[shotgun pipeline README](shotgun_pipeline/README.md) to get started.

## Repository overview

Pipeline scripts are organized in the following directories:

- `install`: Scripts for installing dependencies required by the pipelines.

- `marker_gene_pipeline`: Scripts for an end-to-end marker gene survey pipeline (e.g. 16S amplicon data).

- `shotgun_pipeline`: Scripts for an end-to-end shotgun metagenomics pipeline.

- `lib`: Internal library code used by the pipeline scripts. This code isn't intended to be run or interacted with directly by users.

- `legacy`: Legacy code and documentation that is not currently being used in the production pipelines, but may be useful for future reference.
