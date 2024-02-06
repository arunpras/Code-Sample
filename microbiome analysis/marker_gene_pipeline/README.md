# Marker gene pipeline quickstart guide

This quickstart guide covers everything you need to get started with analyzing
marker gene survey data. The guide will show how to perform an end-to-end
analysis, starting with raw sequence data and ending with consumer report data.

This guide doesn't cover everything you can do with the marker gene pipeline.
It covers the most important steps, considerations, and requirements. Each
script contains additional documentation in its header, so be sure to check
that out if you want to learn more about a particular step in the analysis.
The Python scripts can also be run with `--help` to view help text.

Stay tuned for a full user guide in the future!

## Pipeline organization

The pipeline is organized in the following manner:

- `single_sample`: Scripts for processing a single sample at a time,
  independent of any other samples.

- `merged_samples`: Scripts for performing comparative analyses on a set of
  samples.

- `report_data`: Scripts for generating report data for a client's sample
  compared to baseline samples.

- `run_marker_gene_pipeline.sh`: Wrapper script which executes an end-to-end
  analysis of a client's sample.

- `process_upstream.sh`: Wrapper script for processing a directory of samples
  in parallel (upstream "single-sample" analyses).

- `process_downstream.sh`: Wrapper script for processing samples in parallel
  (downstream "merged samples" and "report data" analyses).

## Pipeline data requirements

The marker gene pipeline has the current requirements:

- Data must be 16S amplicon sequences generated on an Illumina sequencer.
  **Non-Illumina data (e.g. 454, Sanger) are not currently supported.**

- FASTQ Phred quality scores must use Phred offset 33.

## Installation

The marker gene pipeline currently supports Linux and macOS.

1. Install [Miniconda3](https://docs.conda.io/en/latest/miniconda.html).

2. Run `install/create_conda_envs.sh` (script is contained in this repo) to
   install the required conda environments. Do not install, upgrade, or remove
   any packages in these environments.

3. If you intend to use `process_upstream.sh` or `process_downstream.sh` (see
   below), you'll also need the GNU `parallel` tool installed.

## Getting started with raw sequence data

**Tip:** @arunpras has a custom script for downloading raw sequence data from
NCBI or qiita. The download script organizes the data such that it is directly
usable with the pipeline, and is thus an easy way to get started. Please
contact Arun for the script and instructions.

The pipeline has the following data organization requirements:

- Each sample has its own directory, where the name of the directory is the
  sample ID used in all further analyses. The names of the sample's fastq files
  do not determine the sample ID.

- Within the sample directory, there must be a `raw_seqs` subdirectory
  containing either single-end or paired-end fastq files.

- For paired end reads, forward and reverse reads must be distinguished by `_1`
  and `_2`, or `_R1` and `_R2`, in filename suffixes. Single-end reads do not
  have this requirement.

- Reads must have file extension `.fastq` or `.fq`, and may optionally be
  gzipped (`.gz`).

## Processing a batch of samples (upstream)

The pipeline is designed to process each sample independently of other samples,
until it becomes necessary to combine samples for certain downstream analyses
(e.g. building a phylogenetic tree, picrust predictions, etc). This design is
mainly for performance and to support production data flow where a client's
sample is processed independently and compared to a set of baseline samples.

Typically you'll be analyzing more than one sample at a time though, so there
are wrapper scripts to make life easier: `process_upstream.sh` and
`process_downstream.sh`.

`process_upstream.sh` accepts a directory of samples as input, and performs a
number of initial analyses on the raw sequence data:

1. Quality control: primer trimming, adapter trimming, quality score
   trimming/filtering, and FastQC reports.
2. Denoising with Deblur to produce ASV sequences (optionally filtering bloom
   sequences).
3. Generating kmers from denoised reads.
4. Computing Shannon entropy from kmers.

`process_upstream.sh` accepts the following options:

- `--samples_dir`: directory of samples to process. Each sample directory must
  contain a `raw_seqs` subdirectory of single-end or paired-end fastq files.
  The name of each sample directory will be the sample ID used in all
  downstream analyses.

- `--deblur_trim_len`: trim length for deblur. This value will depend on the
  length of your sequences. Note that the pipeline currently only processes
  forward reads (reverse reads are supported, but ignored). A trim length of
  150 is commonly used with qiita datasets, though if your sequences are longer
  you'll typically want to use a longer trim length. If you'll be analyzing
  samples from multiple projects/batches, you'll need to use the same trim
  length in order to make the results comparable across samples (which is why a
  shorter trim length is often used in meta-analyses).

  **Tip:** With Akesogen sequence data, we've found that a trim length of 230
  for V4 data works well (because we have to trim off 19bp of primer sequence;
  see below).

- `--target_subfragment`: target subfragment of the 16S gene each sample's data
  was generated from. The default is `V4`; make sure to modify accordingly if
  you're analyzing non-V4 data. Currently this option only determines whether
  samples can be compared downstream, and a V4-specific taxonomy classifier
  will be used to provide better taxonomic classifications (otherwise a
  full-length 16S taxonomy classifier is used for non-V4 samples).

- `--trim_515f`: if provided, trim the V4 515F primer sequence from each
  forward read (i.e. the first 19 bases). A maximum of one base mismatch is
  allowed. Reads that don't start with the primer sequence, or have more than
  the allowed number of mismatches, will be discarded. **CAUTION:** only use
  this option if reads begin with the V4 515F primer sequence (this applies to
  Akesogen data); if primer sequences have already been removed, do not use
  this option. By default, no bases are trimmed.

  **Tip:** Most studies do not include the primer sequence, so the default
  behavior (no trimming) is usually what you want. Be sure to manually inspect
  the raw forward reads to see if they start with the primer. For example, if
  you're working with V4 data generated from the 515F primer region, search
  the forward reads for this regular expression: `GTG[CT]CAGC[AC]GCCGCGGTAA`.
  If most of the sequences start with this primer sequence, be sure to trim it
  off.

  **Tip:** Akesogen V4 data includes the EMP V4 515F primer at the beginning of
  each read. Be sure to use `--trim_515f` to trim it off.

- `--bloom_seqs`: path to fasta file of bloom ASV sequences to filter from ASV
  table and representative sequences after denoising. A representative sequence
  is identified as bloom if it is an exact match to a bloom sequence in the
  fasta file, or if one is a prefix of the other (i.e. if representative
  sequence and bloom sequence lengths differ). By default, no bloom filtering
  is applied. Only use this option if you're analyzing a dataset that is known
  to have bloomed (e.g. American Gut Project), or if you're performing a
  meta-analysis that includes bloomed samples.

  **Tip:** The American Gut Project published bloom sequences
  (Amir et al. 2017) are available in fasta format
  [here](https://raw.githubusercontent.com/knightlab-analyses/bloom-analyses/master/data/newbloom.all.fna).
  These bloom sequences are specific to 16S V4 data generated from human fecal
  samples. Do not blindly apply a bloom filter without being sure that your
  samples have indeed bloomed. Make sure the bloom sequences are ASV sequences
  generated from the same variable region as your samples, and that they apply
  to the type of samples you're dealing with (e.g. you may use different bloom
  sequences for human vs. animal fecal samples).

- `--jobs`: number of samples to process in parallel (default: 8).

## Processing a batch of samples (downstream)

After running `process_upstream.sh`, it's time to perform downstream analyses
with `process_downstream.sh`. Downstream analyses include:

1. Merging sample data from upstream "single-sample" analyses.
2. Taxonomy classification of ASV sequences.
3. Functional predictions with picrust.
4. Optionally generating consumer report data for each sample (treating all
   other samples as the baseline).
5. Optionally exporting data in common file formats (e.g. TSV, BIOM, FASTA) for
   use with external tools.

`process_downstream.sh` accepts the following options:

- `--samples`: file containing a list of sample directories to process (one
  sample directory path per line; using absolute paths is recommended). Each
  sample directory must contain "single-sample" results (e.g. from running
  `process_upstream.sh` or `run_marker_gene_pipeline.sh`).

- `--ref_dir`: directory containing reference files required by the pipeline.
  Quickest way to obtain the files is to clone the
  [bio-seed](https://github.com/PrimeDiscoveries/bio-seed) repo and using it
  like this: `--ref_dir /path/to/bio-seed/190501/refs`. Otherwise, follow the
  directions in `run_marker_gene_pipeline.sh` to obtain the reference files
  (important if you're processing non-V4 data, because `bio-seed` doesn't
  include those files).

- `--deblur_trim_len`: deblur trim length used in upstream processing.

- `--output_dir`: output directory in which to write results.

- `--target_subfragment`: target subfragment of the 16S gene each sample's data
  was generated from (default is V4). Currently this only affects which
  taxonomy classifier is used (V4-specific vs. full-length 16S classifier).

- `--report_data`: if provided, generate report data for each sample, treating
  all other samples as the "baseline" (default is to not generate report data).

- `--export_data`: if provided, data are exported from qiime2 `.qza` files into
  common file formats (e.g. TSV, BIOM, FASTA) for compatibility with external
  tools (default is to not export data).

- `--jobs`: number of jobs to run in parallel (default: 8).

## Processing a single sample (production use)

Use `run_marker_gene_pipeline.sh` to execute an end-to-end analysis of a single
client's sample compared to a baseline. More docs will be coming in the future
(in the meantime, check out the script's header docs).

## Output files

`TODO` (Eric is familiar with the exported data in the meantime)

## Tutorial/ example commands

`TODO`

In the meantime, check out
`/Volumes/prime/microbiome_downloads/bloom/commands.sh` on the NAS for examples
of processing 3 datasets used in the American Gut Project bloom filter paper.

Each of the pipeline scripts additionally includes usage examples in the header
docs.

## Troubleshooting/ FAQ

**Question:** Why am I getting the following error when running
`process_upstream.sh` or `run_marker_gene_pipeline.sh`?

```
Error: <sample-directory> directory must only contain a raw_seqs subdirectory.
The subdirectory either does not exist, or there are additional
files/directories present (e.g. from previous processing of this sample).
```

**Answer:** This usually happens if the sample has previously been processed
with the pipeline, or when there are unexpected files in the sample directory.
The pipeline tries its best to avoid overwriting data from previous runs to aid
in reproducibility of results.

If you're sure that you want to (re)process the sample, remove all other files
and directories from the sample directory so that only a `raw_seqs` directory
exists, then try rerunning the pipeline.

**Question:** Why am I getting the following error when running a pipeline
script?

```
activate: No such file or directory
```

**Answer:** This is an error that occurs in recent versions of `conda` related
to environment activation/deactivation. To solve it, run `conda deactivate` to
deactivate the active environment, then try rerunning the script.
