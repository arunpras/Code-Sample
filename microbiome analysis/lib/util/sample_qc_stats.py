#!/usr/bin/env python3

__author__ = "Arun Manoharan (arun@primediscoveries.com)"

import collections
import csv
import os
import os.path
import sys

try:
    import biom
    import click
    import pandas as pd
    import qiime2
    from helpers import notify_app, MessageSeverity
except ImportError as e:
    sys.exit(
        "Couldn't import required dependencies. Please ensure that a qiime2 "
        "conda environment is activated prior to running this script (use "
        "install/create_conda_envs.sh to create one).\n\n"
        "Original error message: %s" % e)


@click.command()
@click.option('--qc_dir', required=True,
              type=click.Path(exists=True, file_okay=False, dir_okay=True,
                              readable=True),
              help="Directory containing results from quality_control.sh.")
@click.option('--asv_dir', required=True,
              type=click.Path(exists=True, file_okay=False, dir_okay=True,
                              readable=True),
              help="Directory containing results from generate_asvs.sh.")
@click.option('--output_path', required=True,
              type=click.Path(exists=False, file_okay=True, dir_okay=False,
                              writable=True),
              help="Output path to write TSV file of sample QC stats.")
def sample_qc_stats(qc_dir, asv_dir, output_path):
    """Collect QC stats about a sample.

    Collects various QC stats about a sample into a single TSV file. These
    stats may be useful to review before releasing a client sample report to
    the user, as well as for debugging sample processing issues.

    The following stats are generated for a single sample:

    `raw_reads`: number of raw reads in the sample (currently limited to
    forward reads only)

    `passed_primer_trimming`: number of reads passing primer trimming step in
    quality_control.sh (N/A if primer trimming was not performed)

    `passed_adapter_trimming`: number of reads passing adapter trimming step in
    quality_control.sh

    `passed_qc`: number of reads passing all QC filters in quality_control.sh.
    These QC'd reads are supplied to the denoising process (generate_asvs.sh)

    `phix_reads`: number of PhiX reads detected during denoising

    `chimeric_reads`: number of chimeric reads detected during denoising

    `non_16S_reads`: number of reads that are not from the 16S gene detected
    during denoising

    `bloom_reads`: number of reads identified as bloom (if a bloom filter
    was applied)

    `num_unique_bloom_seqs`: number of unique ASV sequences identifed as bloom
    (if a bloom filter was applied)

    `passed_denoising`: final number of reads passing all denoising (and
    optional bloom-filtering) filters in generate_asvs.sh. This is the total
    number of ASV reads in the sample.

    Example usage:

    \b
    python sample_qc_stats.py
        --qc_dir SAMPLEID/quality_control
        --asv_dir SAMPLEID/asvs/deblur_single_end_trim_len_150
        --output_path qc_stats.tsv

    """
    qc_stats_path = os.path.join(qc_dir, 'single_end_qc', 'read_counts.tsv')
    deblur_stats_path = os.path.join(asv_dir, 'denoised', 'deblur_stats.qza')
    bloom_stats_path = os.path.join(asv_dir, 'denoised',
                                    'filter_blooms_stats.tsv')
    asv_table_path = os.path.join(asv_dir, 'denoised', 'feature_table.qza')

    qc_stats = parse_qc_stats(qc_stats_path)
    deblur_stats = parse_deblur_stats(deblur_stats_path)
    bloom_stats = parse_bloom_stats(bloom_stats_path)
    asv_table = parse_asv_table(asv_table_path)

    sample_id = qc_stats['sample_id']
    if not (sample_id == deblur_stats.name and
            sample_id == asv_table.ids(axis='sample')[0]):
        raise click.ClickException(
            "Sample IDs do not match in %r, %r, and %r files." %
            (qc_stats_path, deblur_stats_path, asv_table_path))

    # Use an OrderedDict so that column order below is preserved in TSV output.
    qc_data = collections.OrderedDict([
        ('sample_id', sample_id),
        ('raw_reads', qc_stats['raw_read_count']),
        ('passed_primer_trimming', qc_stats['primers_trimmed_read_count']),
        ('passed_adapter_trimming', qc_stats['adapters_trimmed_read_count']),
        ('passed_qc', qc_stats['quality_filtered_read_count']),
        ('phix_reads', '%d' % deblur_stats.loc['reads-hit-artifact']),
        ('chimeric_reads', '%d' % deblur_stats.loc['reads-chimeric']),
        ('non_16S_reads', '%d' % deblur_stats.loc['reads-missed-reference']),
        ('bloom_reads', bloom_stats['bloom_read_count']),
        ('num_unique_bloom_seqs', bloom_stats['num_unique_bloom_seqs']),
        # We can sum the entire table to get total denoised read count because
        # we know there is only a single sample in the table.
        ('passed_denoising', '%d' % asv_table.sum())
    ])

    # TODO consider transposing this file format in the future (the current
    # format grows column-wise).
    with open(output_path, 'w', newline='', encoding='utf-8') as fh:
        tsv_writer = csv.DictWriter(fh, fieldnames=qc_data.keys(),
                                    dialect='excel-tab', strict=True)
        tsv_writer.writeheader()
        tsv_writer.writerow(qc_data)

    flag_sample_qc_issues(qc_data)

    fastqc_path = os.path.join(qc_dir, 'single_end_qc', 'fastqc', 'raw_seqs',
                               '%s_fastqc' % sample_id, 'fastqc_data.txt')
    flag_fastqc_issues(fastqc_path)


def parse_qc_stats(path):
    with open(path, 'r', newline='', encoding='utf-8') as fh:
        tsv_reader = csv.DictReader(fh, dialect='excel-tab', strict=True)
        rows = list(tsv_reader)

    if len(rows) != 1:
        raise click.ClickException(
            "%r file has results for %d samples. Expected the file to have "
            "results for a single sample only." % (path, len(rows)))

    return rows[0]


def parse_deblur_stats(path):
    artifact = qiime2.Artifact.load(path)

    if str(artifact.type) != 'DeblurStats':
        raise click.ClickException(
            "%r is type %s. Expected a qiime2 artifact of type DeblurStats." %
            (path, artifact.type))

    stats = artifact.view(pd.DataFrame)

    if len(stats.index) != 1:
        raise click.ClickException(
            "%r file has results for %d samples. Expected the file to have "
            "results for a single sample only." % (path, len(stats.index)))

    return stats.iloc[0]


def parse_bloom_stats(path):
    if os.path.exists(path):
        with open(path, 'r', newline='', encoding='utf-8') as fh:
            tsv_reader = csv.DictReader(fh, dialect='excel-tab', strict=True)
            rows = list(tsv_reader)

        if len(rows) != 1:
            raise click.ClickException(
                "%r file has %d rows. Expected the file to have exactly one "
                "row of data." % (path, len(rows)))

        row = rows[0]
    else:
        row = collections.defaultdict(lambda: 'N/A')

    return row


def parse_asv_table(path):
    artifact = qiime2.Artifact.load(path)

    if str(artifact.type) != 'FeatureTable[Frequency]':
        raise click.ClickException(
            "%r is type %s. Expected a qiime2 artifact of type "
            "FeatureTable[Frequency]." % (path, artifact.type))

    table = artifact.view(biom.Table)

    if table.length(axis='sample') != 1:
        raise click.ClickException(
            "%r file has results for %d samples. Expected the file to have "
            "results for a single sample only." %
            (path, table.length(axis='sample')))

    return table


# Attempt to detect and warn about potential QC issues with the sample.
# Warnings are printed to stdout and sent to the report QC app (if running in
# a production environment). Many of these checks are specific to Akesogen 16S
# V4 data and may not apply to R&D data from other sources (e.g. NCBI SRA).
# These checks may need to be modified in the future if Akesogen data changes.
# The checks are pretty basic and could become more sophisticated in the
# future. It is also possible some of the checks are too stringent. The idea is
# to warn report reviewers about potential issues (erring on the side of false
# positives).
def flag_sample_qc_issues(data):
    if float(data['raw_reads']) < 30000:
        msg = (
            "Sample only contains %s raw forward reads. Akesogen data is "
            "expected to have at least ~40,000 raw forward reads per sample."
            % data['raw_reads']
        )
        click.echo("Warning: %s" % msg)
        notify_app(msg, MessageSeverity.WARNING)

    if data['passed_primer_trimming'] == 'N/A':
        msg = (
            "Sample did not have V4 515F primers trimmed. Akesogen V4 data is "
            "expected to include primers at the beginning of reads, and these "
            "should be trimmed. Please verify whether the raw reads contain "
            "primer sequences."
        )
        click.echo("Warning: %s" % msg)
        notify_app(msg, MessageSeverity.WARNING)
    else:
        prop_passed_primer_trimming = (float(data['passed_primer_trimming']) /
                                       float(data['raw_reads']))
        if prop_passed_primer_trimming < 0.9:
            msg = (
                "Only {:.3%} of the raw reads passed the primer trimming "
                "step. Please verify that the raw reads start with V4 515F "
                "primer sequences. It is possible the data may not be V4 "
                "data, or the data doesn't start with V4 primer sequences."
            ).format(prop_passed_primer_trimming)
            click.echo("Warning: %s" % msg)
            notify_app(msg, MessageSeverity.WARNING)

    if data['phix_reads'] != '0':
        msg = (
            "%s reads may be PhiX. Akesogen data is expected to have already "
            "filtered PhiX reads from the data. If this is Akesogen data, it "
            "may be worth determining why the sample contains PhiX reads."
            % data['phix_reads']
        )
        click.echo("Warning: %s" % msg)
        notify_app(msg, MessageSeverity.WARNING)

    prop_chimeric_reads = (float(data['chimeric_reads']) /
                           float(data['passed_qc']))
    if prop_chimeric_reads > 0.05:
        msg = (
            "{:.3%} of the QC'd reads may be chimeras. This percentage of "
            "chimeric reads may be unusual."
        ).format(prop_chimeric_reads)
        click.echo("Warning: %s" % msg)
        notify_app(msg, MessageSeverity.WARNING)

    prop_non_16S_reads = (float(data['non_16S_reads']) /
                          float(data['passed_qc']))
    if prop_non_16S_reads > 0.05:
        msg = (
            "{:.3%} of the QC'd reads may not be from the 16S gene. This "
            "percentage of non-16S reads may be unusual."
        ).format(prop_non_16S_reads)
        click.echo("Warning: %s" % msg)
        notify_app(msg, MessageSeverity.WARNING)

    prop_passed_denoising = (float(data['passed_denoising']) /
                             float(data['raw_reads']))
    if prop_passed_denoising < 0.2:
        msg = (
            "Only {:.3%} of the raw reads passed QC, denoising, and bloom "
            "filtering (if a bloom filter was applied). It may be worth "
            "investigating why such a high percentage of reads were discarded."
        ).format(prop_passed_denoising)
        click.echo("Warning: %s" % msg)
        notify_app(msg, MessageSeverity.WARNING)


def flag_fastqc_issues(path):
    with open(path, 'r') as f:
        lines = f.readlines()

    encoding = parse_fastqc_metric(lines, 'Encoding')
    if encoding != "Sanger / Illumina 1.9":
        msg = (
            "fastq quality score encoding appears to be phred64. The pipeline "
            "currently supports phred33 encoding only. FastQC detected the "
            "following encoding: %r" % encoding
        )
        click.echo("Warning: %s" % msg)
        notify_app(msg, MessageSeverity.WARNING)

    seq_len = parse_fastqc_metric(lines, 'Sequence length')
    if seq_len != "251":
        msg = (
            "Raw forward reads are %sbp in length. Akesogen reads are "
            "expected to be 251bp in length. If this is Akesogen data, it may "
            "be worth investigating why the read length is different, as this "
            "may affect downstream parameters (e.g. denoiser trim length)."
            % seq_len
        )
        click.echo("Warning: %s" % msg)
        notify_app(msg, MessageSeverity.WARNING)

    per_base_qual = parse_fastqc_metric(lines, '>>Per base sequence quality')
    if per_base_qual != "pass":
        msg = (
            "FastQC check for per-base sequence quality of raw forward reads "
            "did not pass (status: %r). This may indicate an issue with read "
            "quality and it is recommended to review both raw and QC'd FastQC "
            "reports." % per_base_qual
        )
        click.echo("Warning: %s" % msg)
        notify_app(msg, MessageSeverity.WARNING)

    per_seq_qual = parse_fastqc_metric(lines, '>>Per sequence quality scores')
    if per_seq_qual != "pass":
        msg = (
            "FastQC check for per-sequence quality scores of raw forward "
            "reads did not pass (status: %r). This may indicate an issue with "
            "read quality and it is recommended to review both raw and QC'd "
            "FastQC reports." % per_seq_qual
        )
        click.echo("Warning: %s" % msg)
        notify_app(msg, MessageSeverity.WARNING)

    n_content = parse_fastqc_metric(lines, '>>Per base N content')
    if n_content != "pass":
        msg = (
            "FastQC check for per-base N content of raw forward reads did not "
            "pass (status: %r). This may indicate an issue with read quality "
            "and it is recommended to review both raw and QC'd FastQC reports."
            % n_content
        )
        click.echo("Warning: %s" % msg)
        notify_app(msg, MessageSeverity.WARNING)

    adapter_content = parse_fastqc_metric(lines, '>>Adapter Content')
    if adapter_content != "pass":
        msg = (
            "FastQC check for adapter content of raw forward reads did not "
            "pass (status: %r). It may be unusual to have adapter "
            "contamination in marker gene data, and it is recommended to "
            "review both raw and QC'd FastQC reports." % adapter_content
        )
        click.echo("Warning: %s" % msg)
        notify_app(msg, MessageSeverity.WARNING)


def parse_fastqc_metric(lines, key):
    result = None
    for line in lines:
        if line.startswith(key):
            value = line.rstrip('\n').split('\t', maxsplit=1)[1]
            if result is not None:
                raise click.ClickException(
                    "Found more than one line starting with %r in FastQC "
                    "report data." % key)
            result = value

    if result is None:
        raise click.ClickException(
            "Did not find line starting with %r in FastQC report data." % key)

    return result


if __name__ == '__main__':
    sample_qc_stats()
