#!/usr/bin/env python3

__author__ = "Arun Manoharan (arun@primediscoveries.com)"

import csv
import os
import os.path
import sys

try:
    import biom
    import click
    import pandas as pd
    import qiime2
    from qiime2.plugins.feature_table.actions import (
        filter_features, filter_seqs)
    import skbio
    import skbio.io
except ImportError as e:
    sys.exit(
        "Couldn't import required dependencies. Please ensure that a qiime2 "
        "conda environment is activated prior to running this script (use "
        "install/create_conda_envs.sh to create one).\n\n"
        "Original error message: %s" % e)


@click.command()
@click.option('--feature_table', 'feature_table_path', required=True,
              type=click.Path(exists=True, file_okay=True, dir_okay=False,
                              readable=True),
              help="qiime2 .qza file containing feature table to filter. "
                   "Each feature ID in the table must have a corresponding "
                   "representative sequence in --rep_seqs file.")
@click.option('--rep_seqs', 'rep_seqs_path', required=True,
              type=click.Path(exists=True, file_okay=True, dir_okay=False,
                              readable=True),
              help="qiime2 .qza file containing representative sequences to "
                   "filter.")
@click.option('--bloom_seqs', 'bloom_seqs_path', required=True,
              type=click.Path(exists=True, file_okay=True, dir_okay=False,
                              readable=True),
              help="fasta file of bloom sequences to filter from feature "
                   "table and representative sequences. A representative "
                   "sequence is identified as bloom if it is an exact match "
                   "to a bloom sequence, or if one is a prefix of the other "
                   "(i.e. if representative sequence and bloom sequence "
                   "lengths differ).")
@click.option('--output_dir', required=True,
              type=click.Path(exists=False, file_okay=False, dir_okay=True,
                              writable=True),
              help="Output directory in which to write bloom-filtered data.")
def filter_blooms(feature_table_path, rep_seqs_path, bloom_seqs_path,
                  output_dir):
    """Filter bloom sequences from a feature table and representative seqs."""
    if os.path.exists(output_dir):
        raise click.ClickException(
            "%r directory already exists. Please either remove the directory "
            "and rerun this script or specify a different --output_dir." %
            output_dir)

    os.makedirs(output_dir)

    feature_table = qiime2.Artifact.load(feature_table_path)

    if str(feature_table.type) != 'FeatureTable[Frequency]':
        raise click.ClickException(
            "%r is type %s. Please provide --feature_table a qiime2 artifact "
            "of type FeatureTable[Frequency]." %
            (feature_table_path, feature_table.type))

    rep_seqs = qiime2.Artifact.load(rep_seqs_path)

    if str(rep_seqs.type) != 'FeatureData[Sequence]':
        raise click.ClickException(
            "%r is type %s. Please provide --rep_seqs a qiime2 artifact of "
            "type FeatureData[Sequence]." % (rep_seqs_path, rep_seqs.type))

    bloom_seqs = set()
    for seq in skbio.io.read(bloom_seqs_path, format='fasta',
                             constructor=skbio.DNA):
        seq = str(seq)
        if seq in bloom_seqs:
            raise click.ClickException(
                "Found duplicate bloom sequence in %r: %r" %
                (bloom_seqs_path, seq))
        bloom_seqs.add(seq)

    if len(bloom_seqs) < 1:
        raise click.ClickException(
            "--bloom_seqs file must contain at least one sequence.")

    bloom_ids = find_bloom_ids(rep_seqs, bloom_seqs)

    if len(bloom_ids) == 0:
        # Short-circuit if no IDs are to be filtered (otherwise an error is
        # raised by `qiime2.Metadata`).
        filtered_table = feature_table
        filtered_seqs = rep_seqs
    else:
        metadata = qiime2.Metadata(
            pd.DataFrame([], index=pd.Index(bloom_ids, name='feature-id')))

        filtered_table, = filter_features(
            table=feature_table, metadata=metadata, exclude_ids=True)
        filtered_seqs, = filter_seqs(
            data=rep_seqs, metadata=metadata, exclude_ids=True)

    filtered_table.save(os.path.join(output_dir, 'feature_table.qza'))
    filtered_seqs.save(os.path.join(output_dir, 'rep_seqs.qza'))

    # Collect stats about the filtering process.
    table_before = feature_table.view(biom.Table)
    num_features_before = table_before.length(axis='observation')
    read_count_before = table_before.sum()

    table_after = filtered_table.view(biom.Table)
    num_features_after = table_after.length(axis='observation')
    read_count_after = table_after.sum()

    bloom_read_count = read_count_before - read_count_after
    num_unique_bloom_seqs = num_features_before - num_features_after

    stats_path = os.path.join(output_dir, 'stats.tsv')
    with open(stats_path, 'w', newline='', encoding='utf-8') as fh:
        tsv_writer = csv.writer(fh, dialect='excel-tab', strict=True)
        tsv_writer.writerow(['bloom_read_count', 'num_unique_bloom_seqs'])
        tsv_writer.writerow(['%d' % bloom_read_count,
                             '%d' % num_unique_bloom_seqs])


def find_bloom_ids(rep_seqs, bloom_seqs):
    rep_seqs = rep_seqs.view(pd.Series)

    # TODO if performance becomes an issue, this could be rewritten in a
    # smarter way, such as using a prefix tree to identify matches. The current
    # approach seems fine because the list of bloom sequences should always be
    # pretty small (e.g. American Gut only published 20 bloom sequences).
    bloom_ids = set()
    for id_, rep_seq in rep_seqs.iteritems():
        rep_seq = str(rep_seq)

        # Check first if there's an exact match using set lookup (should be
        # quick).
        if rep_seq in bloom_seqs:
            bloom_ids.add(id_)
        else:
            # If not, check each bloom sequence to see if it is a prefix of the
            # rep seq (or if the rep seq is a prefix of it).
            for bloom_seq in bloom_seqs:
                if rep_seq.startswith(bloom_seq) or \
                        bloom_seq.startswith(rep_seq):
                    bloom_ids.add(id_)
                    break

    # Present a more informative error message than the one raised in qiime2
    # filtering methods (this is the case where all data are filtered, which
    # causes issues downstream).
    if len(bloom_ids) == len(rep_seqs):
        raise click.ClickException(
            "All representative sequences were identified as bloom sequences.")
    return bloom_ids


if __name__ == '__main__':
    filter_blooms()
