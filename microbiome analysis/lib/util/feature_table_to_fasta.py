#!/usr/bin/env python3

__author__ = "Arun Manoharan (arun@primediscoveries.com)"

import os
import os.path
import sys

try:
    import click
    import biom
    import pandas as pd
    import qiime2
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
              help="qiime2 .qza file containing feature table to convert. "
                   "Each feature ID in the table must have a corresponding "
                   "representative sequence in --rep_seqs file.")
@click.option('--rep_seqs', 'rep_seqs_path', required=True,
              type=click.Path(exists=True, file_okay=True, dir_okay=False,
                              readable=True),
              help="qiime2 .qza file containing representative sequences "
                   "corresponding to feature IDs in --feature_table file.")
@click.option('--output_dir', required=True,
              type=click.Path(exists=False, file_okay=False, dir_okay=True,
                              writable=True),
              help="Output directory in which to write per-sample fasta "
                   "files.")
def feature_table_to_fasta(feature_table_path, rep_seqs_path, output_dir):
    """Convert feature table and rep seqs to per-sample fasta files.

    Converts a feature table and representative sequences to per-sample fasta
    files. Each feature in a sample has its corresponding representative
    sequence written out N times to that sample's fasta file, where N is the
    feature's abundance within the sample. Each output fasta filename will
    correspond to a sample ID in the feature table. The output file extension
    will be .fasta (e.g. sample ID "Sample1" would have its sequences written
    to "Sample1.fasta").

    Example usage:

    \b
    python feature_table_to_fasta.py
        --feature_table feature_table.qza
        --rep_seqs rep_seqs.qza
        --output_dir per_sample_fasta

    """
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

    biom_table = feature_table.view(biom.Table)
    seqs_series = rep_seqs.view(pd.Series)
    feature_ids = biom_table.ids(axis='observation')

    for feature_id in feature_ids:
        if feature_id not in seqs_series:
            raise click.ClickException(
                "Feature ID %r in feature table does not have a corresponding "
                "representative sequence in %r." % (feature_id, rep_seqs_path))

    for counts, sample_id, _ in biom_table.iter(dense=True, axis='sample'):
        sample_path = os.path.join(output_dir, '%s.fasta' % sample_id)
        click.echo("Writing per-sample fasta file: %s" % sample_path)

        with open(sample_path, 'w') as f:
            counts = counts.astype(int)
            for feature_id, count in zip(feature_ids, counts):
                rep_seq = seqs_series.loc[feature_id]
                for _ in range(count):
                    rep_seq.write(f, format='fasta',
                                  id_whitespace_replacement=None,
                                  description_newline_replacement=None)


if __name__ == '__main__':
    feature_table_to_fasta()
