#!/usr/bin/env python3

__author__ = "Arun Manoharan (arun@primediscoveries.com)"

# Note: these probiotic score calculations are based on @ericproffitt's
# compute_probiotic_score.py script in PR #39

import collections
import csv
import datetime
import os.path
import sys

try:
    import click
    import pandas as pd
    import qiime2
except ImportError as e:
    sys.exit(
        "Couldn't import required dependencies. Please ensure that a qiime2 "
        "conda environment is activated prior to running this script (use "
        "install/create_conda_envs.sh to create one).\n\n"
        "Original error message: %s" % e)


@click.command()
@click.option('--taxonomy_dir', required=True,
              type=click.Path(exists=True, file_okay=False, dir_okay=True,
                              readable=True),
              help="Directory containing taxonomy classification results "
                   "generated by classify_taxonomy.sh.")
@click.option('--sample_id', required=True,
              help="Sample ID for which to compute probiotic scores. This "
                   "sample ID must exist in --taxonomy_dir. The sample's "
                   "probiotic taxa abundances will be compared to probiotic "
                   "taxa abundances in all other samples in the file in order "
                   "to produce probiotic scores.")
@click.option('--probiotic_taxa', 'probiotic_taxa_path', required=True,
              type=click.Path(exists=True, file_okay=True, dir_okay=False,
                              readable=True),
              help="File containing list of probiotic taxa to compute scores "
                   "for. Each line should contain a single probiotic genus "
                   "name in the form g__Lactobacillus. The names must match "
                   "Greengenes genus names (other databases are not supported "
                   "at this time).")
@click.option('--output_dir', required=True,
              type=click.Path(exists=False, file_okay=False, dir_okay=True,
                              writable=True),
              help="Output directory in which to write the sample's probiotic "
                   "scores in TSV format.")
@click.option('--underflow_threshold', required=False, type=float,
              default=0.0, show_default=True,
              help="Threshold below which a sample's probiotic abundance "
                   "should be underflowed to zero. The underflow is applied "
                   "to all samples in --taxonomy_dir (i.e. baseline and "
                   "client samples). By default, no underflow threshold is "
                   "applied. If supplied, should be a positive number "
                   "(e.g. 1e-5).")
def compute_probiotic_scores(taxonomy_dir, sample_id, probiotic_taxa_path,
                             output_dir, underflow_threshold):
    """Compute probiotic scores by comparing sample's abundances to baseline.

    Computes probiotic taxa scores by comparing a sample's probiotic abundances
    to a set of baseline probiotic abundances. A probiotic taxa score is
    computed for each taxon in --probiotic_taxa, as well as an overall score
    across all probiotic taxa. The sample's abundance is ranked by computing
    its percentile (presented as a number between 0.0 and 1.0). The sample's
    abundance is also classified as "low", "average", or "high" based on the
    baseline's lower quartile, inner two quartiles, and upper quartile,
    respectively.

    Example usage:

    \b
    python compute_probiotic_scores.py
        --taxonomy_dir merged_samples/taxonomy/gg-13-8-99-nb-classifier
        --sample_id SAMPLEID
        --probiotic_taxa data/probiotic_taxa.txt
        --output_dir probiotic_scores

    """
    if underflow_threshold < 0.0:
        raise click.ClickException(
            "%r is an invalid value for --underflow_threshold. Please supply "
            "a value greater than or equal to zero." % underflow_threshold)

    if os.path.exists(output_dir):
        raise click.ClickException(
            "%r directory already exists. Please either remove the directory "
            "and rerun this script or specify a different --output_dir." %
            output_dir)

    os.makedirs(output_dir)

    log = os.path.join(output_dir, 'compute_probiotic_scores.log')

    pipeline_started = datetime.datetime.now().strftime('%Y-%m-%d:%H:%M:%S')
    log_message("Pipeline execution started.", log)

    log_message("Writing file of input parameters started.", log)

    params_log = os.path.join(output_dir,
                              "compute_probiotic_scores_params.log")
    with open(params_log, 'w') as f:
        f.write("TAXONOMY_DIR: %s\n" % taxonomy_dir)
        f.write("SAMPLE_ID: %s\n" % sample_id)
        f.write("PROBIOTIC_TAXA: %s\n" % probiotic_taxa_path)
        f.write("OUTPUT_DIR: %s\n" % output_dir)
        f.write("UNDERFLOW_THRESHOLD: %s\n" % underflow_threshold)

    log_message("Writing file of input parameters completed.", log)

    log_message("Loading probiotic taxa file started.", log)

    probiotic_taxa = load_probiotic_taxa_file(probiotic_taxa_path)

    log_message("Loading probiotic taxa file completed.", log)

    table_path = os.path.join(taxonomy_dir, 'level_6_genus',
                              'genus_level_table.qza')

    log_message("Loading genus level taxonomy table %r started." % table_path,
                log)

    if not os.path.exists(table_path):
        raise click.ClickException(
            "Failed to locate genus level taxonomy table %r. Please run "
            "classify_taxonomy.sh with a Greengenes classifier to generate "
            "these results." % table_path)

    table_artifact = qiime2.Artifact.load(table_path)
    table_df = table_artifact.view(pd.DataFrame)

    if sample_id not in table_df.index:
        raise click.ClickException(
            "Sample ID %r is not present in genus level taxonomy table %r." %
            (sample_id, table_path))

    if len(table_df.index) < 2:
        raise click.ClickException(
            "There must be at least one baseline sample to compare to "
            "sample %r." % sample_id)

    log_message("There are %d baseline samples to compare to sample %r." %
                (len(table_df.index) - 1, sample_id), log)

    log_message(
        "Loading genus level taxonomy table %r completed." % table_path, log)

    log_message("Computing probiotic scores started.", log)

    scores = _compute_probiotic_scores(table_df, sample_id, probiotic_taxa,
                                       underflow_threshold)

    log_message("Computing probiotic scores completed.", log)

    log_message("Saving probiotic scores started.", log)

    save_probiotic_scores(scores, output_dir)

    log_message("Saving probiotic scores completed.", log)

    pipeline_completed = datetime.datetime.now().strftime('%Y-%m-%d:%H:%M:%S')
    with open(os.path.join(output_dir, 'pipeline_completed.log'), 'w') as f:
        f.write('%s Pipeline execution started.\n' % pipeline_started)
        f.write('%s Pipeline execution completed.\n' % pipeline_completed)
    log_message("Pipeline execution completed.", log)


def log_message(msg, filepath):
    """Write message with timestamp to log file and also display on stdout.

    Will create log file if it does not already exist.

    """
    if not os.path.exists(filepath):
        with open(filepath, 'w') as f:
            f.write("compute_probiotic_scores.py log file\n\n")

    msg = '%s %s' % (datetime.datetime.now().strftime('%Y-%m-%d:%H:%M:%S'),
                     msg)
    click.echo(msg)
    with open(filepath, 'a') as f:
        f.write('%s\n' % msg)


def load_probiotic_taxa_file(path):
    # Storing as list to maintain order for writing output file. There are list
    # lookups but the list should be small enough (e.g. we're testing with 6
    # elements) that it isn't worth it to use a set.
    probiotic_taxa = []

    with open(path, 'r') as f:
        for line in f:
            taxon = line.strip()

            if not (taxon.startswith('g__') and len(taxon[3:]) > 0):
                raise click.ClickException(
                    "Invalid genus name %r in probiotic taxa file %r. Each "
                    "name must start with 'g__' and correspond to a genus "
                    "name in Greengenes." % (taxon, path))

            if taxon in probiotic_taxa:
                raise click.ClickException(
                    "Duplicate genus name %r in probiotic taxa file %r."
                    % (taxon, path))

            probiotic_taxa.append(taxon)

    if len(probiotic_taxa) == 0:
        raise click.ClickException(
            "Probiotic taxa file %r must contain at least one genus." % path)

    return probiotic_taxa


def _compute_probiotic_scores(table_df, sample_id, probiotic_taxa,
                              underflow_threshold):
    # Drop non-probiotic columns.
    cols_to_drop = []
    for column in table_df.columns:
        genus = column.split(';')[-1].strip()
        if genus not in probiotic_taxa:
            cols_to_drop.append(column)

    table_df = table_df.drop(columns=cols_to_drop)

    # Rename columns to only include genus name with g__ prefix removed (for
    # user display).
    new_col_names = []
    for column in table_df.columns:
        genus = column.split(';')[-1].strip()[3:]
        if genus in new_col_names:
            raise click.ClickException(
                "Probiotic genus %r appears multiple times in genus level "
                "table (the genus name has more than one taxonomic hierarchy "
                "leading to it)." % genus)
        new_col_names.append(genus)

    table_df.columns = new_col_names

    # Ensure all probiotic taxa were found in the genus level table.
    probiotic_taxa = [t[3:] for t in probiotic_taxa]
    missing_taxa = set(probiotic_taxa) - set(table_df.columns)
    if missing_taxa:
        raise click.ClickException(
            "The following probiotic taxa are not present in any of the "
            "samples: %s" % ', '.join(repr(t) for t in sorted(missing_taxa)))

    # Sort the columns in the order provided by the probiotic taxa file.
    table_df = table_df[probiotic_taxa]

    # Apply underflow threshold.
    table_df = table_df * (table_df >= underflow_threshold)

    if len(table_df.columns) > 1:
        overall_column = 'overall'
        if overall_column in table_df.columns:
            # This should never actually happen in practice.
            raise click.ClickException(
                "Column name %r conflicts with a name reserved for computing "
                "an overall score. Please use a different taxon name in the "
                "probiotic taxa file." % overall_column)
        # Insert "overall" column at beginning of dataframe because we want
        # this column listed first in the output scores.
        table_df.insert(0, overall_column, table_df.sum(axis=1))

    # Use an OrderedDict to maintain order of probiotic taxa for writing output
    # file.
    scores = collections.OrderedDict()
    for column in table_df.columns:
        series = table_df[column]

        # Series.pop removes item in-place, so `series` now only contains
        # baseline abundances to compare sample abundance to.
        sample_abundance = series.pop(sample_id)

        rank = (series < sample_abundance).sum() / len(series)

        if rank < 0.25:
            group = 'low'
        elif 0.25 <= rank <= 0.75:
            group = 'average'
        else:
            group = 'high'

        scores[column] = (rank, group)

    return scores


def save_probiotic_scores(scores, output_dir):
    output_path = os.path.join(output_dir, 'scores.tsv')
    with open(output_path, 'w', newline='', encoding='utf-8') as fh:
        tsv_writer = csv.writer(fh, dialect='excel-tab', strict=True)

        header = ('key', 'value', 'class')
        tsv_writer.writerow(header)

        for taxon, (rank, group) in scores.items():
            row = (taxon, '{0:.15g}'.format(rank), group)
            tsv_writer.writerow(row)


if __name__ == '__main__':
    compute_probiotic_scores()
