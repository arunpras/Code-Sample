#!/usr/bin/env python3

__author__ = "Arun Manoharan (arun@primediscoveries.com)"

# Note: these functional score calculations are based on @ericproffitt's
# compute_pathway_score.py script in PR #39

import collections
import csv
import datetime
import os
import os.path
import sys

try:
    import click
    import biom
    import pandas as pd
    import qiime2
    import qiime2.plugins.feature_table.actions
except ImportError as e:
    sys.exit(
        "Couldn't import required dependencies. Please ensure that a qiime2 "
        "conda environment is activated prior to running this script (use "
        "install/create_conda_envs.sh to create one).\n\n"
        "Original error message: %s" % e)


@click.command()
@click.option('--picrust_dir', required=True,
              type=click.Path(exists=True, file_okay=False, dir_okay=True,
                              readable=True),
              help="Directory containing results generated by "
                   "picrust_predictions.sh.")
@click.option('--sample_id', required=True,
              help="Sample ID for which to compute functional scores. This "
                   "sample ID must exist in --picrust_dir. The sample's "
                   "predicted functional abundances will be compared to all "
                   "other samples in the file in order to produce "
                   "functional scores.")
@click.option('--output_dir', required=True,
              type=click.Path(exists=False, file_okay=False, dir_okay=True,
                              writable=True),
              help="Output directory in which to write the sample's "
                   "functional scores in TSV format.")
@click.option('--underflow_threshold', required=False, type=float,
              default=0.0, show_default=True,
              help="Threshold below which a sample's functional abundance "
                   "should be underflowed to zero. The underflow is applied "
                   "to all samples in --picrust_dir (i.e. baseline and client "
                   "samples) after summing the abundances of KOs associated "
                   "with relevant functions. By default, no underflow "
                   "threshold is applied. If supplied, should be a positive "
                   "number (e.g. 1e-5).")
@click.option('--save_abundances', required=False, is_flag=True,
              default=False, show_default=True,
              help="If supplied, save grouped functional abundance tables in "
                   "qiime2 .qza format. This can be useful in R&D for "
                   "determining an appropriate value for "
                   "--underflow_threshold.")
def compute_functional_scores(picrust_dir, sample_id, output_dir,
                              underflow_threshold, save_abundances):
    """Compute functional scores by comparing sample's abundances to baseline.

    Computes functional scores by comparing a sample's functional abundances
    to a set of baseline functional abundances. "Functional abundances" are
    the *predicted functional potential* of a sample as computed by picrust
    (i.e. KO abundances).

    Functional scores are computed for various vitamins, SCFAs, and GABA
    synthesis by summing the abundances of KOs associated with these functions.
    Each summed functional abundance in the sample is converted to relative
    abundance and ranked against baseline abundances by computing its
    percentile (presented as a number between 0.0 and 1.0). The sample's
    abundance is also classified as "low", "average", or "high" based on the
    baseline's lower quartile, inner two quartiles, and upper quartile,
    respectively. A functional score is computed for each of the vitamins,
    SCFAs, and GABA, as well as an overall score within the vitamin and SCFA
    groups.

    Example usage:

    \b
    python compute_functional_scores.py
        --picrust_dir merged_samples/picrust_predictions
        --sample_id SAMPLEID
        --output_dir functional_scores

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

    log = os.path.join(output_dir, 'compute_functional_scores.log')

    pipeline_started = datetime.datetime.now().strftime('%Y-%m-%d:%H:%M:%S')
    log_message("Pipeline execution started.", log)

    log_message("Writing file of input parameters started.", log)

    params_log = os.path.join(output_dir,
                              "compute_functional_scores_params.log")
    with open(params_log, 'w') as f:
        f.write("PICRUST_DIR: %s\n" % picrust_dir)
        f.write("SAMPLE_ID: %s\n" % sample_id)
        f.write("OUTPUT_DIR: %s\n" % output_dir)
        f.write("UNDERFLOW_THRESHOLD: %s\n" % underflow_threshold)
        f.write("SAVE_ABUNDANCES: %s\n" % save_abundances)

    log_message("Writing file of input parameters completed.", log)

    log_message("Loading picrust1 KO metagenomes table started.", log)

    ko_path = os.path.join(picrust_dir, 'ko_metagenomes.biom')
    ko_table = load_picrust1_table(ko_path, sample_id)

    log_message("There are %d baseline samples to compare to sample %r." %
                (ko_table.length(axis='sample') - 1, sample_id), log)

    log_message("Loading picrust1 KO metagenomes table completed.", log)

    data_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                            'data')

    log_message("Computing vitamin synthesis scores started.", log)

    vitamin_path = os.path.join(data_dir, 'vitamin_KOs.tsv')
    vitamin_output_dir = os.path.join(output_dir, 'vitamins')
    os.makedirs(vitamin_output_dir)

    _compute_functional_scores(vitamin_output_dir, ko_table, vitamin_path,
                               sample_id, underflow_threshold, save_abundances)

    log_message("Computing vitamin synthesis scores completed.", log)

    log_message("Computing SCFA synthesis scores started.", log)

    scfa_path = os.path.join(data_dir, 'scfa_KOs.tsv')
    scfa_output_dir = os.path.join(output_dir, 'scfa')
    os.makedirs(scfa_output_dir)

    _compute_functional_scores(scfa_output_dir, ko_table, scfa_path, sample_id,
                               underflow_threshold, save_abundances)

    log_message("Computing SCFA synthesis scores completed.", log)

    log_message("Computing GABA synthesis scores started.", log)

    gaba_path = os.path.join(data_dir, 'gaba_KOs.tsv')
    gaba_output_dir = os.path.join(output_dir, 'gaba')
    os.makedirs(gaba_output_dir)

    _compute_functional_scores(gaba_output_dir, ko_table, gaba_path, sample_id,
                               underflow_threshold, save_abundances)

    log_message("Computing GABA synthesis scores completed.", log)

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
            f.write("compute_functional_scores.py log file\n\n")

    msg = '%s %s' % (datetime.datetime.now().strftime('%Y-%m-%d:%H:%M:%S'),
                     msg)
    click.echo(msg)
    with open(filepath, 'a') as f:
        f.write('%s\n' % msg)


def load_picrust1_table(path, sample_id):
    if not os.path.exists(path):
        raise click.ClickException(
            "Failed to locate picrust1 table %r. Please run "
            "picrust_predictions.sh to generate these results." % path)

    table = biom.load_table(path)

    if not table.exists(sample_id, axis='sample'):
        raise click.ClickException(
            "Sample ID %r is not present in picrust1 table %r."
            % (sample_id, path))

    if table.length(axis='sample') < 2:
        raise click.ClickException(
            "There must be at least one baseline sample in %r to compare to "
            "sample %r." % (path, sample_id))

    return table


def _compute_functional_scores(output_dir, table, groups_path, sample_id,
                               underflow_threshold, save_abundances):
    other_group_name = 'Other'
    id_to_groups, group_names = parse_groups_file(groups_path, table,
                                                  other_group_name)

    grouped_table = group_features_one_to_many(table, id_to_groups,
                                               group_names)
    relative_abundance_table, = \
        qiime2.plugins.feature_table.actions.relative_frequency(grouped_table)

    if save_abundances:
        relative_abundance_table.save(
            os.path.join(output_dir, 'functional_abundances.qza'))

    table_df = relative_abundance_table.view(pd.DataFrame)

    # Drop "Other" column if it exists and apply underflow threshold.
    table_df.drop(columns=[other_group_name], errors='ignore', inplace=True)
    table_df = table_df * (table_df >= underflow_threshold)

    if len(table_df.columns) > 1:
        overall_column = 'overall'
        if overall_column in table_df.columns:
            raise click.ClickException(
                "Group name %r conflicts with a name reserved for computing "
                "an overall score. Please use a different group name in the "
                "groups file." % overall_column)
        # Insert "overall" column at beginning of dataframe because we want
        # this column listed first in the output scores.
        table_df.insert(0, overall_column, table_df.sum(axis=1))

    # Use an OrderedDict to maintain input ordering of groups so that output
    # can match this order.
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

    save_functional_scores(scores, output_dir)


def parse_groups_file(groups_path, table, other_group_name):
    """Parse TSV file defining groups of features.

    Returns a dict mapping each feature ID to the set of group(s) it belongs
    to. Also returns a list of group names in the order they appear in the
    groups file.

    """
    feature_ids = set(table.ids(axis='observation'))

    id_to_groups = {}
    group_names = []

    # Newline settings based on recommendation from csv docs:
    #     https://docs.python.org/3/library/csv.html#id3
    with open(groups_path, 'r', newline='', encoding='utf-8') as fh:
        tsv_reader = (trim_trailing_empty_cells(row)
                      for row in csv.reader(fh, dialect='excel-tab',
                                            strict=True))
        for row in tsv_reader:
            # Skip blank lines.
            if len(row) == 0:
                continue

            group_name = row[0]
            ids = row[1:]

            if group_name == '':
                raise click.ClickException(
                    "Groups file %r has a group without a name. The first "
                    "cell of each line must contain a group name." %
                    groups_path)

            if group_name == other_group_name:
                raise click.ClickException(
                    "Groups file %r has a group name (%r) conflicting with a "
                    "name reserved for features that don't fall within any "
                    "other groups." % (groups_path, other_group_name))

            # `group_names` is a list, so this lookup could be slow if there's
            # a large number of groups. Consider using a set or something with
            # quicker lookup time if performance becomes an issue in the
            # future (but note that we still need to maintain the order of
            # `group_names` as they appear in the file for sorting purposes).
            if group_name in group_names:
                raise click.ClickException(
                    "Groups file %r includes a duplicate group name %r. "
                    "Please ensure group names are unique." %
                    (groups_path, group_name))

            group_names.append(group_name)

            if len(ids) == 0:
                raise click.ClickException(
                    "Groups file %r has a group %r without any features."
                    % (groups_path, group_name))

            group_exists_in_table = False
            for id_ in ids:
                if id_ == '':
                    raise click.ClickException(
                        "Groups file %r has a group %r containing a feature "
                        "ID without a name." % (groups_path, group_name))

                # picrust1 outputs all KOs that it could have predicted
                # abundances for (even if it didn't detect the KO in the
                # samples). If a KO isn't in the table, that means
                # picrust1 doesn't know about it, and we should error since
                # that KO will never appear in the results.
                if id_ not in feature_ids:
                    raise click.ClickException(
                        "Groups file %r has a group %r containing a feature "
                        "ID %r that is not a KO predicted by picrust1." %
                        (groups_path, group_name, id_))

                if id_ not in id_to_groups:
                    id_to_groups[id_] = set()

                if group_name in id_to_groups[id_]:
                    raise click.ClickException(
                        "Groups file %r has a group %r containing a duplicate "
                        "feature ID %r. Each feature ID in a group must be "
                        "unique." % (groups_path, group_name, id_))

                id_to_groups[id_].add(group_name)

                if table.data(id_, axis='observation').sum() > 0.0:
                    group_exists_in_table = True

            if not group_exists_in_table:
                raise click.ClickException(
                    "Groups file %r has a group %r where none of its feature "
                    "IDs have abundances greater than zero." %
                    (groups_path, group_name))

    # Group all remaining features into an "Other" group.
    has_other_group = False
    for id_ in feature_ids:
        if id_ not in id_to_groups:
            id_to_groups[id_] = {other_group_name}
            has_other_group = True

    if has_other_group:
        group_names.append(other_group_name)

    return id_to_groups, group_names


def trim_trailing_empty_cells(row):
    """Trim empty cells from the righthand side of the row.

    Since the groups file is a jagged TSV file, some programs (e.g. Excel) will
    export the TSV file with trailing empty cells on each row, which can be
    ignored.

    This function also trims leading and trailing whitespace from each cell.

    """
    row = [cell.strip() for cell in row]

    data_extent = -1
    for idx, cell in enumerate(row):
        if cell != '':
            data_extent = idx
    return row[:data_extent+1]


def group_features_one_to_many(biom_table, id_to_groups, group_names):
    """Group features by summing their abundances, with one-to-many support.

    Note: this function operates similarly to `qiime feature-table group`,
    with the exception that qiime2 doesn't support one-to-many relationships.

    """
    def group_function(id_, _):
        for group in id_to_groups[id_]:
            yield group, group

    # `biom.Table.collapse` requires the observations to have metadata, so fill
    # in empty metadata here.
    feature_metadata = {id_: {} for id_ in biom_table.ids(axis='observation')}
    biom_table.add_metadata(feature_metadata, axis='observation')

    # Group features by summing their abundances with support for one-to-many
    # relationships.
    grouped_table = biom_table.collapse(
        group_function, one_to_many=True, norm=False,
        include_collapsed_metadata=False, one_to_many_mode='add', strict=True,
        axis='observation')

    # Reorder the grouped features (rows) based on their order in the groups
    # file.
    sorted_grouped_table = grouped_table.sort_order(group_names,
                                                    axis='observation')

    return qiime2.Artifact.import_data('FeatureTable[Frequency]',
                                       sorted_grouped_table)


def save_functional_scores(scores, output_dir):
    if len(scores) == 1:
        output_path = os.path.join(output_dir, 'score.tsv')
        header = ('value', 'class')
    else:
        output_path = os.path.join(output_dir, 'scores.tsv')
        header = ('key', 'value', 'class')

    with open(output_path, 'w', newline='', encoding='utf-8') as fh:
        tsv_writer = csv.writer(fh, dialect='excel-tab', strict=True)
        tsv_writer.writerow(header)

        for key, (rank, group) in scores.items():
            if len(scores) == 1:
                row = ('{0:.15g}'.format(rank), group)
            else:
                row = (key, '{0:.15g}'.format(rank), group)

            tsv_writer.writerow(row)


if __name__ == '__main__':
    compute_functional_scores()