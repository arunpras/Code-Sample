#!/usr/bin/env python3

__author__ = "Arun Manoharan (arun@primediscoveries.com)"

import csv
import datetime
import glob
import os.path
import re
import sys

try:
    import click
    import numpy as np
    import scipy.sparse
    import scipy.stats
except ImportError as e:
    sys.exit(
        "Couldn't import required dependencies. Please ensure that a qiime2 "
        "conda environment is activated prior to running this script (use "
        "install/create_conda_envs.sh to create one).\n\n"
        "Original error message: %s" % e)


@click.command()
@click.option('--kmer_dir', required=True,
              type=click.Path(exists=True, file_okay=False, dir_okay=True,
                              readable=True, writable=True),
              help="Directory containing single-sample kmer results for which "
                   "to compute entropy. This directory should be the output "
                   "from running generate_kmers.sh.")
def compute_kmer_entropy(kmer_dir):
    """Compute Shannon entropy for a single sample's kmers.

    Computes Shannon entropy (using natural logarithm base e) for a single
    sample's kmers. The sample must first be processed with generate_kmers.sh
    to generate the necessary kmer results. Entropy is computed for each set of
    kmers present (i.e. kmer length and subsampling size), and results are
    written in TSV format.

    Example usage:

    \b
    python compute_kmer_entropy.py
        --kmer_dir SAMPLEID/asvs/deblur_single_end_trim_len_150/kmers

    Outputs:

    <KMER_DIR>/entropy: directory containing entropies stored in TSV format

    """
    output_dir = os.path.join(kmer_dir, 'entropy')

    if os.path.exists(output_dir):
        raise click.ClickException(
            "%r directory already exists. Please either remove the directory "
            "and rerun this script or specify a different --kmer_dir." %
            output_dir)

    os.makedirs(output_dir)

    log = os.path.join(output_dir, 'compute_kmer_entropy.log')

    # Save pipeline start timestamp for writing to pipeline_completed.log file
    # at end of script.
    pipeline_started = datetime.datetime.now().strftime('%Y-%m-%d:%H:%M:%S')
    log_message("Pipeline execution started.", log)

    # Export a copy of input parameters.
    log_message("Writing file of input parameters started.", log)

    params_log = os.path.join(output_dir, "compute_kmer_entropy_params.log")
    with open(params_log, 'w') as f:
        f.write("KMER_DIR: %s\n" % kmer_dir)

    log_message("Writing file of input parameters completed.", log)

    log_message("Loading single-sample kmer results started.", log)

    sample_id = get_sample_id(kmer_dir)
    kmers = load_single_sample_kmers(kmer_dir, sample_id)

    log_message("Loading single-sample kmer results completed.", log)

    if len(kmers) == 0:
        log_message("Did not find kmer results for sample %r. Entropy will "
                    "not be computed for this sample." % sample_id, log)
    else:
        log_message("Computing kmer entropy started.", log)

        kmer_entropy = _compute_kmer_entropy(kmers)

        log_message("Computing kmer entropy completed.", log)

        log_message("Saving kmer entropy results started.", log)

        save_kmer_entropy_results(kmer_entropy, sample_id, output_dir)

        log_message("Saving kmer entropy results completed.", log)

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
            f.write("compute_kmer_entropy.py log file\n\n")

    msg = '%s %s' % (datetime.datetime.now().strftime('%Y-%m-%d:%H:%M:%S'),
                     msg)
    click.echo(msg)
    with open(filepath, 'a') as f:
        f.write('%s\n' % msg)


def get_sample_id(kmer_dir):
    # normpath removes trailing slashes from directory path. This is
    # necessary prior to using dirname and basename.
    current_dir = os.path.normpath(kmer_dir)
    for _ in range(3):
        current_dir = os.path.dirname(current_dir)
    return os.path.basename(current_dir)


# TODO `load_single_sample_kmers` and `load_kmer_matrix` are mostly duplicated
# from merge_samples.py. Consider centralizing kmer I/O code in `util`
# directory in the future.

def load_single_sample_kmers(kmer_dir, sample_id):
    kmers = {}
    for npz_path in glob.glob(os.path.join(kmer_dir, '*.npz')):
        # .npz filename should match pattern:
        #
        # <sample-id>_<kmer-len>-mers_<subsample-size>.npz
        #
        # Escape the sample ID in case it has special regex characters.
        # Also, <subsample-size> can be -1.
        regex = r'^%s_(\d+)-mers_(-?\d+)\.npz$' % re.escape(sample_id)
        match = re.match(regex, os.path.basename(npz_path))

        if match is None:
            raise click.ClickException(
                "%r does not match the expected kmer pipeline output "
                "filename format. Expected the filename to match the "
                "following pattern:\n\n"
                "%s_<kmer-len>-mers_<subsample-size>.npz\n\n"
                "Please run generate_kmers.sh to generate kmers for this "
                "sample in the expected format." % (npz_path, sample_id))

        kmer_len = int(match.group(1))
        subsample_size = int(match.group(2))

        sparse_matrix = load_kmer_matrix(npz_path, sample_id, kmer_len)
        kmers[(kmer_len, subsample_size)] = sparse_matrix

    return kmers


def load_kmer_matrix(npz_path, sample_id, kmer_len):
    """Load a single-sample kmer matrix from .npz/_meta files."""
    meta_path = os.path.splitext(npz_path)[0] + '_meta'

    if not os.path.exists(meta_path):
        raise click.ClickException(
            "%r file does not exist. %r must have a corresponding _meta file "
            "produced by the kmer pipeline." % (meta_path, npz_path))

    with open(meta_path, 'r') as fh:
        lines = fh.readlines()

    if len(lines) != 1:
        raise click.ClickException(
            "%r file must have exactly one line, but found %d lines. kmer "
            "pipeline results are expected to be single-sample."
            % (meta_path, len(lines)))

    line = lines[0].rstrip('\n')
    if os.path.basename(line) != '%s.fasta' % sample_id:
        raise click.ClickException(
            "%r file must contain a fasta file path for this sample: %s.fasta"
            % (meta_path, sample_id))

    with np.load(npz_path) as npz_fh:
        data = npz_fh['data']
        indices = npz_fh['indices']
        indptr = npz_fh['indptr']
        shape = npz_fh['shape']

    expected_shape = np.array([1, 4**kmer_len])
    if not np.array_equal(shape, expected_shape):
        raise click.ClickException(
            "%r file contains matrix with shape %r. Expected shape to be %r, "
            "which corresponds to %d-mers for a single sample."
            % (npz_path, shape.tolist(), expected_shape.tolist(), kmer_len))

    return scipy.sparse.csr_matrix((data, indices, indptr), shape=shape)


def _compute_kmer_entropy(kmers):
    kmer_entropy = {}
    for (kmer_len, subsample_size), sparse_matrix in kmers.items():
        # Convert to dense array and squeeze to change shape from (1, N) to
        # (N,).
        dense_array = np.squeeze(sparse_matrix.toarray())
        entropy = scipy.stats.entropy(dense_array)
        kmer_entropy[(kmer_len, subsample_size)] = entropy
    return kmer_entropy


def save_kmer_entropy_results(kmer_entropy, sample_id, output_dir):
    header = ['sample-id', 'shannon-entropy']

    for (kmer_len, subsample_size), entropy in kmer_entropy.items():
        tsv_path = os.path.join(output_dir,
                                '%s_%d-mers_%d_entropy.tsv' %
                                (sample_id, kmer_len, subsample_size))

        # Newline settings based on recommendation from csv docs:
        #     https://docs.python.org/3/library/csv.html#id3
        with open(tsv_path, 'w', newline='', encoding='utf-8') as fh:
            tsv_writer = csv.writer(fh, dialect='excel-tab', strict=True)
            tsv_writer.writerow(header)

            # Format floats with *up to* 15 digits of total precision.
            # Scientific notation may be used if appropriate. Integers won't
            # have trailing zeros included (e.g. 0.0 will be written as "0";
            # 42.0 will be written as "42", etc).
            row = [sample_id, '{0:.15g}'.format(entropy)]
            tsv_writer.writerow(row)


if __name__ == '__main__':
    compute_kmer_entropy()
