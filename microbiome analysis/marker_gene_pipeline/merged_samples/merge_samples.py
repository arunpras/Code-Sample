#!/usr/bin/env python3

__author__ = "Arun Manoharan (arun@primediscoveries.com)"

import csv
import datetime
import itertools
import os
import os.path
import sys

try:
    import click
    import biom
    import numpy as np
    import qiime2
    import q2_types.feature_data
    import q2_types.feature_table
    import scipy.sparse
except ImportError as e:
    sys.exit(
        "Couldn't import required dependencies. Please ensure that a qiime2 "
        "conda environment is activated prior to running this script (use "
        "install/create_conda_envs.sh to create one).\n\n"
        "Original error message: %s" % e)


@click.command()
@click.option('--sample_list', required=True,
              type=click.Path(exists=True, file_okay=True, dir_okay=False,
                              readable=True),
              help="File containing sample directories to merge. Each line "
                   "must be an absolute or relative path to a sample "
                   "directory.")
@click.option('--denoiser', required=True,
              type=click.Choice(['deblur', 'dada2']),
              help="Type of denoiser results to merge.")
@click.option('--trim_len', required=True, type=int,
              help="Trim length of denoiser results to merge.")
@click.option('--KN', 'KN', required=True, type=str,
              help="Type of kmer results to merge (e.g. 6:5000).")
@click.option('--output_dir', required=True,
              type=click.Path(exists=False, file_okay=False, dir_okay=True,
                              writable=True),
              help="Output directory in which to write merged sample data.")
def merge_samples(sample_list, denoiser, trim_len, KN, output_dir):
    """Merge single-sample results for further comparative analyses.

    Merges single-sample ASV feature tables and representative sequences into a
    single ASV feature table and set of representative sequences. The merged
    feature table and representative sequences can be used in downstream
    analyses requiring groups of samples (e.g. comparative analyses such as
    beta diversity, building a phylogenetic tree from a single set of
    representative sequences, etc).

    Additionally, single-sample kmer results and kmer entropies will be merged.
    The merged kmer results will be stored in .npz/_meta format (i.e. the
    format used by the kmer pipeline).

    Single-sample results may only be merged if they are from the same target
    gene, target subfragment, and denoiser/kmer parameters. This ensures that
    the ASVs and kmers being combined across samples are compatible for
    meta-analysis.

    Example output:

    \b
    output_dir/
        asvs/
            rep_seqs.qza
            feature_table.qza
        kmers/
            6-mers_5000.npz
            6-mers_5000_meta
            entropy/
                6-mers_5000_entropy.tsv

    Example usage:

    \b
    python merge_samples.py
        --sample_list samples_to_merge.txt
        --denoiser deblur
        --trim_len 150
        --KN 6:5000
        --output_dir merged_samples

    """
    try:
        kmer_len, subsample_size = KN.split(':')
        kmer_len = int(kmer_len)
        subsample_size = int(subsample_size)
    except ValueError:
        raise click.ClickException(
            "%r is an invalid value for --KN. Please specify a value "
            "following the form <kmer>:<subsample-size> (e.g. 6:5000)." % KN)

    if os.path.exists(output_dir):
        raise click.ClickException(
            "%r directory already exists. Please either remove the directory "
            "and rerun this script or specify a different --output_dir." %
            output_dir)

    os.makedirs(output_dir)

    log = os.path.join(output_dir, 'merge_samples.log')

    # Save pipeline start timestamp for writing to pipeline_completed.log file
    # at end of script.
    pipeline_started = datetime.datetime.now().strftime('%Y-%m-%d:%H:%M:%S')
    log_message("Pipeline execution started.", log)

    # Export a copy of input parameters.
    log_message("Writing file of input parameters started.", log)

    with open(os.path.join(output_dir, "merge_samples_params.log"), 'w') as f:
        f.write("SAMPLE_LIST: %s\n" % sample_list)
        f.write("DENOISER: %s\n" % denoiser)
        f.write("TRIM_LEN: %s\n" % trim_len)
        f.write("KN: %s\n" % KN)
        f.write("OUTPUT_DIR: %s\n" % output_dir)

    log_message("Writing file of input parameters completed.", log)

    log_message("Validating sample directories started.", log)

    sample_data = parse_sample_list(sample_list, denoiser, trim_len, kmer_len,
                                    subsample_size)

    log_message("Validating sample directories completed.", log)

    asv_merge_dir = os.path.join(output_dir, 'asvs')

    log_message(
        "Merging sample ASV data into %r output directory started."
        % asv_merge_dir, log)

    os.makedirs(asv_merge_dir)

    merge_sample_asv_data(sample_data, asv_merge_dir, log)

    log_message(
        "Merging sample ASV data into %r output directory completed."
        % asv_merge_dir, log)

    kmer_merge_dir = os.path.join(output_dir, 'kmers')

    log_message(
        "Merging sample kmer data and entropies into %r output directory "
        "started." % kmer_merge_dir, log)

    os.makedirs(kmer_merge_dir)

    merge_sample_kmer_data(sample_data, kmer_len, subsample_size,
                           kmer_merge_dir)

    log_message(
        "Merging sample kmer data and entropies into %r output directory "
        "completed." % kmer_merge_dir, log)

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
            f.write("merge_samples.py log file\n\n")

    msg = '%s %s' % (datetime.datetime.now().strftime('%Y-%m-%d:%H:%M:%S'),
                     msg)
    click.echo(msg)
    with open(filepath, 'a') as f:
        f.write('%s\n' % msg)


def parse_sample_list(sample_list, denoiser, trim_len, kmer_len,
                      subsample_size):
    """Parse and validate sample directories in sample list."""
    sample_data = []
    sample_ids = set()
    target_subfragments = set()
    bloom_seqs_md5s = set()

    with open(sample_list, 'r') as f:
        for line in f:
            sample_dir = line.strip()

            # Skip blank lines.
            if not sample_dir:
                continue

            sample_id, target_gene, target_subfragment = \
                get_sample_dir_info(sample_dir)

            if target_gene != '16S':
                raise click.ClickException(
                    "Sample %r is %r target gene data. The only supported "
                    "target gene at this time is '16S'." %
                    (sample_id, target_gene))

            if sample_id in sample_ids:
                raise click.ClickException(
                    "Sample list includes a duplicate sample ID %r. Please "
                    "ensure the sample list contains unique sample IDs." %
                    sample_id)

            sample_ids.add(sample_id)
            target_subfragments.add(target_subfragment)

            denoiser_dir = os.path.join(
                sample_dir, 'asvs',
                '%s_single_end_trim_len_%d' % (denoiser, trim_len))
            if not os.path.exists(denoiser_dir):
                raise click.ClickException(
                    "%r directory does not exist. Please run generate_asvs.sh "
                    "with '--denoiser %s' and '--trim_len %d' prior to "
                    "running this script, or specify different "
                    "--denoiser/--trim_len options." %
                    (denoiser_dir, denoiser, trim_len))

            rep_seqs_path = os.path.join(denoiser_dir, 'denoised',
                                         'rep_seqs.qza')
            if not os.path.exists(rep_seqs_path):
                raise click.ClickException(
                    "%r file does not exist. Please run generate_asvs.sh "
                    "prior to running this script." % rep_seqs_path)

            asv_table_path = os.path.join(denoiser_dir, 'denoised',
                                          'feature_table.qza')
            if not os.path.exists(asv_table_path):
                raise click.ClickException(
                    "%r file does not exist. Please run generate_asvs.sh "
                    "prior to running this script." % asv_table_path)

            bloom_seqs_md5 = get_bloom_seqs_md5(denoiser_dir)
            bloom_seqs_md5s.add(bloom_seqs_md5)
            if len(bloom_seqs_md5s) > 1:
                raise click.ClickException(
                    "Sample list contains samples with different bloom "
                    "sequence filters applied to the denoised results being "
                    "merged. Samples may only be merged if the same bloom "
                    "sequence filter was applied to all samples, or if no "
                    "bloom sequence filter was applied to all samples.")

            kmer_dir = os.path.join(denoiser_dir, 'kmers')
            kmer_path = os.path.join(
                kmer_dir,
                '%s_%d-mers_%d.npz' % (sample_id, kmer_len, subsample_size))
            if not os.path.exists(kmer_path):
                raise click.ClickException(
                    "%r file does not exist. Please run generate_kmers.sh "
                    "with '--KN %d:%d' prior to running this script, or "
                    "specify a different --KN option." %
                    (kmer_path, kmer_len, subsample_size))

            entropy_path = os.path.join(
                kmer_dir, 'entropy',
                '%s_%d-mers_%d_entropy.tsv' %
                (sample_id, kmer_len, subsample_size))

            if not os.path.exists(entropy_path):
                raise click.ClickException(
                    "%r file does not exist. Please run "
                    "compute_kmer_entropy.py prior to running this script." %
                    entropy_path)

            sample_data.append({
                'id': sample_id,
                'rep_seqs': rep_seqs_path,
                'asv_table': asv_table_path,
                'kmers': kmer_path,
                'entropy': entropy_path
            })

    if len(target_subfragments) > 1:
        raise click.ClickException(
            "Sample list contains samples from different target subfragments. "
            "Samples may only be merged from the same target subfragment. "
            "Found the following target subfragments in the sample list: %r" %
            target_subfragments)

    if len(sample_data) == 0:
        raise click.ClickException(
            "Sample list is empty. Please specify at least one sample "
            "directory in the file.")

    return sample_data


def get_sample_dir_info(sample_dir):
    """Extract info from sample directory necessary for merging samples."""
    if not os.path.exists(sample_dir):
        raise click.ClickException(
            "%r sample directory does not exist or is not accessible. "
            "Please check that the sample directory path is accessible on "
            "this system." % sample_dir)

    sample_id = get_sample_id(sample_dir)

    qc_dir = os.path.join(sample_dir, 'quality_control')
    if not os.path.exists(qc_dir):
        raise click.ClickException(
            "%r directory does not exist. Please run quality_control.sh prior "
            "to running this script." % qc_dir)

    target_gene, target_subfragment = get_target_info(qc_dir)

    return sample_id, target_gene, target_subfragment


def get_sample_id(sample_dir):
    """Return sample ID based on sample directory name."""
    # normpath removes trailing slashes from directory path. This is
    # necessary prior to using basename.
    return os.path.basename(os.path.normpath(sample_dir))


def get_target_info(qc_dir):
    qc_params_path = os.path.join(qc_dir, 'qc_params.log')

    if not os.path.exists(qc_params_path):
        raise click.ClickException(
            "%r file does not exist. Please run quality_control.sh prior to "
            "running this script." % qc_params_path)

    target_gene = None
    target_subfragment = None
    with open(qc_params_path, 'r') as f:
        for line in f:
            line = line.rstrip('\n')

            if line.startswith('TARGET_GENE:'):
                _, target_gene = line.split()
            elif line.startswith('TARGET_SUBFRAGMENT:'):
                _, target_subfragment = line.split()

    if target_gene is None:
        raise click.ClickException(
            "Could not locate TARGET_GENE parameter in %r" % qc_params_path)
    if target_subfragment is None:
        raise click.ClickException(
            "Could not locate TARGET_SUBFRAGMENT parameter in %r" %
            qc_params_path)

    return target_gene, target_subfragment


def get_bloom_seqs_md5(denoiser_dir):
    params_path = os.path.join(denoiser_dir, 'generate_asvs_params.log')

    if not os.path.exists(params_path):
        raise click.ClickException(
            "%r file does not exist. Please run generate_asvs.sh prior to "
            "running this script." % params_path)

    bloom_seqs_md5 = None
    with open(params_path, 'r') as f:
        for line in f:
            line = line.rstrip('\n')

            if line.startswith('BLOOM_SEQS_MD5:'):
                tokens = line.split()
                if len(tokens) == 1:
                    bloom_seqs_md5 = ''
                else:
                    _, bloom_seqs_md5 = tokens

    if bloom_seqs_md5 is None:
        raise click.ClickException(
            "Could not locate BLOOM_SEQS_MD5 parameter in %r" % params_path)

    return bloom_seqs_md5


def merge_sample_asv_data(sample_data, output_dir, log):
    # We're not using `qiime feature-table merge` and `qiime feature-table
    # merge-seqs` because we ran into serious performance issues when merging a
    # large number of samples (American Gut 21k samples required ~2TB RAM and
    # took over 3 days to run).
    #
    # Merge strategy:
    #
    # 1. Read each sample's representative sequences file and populate a master
    #    dictionary mapping sequence ID to DNA sequence. This dictionary stores
    #    the complete set of ASVs found across the samples.
    #
    # 2. Write out the complete set of ASV representative sequences.
    #
    # 3. Read each sample's ASV feature table and populate a list of scipy
    #    sparse matrices. Pad each sample's matrix with zeros for any ASVs that
    #    are found in the complete set but not in the sample. Sort the sample's
    #    padded data to match the ordering of ASVs in the complete set.
    #
    # 4. Stack all of the per-sample sparse matrices into a single sparse
    #    matrix. Create a biom table from the single sparse matrix for writing
    #    output.
    #
    # Note: there are likely some additional performance optimizations that can
    # be made, such as using multiprocessing to merge smaller groups of samples
    # in parallel. Other parts of the code below may benefit from more in-depth
    # profiling to identify bottlenecks (only high-level profiling was done to
    # fix the most obvious bottlenecks).
    merged_rep_seqs = {}
    for sample in sample_data:
        log_message(
            "Processing representative sequences for sample %r..."
            % sample['id'], log)

        rep_seqs_artifact = qiime2.Artifact.load(sample['rep_seqs'])

        # Manually parse the rep seqs fasta file because it's much faster than
        # using qiime2 to read it into a pandas Series.
        #
        # Note: calling str() on DNAFASTAFormat object returns absolute path to
        # the file.
        rep_seqs_fasta = str(rep_seqs_artifact.view(
            q2_types.feature_data.DNAFASTAFormat))

        for seq_id, seq in iter_fasta_records(rep_seqs_fasta):
            if seq_id not in merged_rep_seqs:
                merged_rep_seqs[seq_id] = seq

    # Sort the complete set of representative sequences, as we'll match this
    # order when padding and merging each sample's ASV table below. Technically
    # we don't need to write the representative sequences in any particular
    # order (i.e. they don't need to match the order of the ASVs in the merged
    # table), but we might as well write the sequences in this order for
    # tidiness.
    merged_asv_ids = sorted(merged_rep_seqs)

    # Manually write the complete set of representative sequences in fasta
    # format, then import into qiime2 artifact (this is much faster than using
    # qiime2 to write a pandas Series directly to an artifact).
    merged_rep_seqs_fasta = os.path.join(output_dir, 'rep_seqs.fasta')
    with open(merged_rep_seqs_fasta, 'w') as f:
        for seq_id in merged_asv_ids:
            f.write('>%s\n' % seq_id)
            f.write('%s\n' % merged_rep_seqs[seq_id])

    merged_rep_seqs_artifact = qiime2.Artifact.import_data(
        q2_types.feature_data.FeatureData[q2_types.feature_data.Sequence],
        merged_rep_seqs_fasta)

    merged_rep_seqs_artifact.save(os.path.join(output_dir, 'rep_seqs.qza'))

    # We no longer need the fasta file because it has a copy saved in
    # rep_seqs.qza.
    os.remove(merged_rep_seqs_fasta)

    # Store all ASV IDs as a set for efficiency because we'll be doing set
    # operations below in the loop.
    merged_asv_ids_set = set(merged_rep_seqs)

    # Build up a list of sample IDs and per-sample sparse matrices. Each sparse
    # matrix will be padded with zeros for any ASVs found in the complete set
    # but not in the sample. Each matrix will have its padded ASVs sorted to
    # match the complete ASV ordering (this makes the matrices stackable by
    # giving them a consistent shape and ordering).
    sample_ids = []
    sparse_matrices = []
    for sample in sample_data:
        sample_id = sample['id']
        log_message("Processing ASV table for sample %r..." % sample_id, log)
        sample_ids.append(sample_id)

        biom_table = qiime2.Artifact.load(sample['asv_table']).view(biom.Table)
        assert list(biom_table.ids(axis='sample')) == [sample_id]

        asv_ids = list(biom_table.ids(axis='observation'))
        asv_data = biom_table.data(sample_id, axis='sample', dense=True)

        # Pad the righthand side of the array with zeros for each missing ASV.
        missing_ids = merged_asv_ids_set - set(asv_ids)
        padded_asv_data = np.pad(asv_data, (0, len(missing_ids)),
                                 mode='constant', constant_values=0.0)
        asv_ids.extend(missing_ids)
        assert len(asv_ids) == len(merged_asv_ids)
        sorted_asv_data = padded_asv_data[np.argsort(asv_ids)]

        sparse_matrices.append(scipy.sparse.csr_matrix([sorted_asv_data]))

    log_message(
        "Merging ASV tables, converting to biom.Table object, importing into "
        "qiime2 artifact...", log)

    # Merge the sparse matrices by stacking them. Convert to biom Table object
    # for writing output.
    merged_asv_data = scipy.sparse.vstack(
        sparse_matrices, format='csr').transpose(copy=False)

    merged_biom_table = biom.Table(merged_asv_data,
                                   observation_ids=merged_asv_ids,
                                   sample_ids=sample_ids)
    merged_table_artifact = qiime2.Artifact.import_data(
        q2_types.feature_table.FeatureTable[q2_types.feature_table.Frequency],
        merged_biom_table)

    log_message("Saving qiime2 merged feature table artifact...", log)

    merged_table_artifact.save(os.path.join(output_dir, 'feature_table.qza'))


def iter_fasta_records(path):
    """Simple fasta parser yielding (sequence ID, sequence) tuples.

    Only very basic validation is performed to keep the parser fast. When
    parsing the sequence header, only the sequence ID is retained (sequence
    description field is excluded).

    """
    with open(path, 'r') as f:
        for line1, line2 in itertools.zip_longest(*[f]*2):
            line1 = line1.rstrip()

            if not line1.startswith('>'):
                raise click.ClickException(
                    "%r file does not appear to be a valid fasta file. "
                    "Expected to find sequence header line starting with '>' "
                    "character, but found: %r" % (path, line1))

            seq_id = line1[1:].split(maxsplit=1)[0]

            if line2 is None:
                raise click.ClickException(
                    "%r file does not appear to be a valid fasta file. Found "
                    "an incomplete record (sequence header line without a "
                    "sequence)." % path)

            seq = line2.rstrip()
            yield seq_id, seq


def merge_sample_kmer_data(sample_data, kmer_len, subsample_size, output_dir):
    sample_ids = []
    sparse_matrices = []
    kmer_entropies = []

    for sample in sample_data:
        sample_id = sample['id']
        sparse_matrix = load_kmer_matrix(sample['kmers'], sample_id, kmer_len)
        kmer_entropy = load_kmer_entropy_results(sample['entropy'], sample_id)

        sample_ids.append(sample_id)
        sparse_matrices.append(sparse_matrix)
        kmer_entropies.append(kmer_entropy)

    # Merge single-sample sparse matrices (each 1xN matrix becomes a row).
    merged_matrix = scipy.sparse.vstack(sparse_matrices, format='csr')

    write_kmer_matrix(merged_matrix, sample_ids, kmer_len, subsample_size,
                      output_dir)

    output_entropy_dir = os.path.join(output_dir, 'entropy')
    os.makedirs(output_entropy_dir)

    output_entropy_path = os.path.join(
        output_entropy_dir,
        '%d-mers_%d_entropy.tsv' % (kmer_len, subsample_size))

    save_kmer_entropy_results(sample_ids, kmer_entropies, output_entropy_path)


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
            "pipeline results are expected to be single-sample when merging."
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


def write_kmer_matrix(sparse_matrix, sample_ids, kmer_len, subsample_size,
                      output_dir):
    """Write kmer sparse matrix to disk in .npz/_meta format."""
    path_prefix = os.path.join(output_dir,
                               '%d-mers_%d' % (kmer_len, subsample_size))

    npz_path = path_prefix + '.npz'
    meta_path = path_prefix + '_meta'

    np.savez(npz_path, data=sparse_matrix.data, indices=sparse_matrix.indices,
             indptr=sparse_matrix.indptr, shape=sparse_matrix.shape)

    with open(meta_path, 'w') as f:
        for sample_id in sample_ids:
            f.write('%s\n' % sample_id)


def load_kmer_entropy_results(entropy_path, sample_id):
    data = []

    # Newline settings based on recommendation from csv docs:
    #     https://docs.python.org/3/library/csv.html#id3
    with open(entropy_path, 'r', newline='', encoding='utf-8') as fh:
        tsv_reader = csv.reader(fh, dialect='excel-tab', strict=True)

        header = next(tsv_reader)
        if header != ['sample-id', 'shannon-entropy']:
            raise click.ClickException(
                "%r file does not have the expected header of a kmer entropy "
                "file. Expected column headers are 'sample-id' and "
                "'shannon-entropy'." % entropy_path)

        for row in tsv_reader:
            id_, entropy = row
            data.append((id_, float(entropy)))

    if len(data) != 1:
        raise click.ClickException(
            "%r file has entropy results for %d samples. Expected to find "
            "results for a single sample only." % (entropy_path, len(data)))

    if data[0][0] != sample_id:
        raise click.ClickException(
            "%r file has entropy results for sample %r. Expected to find "
            "results for sample %r." % (entropy_path, data[0][0], sample_id))

    return data[0][1]


# TODO `save_kmer_entropy_results` is similar to the function in
# compute_kmer_entropy.py. Consider centralizing this code in `util` directory
# in the future.

def save_kmer_entropy_results(sample_ids, kmer_entropies, output_path):
    header = ['sample-id', 'shannon-entropy']

    # Newline settings based on recommendation from csv docs:
    #     https://docs.python.org/3/library/csv.html#id3
    with open(output_path, 'w', newline='', encoding='utf-8') as fh:
        tsv_writer = csv.writer(fh, dialect='excel-tab', strict=True)
        tsv_writer.writerow(header)

        for sample_id, kmer_entropy in zip(sample_ids, kmer_entropies):
            # Format floats with *up to* 15 digits of total precision.
            # Scientific notation may be used if appropriate. Integers won't
            # have trailing zeros included (e.g. 0.0 will be written as "0";
            # 42.0 will be written as "42", etc).
            tsv_writer.writerow((sample_id, '{0:.15g}'.format(kmer_entropy)))


if __name__ == '__main__':
    merge_samples()
