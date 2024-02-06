#!/usr/bin/env python3

__author__ = "Jai Ram Rideout (jai@primediscoveries.com)"

import gzip
import itertools
import os
import os.path
import re
import sys


_whitespace_regex = re.compile(r'\s')


def concat_paired_end(r1_path, r2_path, output_path):
    """Concatenate paired-end R1 and R2 reads into a single fastq file.

    Concatenates paired-end R1 and R2 .fastq.gz files into a single .fastq.gz
    file. Sequence headers in the concatenated file will have all whitespace
    characters replaced with underscores, and '_R1' or '_R2' will be appended
    to the end of the header to differentiate between R1 and R2 reads.

    This script mainly exists for compatibility with humann2/metaphlan2. The
    tools require paired-end reads to be concatenated into a single file, and
    also require sequence headers to be unique. There are cases (e.g. with NCBI
    SRA data) where duplicate sequence headers exist after concatenating
    paired-end reads, which produces incorrect results with humann2/metaphlan2.

    See the following references/issues for details:

    https://forum.biobakery.org/t/same-id-for-paired-end-reads/65

    https://groups.google.com/d/topic/humann-users/QSoO_7C-Bvw/discussion

    https://bitbucket.org/biobakery/humann2/wiki/Home#markdown-header-humann2-
    and-paired-end-sequencing-data

    Note: This is an internal script for the shotgun pipeline and is not
    intended to be executed directly.

    """
    if os.path.exists(output_path):
        raise Exception("Output path %r already exists." % output_path)

    headers = set()
    with gzip.open(output_path, 'wt') as output_file:
        for input_path, read_indicator in (r1_path, '_R1'), (r2_path, '_R2'):
            for header, seq, qual in iter_fastq_records(input_path):
                new_header = (
                    _whitespace_regex.sub('_', header) + read_indicator)

                if new_header in headers:
                    # This error shouldn't happen in practice.
                    raise Exception(
                        "Found duplicate sequence header %r. The duplicate "
                        "either exists in the input fastq file, or a "
                        "duplicate exists after replacing whitespace with "
                        "underscores and appending %r to the header." %
                        (new_header, read_indicator))
                else:
                    headers.add(new_header)

                output_file.write('%s\n' % new_header)
                output_file.write('%s\n' % seq)
                output_file.write('+\n')
                output_file.write('%s\n' % qual)


def iter_fastq_records(path):
    with gzip.open(path, 'rt') as fh:
        for lines in itertools.zip_longest(*[fh] * 4):
            if any(line is None for line in lines):
                raise Exception(
                    "Number of lines in FASTQ file must be multiple of four "
                    "(i.e. each record must be exactly four lines long).")

            seq_header, seq, qual_header, qual = lines

            if not seq_header.startswith('@'):
                raise Exception(
                    "Invalid FASTQ sequence header: %r" % seq_header)
            if not qual_header.startswith('+'):
                raise Exception(
                    "Invalid FASTQ quality header: %r" % qual_header)
            if qual == '\n':
                raise Exception("FASTQ record is missing quality scores.")

            yield seq_header.rstrip('\n'), seq.rstrip('\n'), qual.rstrip('\n')


if __name__ == '__main__':
    # Could use something fancier like argparse or Click, but this works fine
    # considering it's an internal script.
    try:
        r1_path, r2_path, output_path = sys.argv[1:]
    except Exception:
        sys.exit(
            "Could not parse script arguments.\n\n"
            "Usage: concat_paired_end.py <input R1 .fastq.gz> "
            "<input R2 .fastq.gz> <output .fastq.gz>")

    concat_paired_end(r1_path, r2_path, output_path)
