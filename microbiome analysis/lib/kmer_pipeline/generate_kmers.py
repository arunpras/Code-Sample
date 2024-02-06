#!/usr/bin/env python3

import datetime
import os
import argparse
import sys
import os.path
import fnmatch

try:
    from utility.file_utility import FileUtility
    from make_representations.representation_maker import Metagenomic16SRepresentation
except ImportError as e:
    sys.exit(
        "Couldn't import required dependencies for kmer pipeline. Please "
        "ensure that a kmer-pipeline conda environment is activated prior to "
        "running this script (use install/create_conda_envs.sh to create one)."
        "\n\nOriginal error message: %s" % e)


def representation_creation_dir(
        inp_dir,
        out_dir,
        dataset_name,
        num_p,
        filetype='fastq.gz',
        sampling_dict={
            3: [20],
            4: [100],
            5: [500],
            6: [100, 1000, 2000, 5000, 10000, -1],
            7: [5000],
            8: [8000]
        },
        include_low_depth_files=False):

    pipeline_started = datetime.datetime.now().strftime('%Y-%m-%d:%H:%M:%S')
    completed_file = os.path.join(out_dir, 'pipeline_completed.log')
    if os.path.exists(completed_file):
        os.remove(completed_file)

    fasta_files, mapping = FileUtility.read_fasta_directory(inp_dir, filetype)

    for k in sampling_dict.keys():
        for N in sampling_dict[k]:
            print(k, '-mers with sampling size ', N)
            RS = Metagenomic16SRepresentation(
                fasta_files, mapping, N, num_p,
                include_low_depth_files=include_low_depth_files)
            # path to save the generated files
            RS.generate_kmers_all(
                k,
                save=os.path.join(out_dir, '_'.join(
                    [dataset_name, str(k) + '-mers',
                     str(N)])))

    pipeline_completed = datetime.datetime.now().strftime('%Y-%m-%d:%H:%M:%S')
    with open(completed_file, 'w') as f:
        f.write('%s Pipeline execution started.\n' % pipeline_started)
        f.write('%s Pipeline execution completed.\n' % pipeline_completed)


def check_args(args):

    err = ""
    # Using the argument parser in case of -h or wrong usage the correct argument usage
    # will be prompted
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--inaddr',
        action='store',
        dest='genrep_input_addr',
        default=False,
        type=str,
        help=
        'genkmer: Generate representations for input fasta file or directory of 16S rRNA samples',
        required='--genkmer' in sys.argv)

    parser.add_argument(
        '--filetype',
        action='store',
        dest='filetype',
        type=str,
        default='fastq',
        help='fasta fsa fastq etc')

    parser.add_argument(
        '--out',
        action='store',
        dest='output_addr',
        type=str,
        default='out',
        help='Out put directory')

    parser.add_argument(
        '--name',
        action='store',
        dest='data_name',
        type=str,
        default=None,
        help='name of the dataset')

    parser.add_argument(
        '--KN',
        action='store',
        dest='K_N',
        default=None,
        type=str,
        help=
        'pair of comma separated Kmer:sub-sample-size ==> 2:100,6:-1 (N=-1 means using all sequences)'
    )
    parser.add_argument(
        '--include-low-depth-files',
        action='store_true',
        help='Include files that have fewer sequences than the subsampling '
             'size specified in --KN. If this flag is supplied, all sequences '
             'in the file will be used if it has fewer sequences than the '
             'subsampling size. By default, files with fewer sequences than '
             'the subsampling size will not be processed. Be cautious with '
             'using this flag, as it allows kmers to potentially be generated '
             'from an uneven number of sequences.'
    )
    parser.add_argument(
        '--cores',
        action='store',
        dest='cores',
        default=4,
        type=int,
        help='Number of cores to be used')

    parsedArgs = parser.parse_args()

    if (not os.access(parsedArgs.genrep_input_addr, os.F_OK)):
        err = err + "\nError: Permission denied or could not find the directory!"
        return err
    elif os.path.isdir(parsedArgs.genrep_input_addr):
        print('Representation creation requested for directory ' +
              parsedArgs.genrep_input_addr + '\n')
        try:
            os.stat(parsedArgs.output_addr)
        except:
            os.mkdir(parsedArgs.output_addr)

        if len(
                FileUtility.recursive_glob(parsedArgs.genrep_input_addr,
                                           '*' + parsedArgs.filetype)) == 0:
            err = err + "\nThe filetype " + parsedArgs.filetype + " could not find the directory!"
            return err

        if not parsedArgs.data_name:
            parsedArgs.data_name = parsedArgs.genrep_input_addr.split('/')[-1]

        try:
            sampling_dict = dict()
            for x in parsedArgs.K_N.split(','):
                k, n = x.split(':')
                k = int(k)
                n = int(n)
                if k in sampling_dict:
                    sampling_dict[k].append(n)
                else:
                    sampling_dict[k] = [n]
        except:
            err = err + "\nWrong format for KN (k-mer sample sizes)!"
            return err

        representation_creation_dir(
            parsedArgs.genrep_input_addr,
            parsedArgs.output_addr,
            parsedArgs.data_name,
            parsedArgs.cores,
            filetype=parsedArgs.filetype,
            sampling_dict=sampling_dict,
            include_low_depth_files=parsedArgs.include_low_depth_files)
    else:
        print('Representation creation requested for file ' +
              parsedArgs.genrep_input_addr + '\n')


if __name__ == '__main__':
    err = check_args(sys.argv)
    if err:
        print(err)
        exit()
