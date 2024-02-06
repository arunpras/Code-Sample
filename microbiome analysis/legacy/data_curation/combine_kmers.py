"""Written by Arun Manoharan
arun@primediscoveries.com

Combines two kmer files

usage:
    python combine_kmers,py
        --kmer_1            #first kmer file
        --kmer_2            #additional kmers to add to kmer_1
        --num_add_x         #how many samples from kmer_1 to add to kmer_2. Useful when you only want a select number from an additional dataset
        --filename          #name of output file
    

example usage:
    cd /home/arun/MicroPheno/test_kmer_combine
    python ../combine_kmers.py --kmer_1='../449_kmer/kmers/449_6-mers_1000.npz' --kmer_2='../1939_kmer/1939_6-mers_1000.npz' --num_add_x='200' --filename='test'

"""
import argparse
import sys
import numpy as np
from scipy import sparse


'''Checks and parses the system arguments

Args:
    args: sys.argv
'''
def checkArgs(args):
    err = ""
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--kmer_1',
        action='store',
        dest='kmer_1',
        help='first kmer file',
        default=False,
        type=str)

    parser.add_argument(
        '--kmer_2',
        action='store',
        dest='kmer_2',
        help='additional kmers to add to kmer_1',
        default=False,
        type=str)

    parser.add_argument(
        '--num_add_x',
        action='store',
        dest='num_add_X',
        help='how many samples from kmer_1 to add to kmer_2. Useful when you only want a select number from an additional dataset',
        default='-1',
        type=int)

    parser.add_argument(
        '--filename',
        action='store',
        dest='filename',
        help='name of output file',
        type=str)
    return parser.parse_args()


'''Loads a sparse csr

Args:
    filename: name of input file
'''
def load_sparse_csr(filename):
    loader = np.load(filename)
    return sparse.csr_matrix(
        (loader['data'], loader['indices'], loader['indptr']),
        shape=loader['shape'])


'''From the npz file, loads it into a csr

Args:
    X_file: main kmer file
    add_X_file: kmer file to add to X_file. Must have same col len as X_file

'''
def combine_kmers(X_file, add_X_file, num_add_X):
    X = load_sparse_csr(X_file).toarray()
    add_X = load_sparse_csr(add_X_file).toarray()
    if num_add_X != -1:
        add_X=add_X[0:num_add_X]
    new_X = sparse.csr_matrix(sparse.vstack([X, add_X]).toarray())
    return new_X


'''Saves an array into filename

Args:
    filename: filename to save into
    array: array to save
'''
def save_sparse_csr(filename, array):
    np.savez(
        filename,
        data=array.data,
        indices=array.indices,
        indptr=array.indptr,
        shape=array.shape)

parsed_args = checkArgs(sys.argv)
new_X = combine_kmers(parsed_args.kmer_1, parsed_args.kmer_2, parsed_args.num_add_X)
save_sparse_csr(parsed_args.filename, new_X)

