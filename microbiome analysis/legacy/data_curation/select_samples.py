"""Written by Arun Manoharan
arun@primediscoveries.com

Selects the samples you want
"""
import csv
import numpy as np
from argparse import ArgumentParser
from scipy import sparse

def load_kmers(filename):
    loader = np.load(filename)
    csr_matrix = sparse.csr_matrix(
        (loader['data'], loader['indices'], loader['indptr']),
        shape=loader['shape'])
    kmers = csr_matrix.toarray()
    return kmers

def save_kmers(filename, array):
    np.savez(
        filename,
        data=array.data,
        indices=array.indices,
        indptr=array.indptr,
        shape=array.shape)

def get_metadata_dict(sample_metadata, column_of_interest, desired_labels):
    """
    Iterates through the sample metadata, and returns a dictionary of samplename to label of
    the interested column
    :param sample_metadata: the sample metadata file
    :param column_of_interest: the column you are interested in (eg. body_habitat)
    :return: a dictionary of the sample_name to the label
    """
    print('Selecting following labels from: ' + column_of_interest)
    for label in desired_labels:
        print(label)
    Y = np.genfromtxt(sample_metadata, delimiter='\t', dtype='str',autostrip=True)
    #kinda convoluted line because np.where returns a tuple of arrays
    col = np.where(Y[0]==column_of_interest)[0][0]
    samplename_interestedcol = Y[1:, [0,col]]
    sample_to_meta = {}
    for row in samplename_interestedcol:
        sample_name = row[0]
        label = row[1]
        if desired_labels == 'all':
            sample_to_meta[sample_name] = label
        elif label in desired_labels:
            sample_to_meta[sample_name] = label

    return sample_to_meta

def select_and_save_kmers(metadata_dict, kmer_metadata, kmer_filename, kmer_save_filename, meta_save_filename):
    """
    Selects the kmers with the desired labels, then saves them along with the
    :param metadata_dict: dictionary from sample_name to metadata label
    :param kmer_metadata: filename to the kmer metadata
    :param kmer_filename: filename to the kmers
    :param save_filename: filename to save to
    :return:
    """
    kmers = load_kmers(kmer_filename)
    kmer_metadata = np.genfromtxt(kmer_metadata, delimiter='\t', dtype='str',autostrip=True)
    kmer_indices = []
    metadata = []
    for index, sample_name in enumerate(kmer_metadata):
        if sample_name in metadata_dict:
            kmer_indices.append(index)
            metadata.append([sample_name, metadata_dict[sample_name]])

    with open(meta_save_filename, "w+") as mycsv:
        csvWriter = csv.writer(mycsv, delimiter =',')
        csvWriter.writerows(metadata)

    selected_kmers = kmers[kmer_indices]
    selected_kmers= sparse.csr_matrix(selected_kmers)
    save_kmers(kmer_save_filename, selected_kmers)

    print('Saved '+str(selected_kmers.shape[0]) +' kmers')


parser = ArgumentParser()
parser.add_argument("--kmer_file", dest= "kmers", help='filepath to kmers')
parser.add_argument('--kmer_meta', dest ='kmer_meta', help ='filepath to kmer metadata')
parser.add_argument('--sample_meta', dest='sample_meta', help ='filepath to sample metadata')
parser.add_argument('--meta_col', dest='column_of_interest', help='metadata column you want')
parser.add_argument('--desired_features', dest='desired_features', help='features you want, use all for all')
parser.add_argument('--meta_save_dest', dest='meta_save_dest', help='destination to save metadata to')
parser.add_argument('--kmer_save_dest', dest='kmer_save_dest', help='destination to save kmers to')

args = parser.parse_args()

kmers = args.kmers
kmer_metadata = args.kmer_meta
sample_metadata = args.sample_meta
column_of_interest = args.column_of_interest
desired_features = args.desired_features.split(',')
meta_save_file = args.meta_save_dest
kmer_save_file = args.kmer_save_dest

sample_to_meta = get_metadata_dict(sample_metadata=sample_metadata, column_of_interest=column_of_interest,
                                   desired_labels=desired_features)
select_and_save_kmers(metadata_dict=sample_to_meta, kmer_metadata=kmer_metadata, kmer_filename=kmers, kmer_save_filename = kmer_save_file, meta_save_filename=meta_save_file)

