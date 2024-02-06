__author__ = "Arun Manoharan"
__license__ = "Apache 2"
__version__ = "2.0.0"
__maintainer__ = "Arun Manoharan"
__email__ = "arun@primediscoveries.com"

import codecs
import numpy as np
from scipy import sparse
import fnmatch
import os
import _pickle as pickle
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from biom import example_table
from biom import load_table

class FileUtility(object):
    def __init__(self):
        print ('File utility object created..')

    @staticmethod
    def create_fasta_file(file_address, corpus, label):
        seq_id_pairs=[('.'.join([str(idx+1),label[idx]]),x) for idx, x in enumerate(corpus)]
        seq_recs=[ SeqRecord(Seq(seq,generic_dna),id=id, description='') for id,seq in seq_id_pairs]
        SeqIO.write(seq_recs, file_address, "fasta")


    @staticmethod
    def read_fasta_directory(file_directory, file_extenstion, only_files=[]):
        '''
        :param file_directory:
        :param file_extenstion:
        :param only_files:
        :return: list of fasta files, and a dic to map file to index
        '''
        if len(only_files) > 0:
            fasta_files = [x for x in FileUtility.recursive_glob(file_directory, '*.' + file_extenstion) if
                           x.split('/')[-1] in only_files]
        else:
            fasta_files = [x for x in FileUtility.recursive_glob(file_directory, '*.' + file_extenstion)]

        fasta_files.sort()
        mapping = {v: k for k, v in enumerate(fasta_files)}
        return fasta_files, mapping

    @staticmethod
    def read_OTU_format(biom_file):
        '''
        return OTU content
        :param biom_file:
        :return:
        '''
        table = load_table(biom_file)
        X_otu=table.matrix_data
        OTU_ID_Mapping={x.split('.')[1]:idx for idx,x in enumerate(list(table.ids()))}
        return X_otu, OTU_ID_Mapping

    @staticmethod
    def save_obj(filename, value):
        with open(filename + '.pickle', 'wb') as f:
            pickle.dump(value, f)

    @staticmethod
    def load_obj(filename):
        return pickle.load(open(filename, "rb"))

    @staticmethod
    def save_list(filename, list_names):
        f = codecs.open(filename, 'w', 'utf-8')
        for x in list_names:
            f.write(x + '\n')
        f.close()

    @staticmethod
    def load_list(filename):
        return [line.strip() for line in codecs.open(filename, 'r', 'utf-8').readlines()]

    @staticmethod
    def save_sparse_csr(filename, array):
        np.savez(filename, data=array.data, indices=array.indices,
                 indptr=array.indptr, shape=array.shape)

    @staticmethod
    def load_sparse_csr(filename):
        loader = np.load(filename)
        return sparse.csr_matrix((loader['data'], loader['indices'], loader['indptr']), shape=loader['shape'])

    @staticmethod
    def _float_or_zero(value):
        try:
            return float(value)
        except:
            return 0.0

    @staticmethod
    def recursive_glob(treeroot, pattern):
        '''
        :param treeroot: the path to the directory
        :param pattern:  the pattern of files
        :return:
        '''
        results = []
        for base, dirs, files in os.walk(treeroot):
            good_files = fnmatch.filter(files, pattern)
            results.extend(os.path.join(base, f) for f in good_files)
        return results
