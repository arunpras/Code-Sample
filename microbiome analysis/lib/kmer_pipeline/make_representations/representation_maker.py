__author__ = "Arun Manoharan"
__license__ = "Apache 2"
__version__ = "2.0.0"
__maintainer__ = "Arun Manoharan"
__email__ = "arun@primediscoveries.com"


import gzip
import sys
sys.path.append('../')
from sklearn.preprocessing import normalize
from sklearn.feature_extraction.text import TfidfVectorizer
import itertools
import numpy as np
from multiprocessing import Pool
import tqdm
import random
from scipy import sparse
from utility.file_utility import FileUtility
from Bio import SeqIO


class Metagenomic16SRepresentation:
    '''
        Make k-mer from directory of fasta files
    '''

    def __init__(self, fasta_files, indexing, sampling_number=3000, num_p=20,
                 include_low_depth_files=False):
        '''
        :param fasta_files: list of fasta files
        :param indexing: the index
        :param sampling_number:
        :param num_p:
        '''
        self.fasta_files=fasta_files
        self.num_p=num_p
        self.sampling_number=sampling_number
        self.indexing=indexing
        self.include_low_depth_files = include_low_depth_files

    def generate_kmers_all(self, k, save=False):
        '''
        :param k:
        :param save:
        :return:
        '''
        self.k=k
        self.vocab = [''.join(xs) for xs in itertools.product('atcg', repeat=k)]
        self.vocab.sort()
        self.vectorizer = TfidfVectorizer(use_idf=False, vocabulary=self.vocab, analyzer='char', ngram_range=(k, k),
                                          norm=None, stop_words=[], lowercase=True, binary=False)

        data = np.zeros((len(self.fasta_files), len(self.vocab))).astype(np.float64)

        # multi processing extraction of k-mer distributions
        processed_indices = []
        s_steps=[]
        pool = Pool(processes=self.num_p)
        if save:
            log_file = open(save + '_log', 'w')

        for ky, (v,s) in tqdm.tqdm(pool.imap_unordered(self.get_kmer_distribution, self.fasta_files, chunksize=1),
                                   total=len(self.fasta_files)):
            if v is None:
                msg = (
                    "Ignoring file %r while generating %s-mers with sampling "
                    "size %d. The file must have at least as many sequences "
                    "as the sampling size to be included in kmer generation, "
                    "but the file only has %d sequences.\n" %
                    (ky, k, self.sampling_number, s))
                print(msg)
                if save:
                    log_file.write(msg + '\n')
            else:
                index = self.indexing[ky]
                data[index, :] = v
                processed_indices.append(index)
                s_steps.append(s)

        if len(processed_indices) == 0:
            msg = (
                "All files were ignored while generating %d-mers with "
                "sampling size %d. None of the files had at least as many "
                "sequences as the sampling size.\n" %
                (k, self.sampling_number))
            print(msg)
            if save:
                log_file.write(msg + '\n')

            fasta_files = None
            data = None
        else:
            # Filter list of fasta files and data rows to include only files
            # that were processed.
            processed_indices.sort()
            fasta_files = [self.fasta_files[i] for i in processed_indices]
            data = data[processed_indices]

            # normalize the frequencies
            data = normalize(data, axis=1, norm='l1')
            data = sparse.csr_matrix(data)

            if save:
                FileUtility.save_sparse_csr(save, data)
                FileUtility.save_list(save+'_meta',fasta_files)
                log_file.write('mean_size: %s\n' % np.mean(s_steps))
                log_file.write('std_size: %s\n' % np.std(s_steps))

        if save:
            log_file.close()

        return fasta_files, data

    def open_seq_file(self, file_name):
        if file_name.endswith('gz'):
            f = gzip.open(file_name, 'rt')

        else:
            f = open(file_name, 'rt')

        return f

    def fragment_sequence(self, seq, kmer_dict):
        for i in range(len(seq) - self.k + 1):

            fragment = seq[i:i+self.k]
            if 'n' in fragment:
                continue

            kmer_dict[fragment] += 1

        return kmer_dict

    def get_kmer_distribution(self, file_name):
        '''

        :param file_name:
        :return:
        '''

        vector = None
        kmer_dict = dict(zip(self.vocab, [0] * len(self.vocab)))

        corpus = []
        with self.open_seq_file(file_name) as f:

            if self.sampling_number == -1:
                if file_name.endswith('fastq') or file_name.endswith('fastq.gz'): 
                    for n, line in enumerate(f):
                        if n % 4 == 1:
                            kmer_dict = self.fragment_sequence(line.strip().lower(), kmer_dict)

                    tot_size = int((n + 1) / 4)

                else:
                    for n, line in enumerate(f):
                        if n % 2 == 1:
                            kmer_dict = self.fragment_sequence(line.strip().lower(), kmer_dict)

                    tot_size = int((n + 1) / 2)

                vector = np.array([kmer_dict[v] for v in self.vocab])

            else:
                with self.open_seq_file(file_name) as fcount:
                    for n, line in enumerate(fcount):
                        pass

                if file_name.endswith('fastq') or file_name.endswith('fastq.gz'):
                    seq_indices = range(1, n + 1, 4)

                else:
                    seq_indices = range(1, n + 1, 2)

                tot_size = len(seq_indices)

                if self.include_low_depth_files or (tot_size >= self.sampling_number):
                    sampled_seq_indices = random.sample(seq_indices, min(self.sampling_number, tot_size))
                    sampled_seq_indices.sort(reverse=True)

                    for n, line in enumerate(f):
                        if n == sampled_seq_indices[-1]:
                            kmer_dict = self.fragment_sequence(line.strip().lower(), kmer_dict)
                            sampled_seq_indices.pop()

                            if len(sampled_seq_indices) == 0:
                                break

                    vector = np.array([kmer_dict[v] for v in self.vocab])

        return file_name, (vector, tot_size)

class FastaRepresentations(object):
    '''
        Make k-mer from single fasta file
        where the headers contain info about the label
    '''
    def __init__(self, fasta_address, label_modifying_func=str):
        '''
        :param fasta_address:
        :param label_modifying_func: extract label from the header
        '''
        self.labels=[]
        self.corpus=[]
        for cur_record in SeqIO.parse(fasta_address, 'fasta'):
            self.corpus.append(str(cur_record.seq).lower())
            self.labels.append(str(cur_record.id).lower())
        self.labels=[label_modifying_func(l) for l in self.labels]

    def get_samples(self, envs, N):
        '''
        :param envs: list of labels
        :param N: sample size
        :return: extract stratified with size N corpus and label list
        '''
        labels=[]
        corpus=[]
        for env in envs:
            selected=[idx for idx,v in enumerate(self.labels) if env==v]
            if N==-1:
                N=len(selected)
            idxs=random.sample(selected, N)
            corpus=corpus+[self.corpus[idx] for idx in idxs]
            labels=labels+[self.labels[idx] for idx in idxs]
        return corpus, labels

    def get_vector_rep(self, corpus, k, restricted=True):
        '''
        :param corpus:
        :param k: k-mer size
        :param restricted: restricted to known values
        :return:
        '''
        if restricted:
            vocab = [''.join(xs) for xs in itertools.product('atcg', repeat=k)]
            tf_vec = TfidfVectorizer(use_idf=True, vocabulary=vocab, analyzer='char', ngram_range=(k, k),
                                                  norm='l1', stop_words=[], lowercase=True, binary=False)
        else:
            tf_vec = TfidfVectorizer(use_idf=True, analyzer='char', ngram_range=(k, k),
                                                  norm='l1', stop_words=[], lowercase=True, binary=False)
        return tf_vec.fit_transform(corpus)

if __name__=='__main__':
    FR=FastaRepresentations('sample.fasta')
    MR=Metagenomic16SRepresentation('16ssamples/')

