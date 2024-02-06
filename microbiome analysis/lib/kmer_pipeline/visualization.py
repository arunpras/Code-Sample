import scipy
import numpy as np
from scipy import sparse
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from classifier.DNN import DNNMutliclass16S
from utility.file_utility import FileUtility
import matplotlib
import matplotlib.pyplot as plt
import getopt
import codecs
import sys
import ast
import types
import yaml
from collections import OrderedDict


def plot_scatter(X, Y, x_label, y_label, to_annotate, title):
    global color_schemes
    color_list = color_schemes

    target = list(set(Y))

    target.sort()
    color_idx = [target.index(x) for x in Y]

    for i, value in enumerate(X):
        plt.scatter(
            value[0],
            value[1],
            label=target[color_idx[i]],
            c=color_list[color_idx[i]],
            cmap='viridis',
            alpha=0.4,
            edgecolors=None)

    for dataset in to_annotate:
        start = to_annotate[dataset]['start']
        end = to_annotate[dataset]['stop']
        range_to_annotate = list(range(start, end + 1))
        for i in range_to_annotate:
            plt.annotate(dataset, (X[i - 1, 0], X[i - 1, 1]), size=10)

    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = OrderedDict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys(), loc=0)

    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.xticks([])
    plt.yticks([])


def plot_pca_tsne_nn(X_pca, X_tsne, X_tsne_NN, Y, filename, to_annotate):

    myplot = plt.figure(figsize=(12, 12))
    #PCA plot
    #plot_scatter(ax, X_pca, Y, 'PCA_1', 'PCA_0', '(i) PCA over 6-mer representations')

    #regular TSNE plot
    #plot_scatter(ax, X_tsne, Y, 't-SNE_1', 't-SNE_0', '(ii) t-SNE over 6-mer representations')

    plot_scatter(
        X_tsne_NN, Y, 't-SNE_1', 't-SNE_0', to_annotate,
        '(i) t-SNE / activation function of the last layer of neural network')

    plt.savefig(filename + '.pdf')


def get_pca_tsne(X, X_NN, random_seed=0):
    X_pca = PCA(n_components=30).fit_transform(X.toarray())
    X_tsne = TSNE(
        n_components=2,
        perplexity=40,
        verbose=2,
        learning_rate=10,
        random_state=random_seed).fit_transform(X.toarray())

    X_tsne_NN = TSNE(
        n_components=2,
        perplexity=40,
        verbose=2,
        learning_rate=10,
        random_state=random_seed).fit_transform(X_NN)

    return X_pca, X_tsne, X_tsne_NN


global color_schemes
color_schemes = [
    '#ff0505', '#f2a041', '#cdff05', '#04d9cb', '#45a8ff', '#8503a6',
    '#590202', '#734d02', '#4ab304', '#025359', '#0454cc'
]

yaml_file = sys.argv[1]

with open(yaml_file, 'r') as stream:
    try:
        file = yaml.load(stream)

        kmer_file = file['kmer']
        Y_file = file['metadata']
        X_NN_file = file['trained_NN']
        to_annotate = file['annotation']
        output_filename = file['output']
    except yaml.YAMLError as exc:
        print(exc)

Y_filename = Y_file
loader = np.load(kmer_file)

X = sparse.csr_matrix(
    (loader['data'], loader['indices'], loader['indptr']),
    shape=loader['shape'])

X_NN = DNNMutliclass16S.make_activation_function(X_NN_file, X.toarray())
X_pca_bs, X_tsne_bs, X_tsne_NN_bs = get_pca_tsne(X, X_NN)
Y = [
    line.strip() for line in codecs.open(Y_filename, 'r', 'utf-8').readlines()
]

plot_pca_tsne_nn(X_pca_bs, X_tsne_bs, X_tsne_NN_bs, Y, output_filename,
                 to_annotate)
