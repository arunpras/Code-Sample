#Written By arun@primediscoveries.com
import os
import argparse
import sys
import os.path
import fnmatch
import keras

from classifier.DNN import DNNMutliclass16S
from utility.file_utility import FileUtility
from classifier.classical_classifiers import RFClassifier, SVM


def DNN_classifier(X_file, Y_file, arch, out_dir, dataset_name, gpu_id, epochs,
                   batch_size):
    # k-mer data
    X = FileUtility.load_sparse_csr(X_file).toarray()
    print(X)
    print(X.shape)
    # labels
    Y = FileUtility.load_list(Y_file)
    DNN = DNNMutliclass16S(X, Y, model_arch=arch)
    DNN.cross_validation(
        out_dir + 'nn_classification_results_' + dataset_name,
        gpu_dev=gpu_id,
        n_fold=2,
        epochs=epochs,
        batch_size=batch_size,
        model_strct='mlp')


def classical_classifier(X_file, Y_file, model, out_dir, dataset_name, cores):
    #
    X = FileUtility.load_sparse_csr(X_file)
    # labels
    Y = FileUtility.load_list(Y_file)

    if model == 'RF':
        #### Random Forest classifier
        MRF = RFClassifier(X, Y)
        # results containing the best parameter, confusion metrix, best estimator, results on fold will be stored in this address
        MRF.tune_and_eval(
            out_dir + '/classification_results_' + dataset_name, n_jobs=cores)
    else:
        #### Support Vector Machine classifier
        MSVM = SVM(X, Y)
        # results containing the best parameter, confusion metrix, best estimator, results on fold will be stored in this address
        MSVM.tune_and_eval(
            out_dir + '/classification_results_' + dataset_name, n_jobs=cores)


def check_args(args):

    err = ""
    # Using the argument parser in case of -h or wrong usage the correct argument usage
    # will be prompted
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--model',
        action='store',
        dest='model',
        type=str,
        default=False,
        choices=[False, 'RF', 'SVM', 'DNN'],
        help='train_predictor: choice of classifier from RF, SVM, DNN')

    parser.add_argument(
        '--x',
        action='store',
        dest='X',
        type=str,
        default=False,
        help=
        'train_predictor: The data in the npy format rows are instances and columns are features'
    )

    parser.add_argument(
        '--y',
        action='store',
        dest='Y',
        type=str,
        default=False,
        help=
        'train_predictor: The labels associated with the rows of classifyX, each line is a associated with a row'
    )

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
        '--epochs',
        action='store',
        dest='epochs',
        type=int,
        default=100,
        help='train_predictor-model/DNN: number of epochs for deep learning')

    parser.add_argument(
        '--arch',
        action='store',
        dest='dnn_arch',
        type=str,
        default='1024,0.2,512',
        help=
        'train_predictor-model/DNN: The comma separated definition of neural network layers connected to eahc other, you do not need to specify the input and output layers, values between 0 and 1 will be considered as dropouts'
    )
    parser.add_argument(
        '--gpu_id',
        action='store',
        dest='gpu_id',
        type=str,
        default='0',
        help='train_predictor-model/DNN: GPU id for deep learning')

    parser.add_argument(
        '--cores',
        action='store',
        dest='cores',
        default=4,
        type=int,
        help='Number of cores to be used')

    parser.add_argument(
        '--batchsize',
        action='store',
        dest='batch_size',
        type=int,
        default=10,
        help='train_predictor-model/DNN: batch size for deep learning')

    parsedArgs = parser.parse_args()

    print('Classification and parameter tuning requested..\n')
    if not parsedArgs.model:
        err = err + "\nNo classification model is specified"
    if (not os.access(parsedArgs.X, os.F_OK)):
        err = err + "\nError: Permission denied or could not find the X!"
        return err
    if (not os.access(parsedArgs.Y, os.F_OK)):
        err = err + "\nError: Permission denied or could not find the Y!"
        return err
    else:
        try:
            os.stat(parsedArgs.output_addr)
        except:
            os.mkdir(parsedArgs.output_addr)
            print(parsedArgs.output_addr, ' directory created')

    if not parsedArgs.data_name:
        parsedArgs.data_name = parsedArgs.X.split('/')[-1].split('.')[0]

    if parsedArgs.model == 'DNN':
        '''
            Deep learning
        '''
        arch = [
            int(layer) if float(layer) > 1 else float(layer)
            for layer in parsedArgs.dnn_arch.split(',')
        ]
        DNN_classifier(parsedArgs.X, parsedArgs.Y, arch,
                       parsedArgs.output_addr, parsedArgs.data_name,
                       parsedArgs.gpu_id, parsedArgs.epochs,
                       parsedArgs.batch_size)
    else:
        '''
            SVM and Random Forest
        '''
        if parsedArgs.model in ['SVM', 'RF']:
            classical_classifier(parsedArgs.X, parsedArgs.Y, parsedArgs.model,
                                 parsedArgs.output_addr, parsedArgs.data_name,
                                 parsedArgs.cores)
        else:
            return "\nNot able to recognize the model!"


if __name__ == '__main__':
    err = check_args(sys.argv)
    if err:
        print(err)
        exit()
