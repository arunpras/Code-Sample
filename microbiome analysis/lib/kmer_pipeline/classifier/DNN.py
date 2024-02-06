__author__ = "Arun Manoharan"
__license__ = "Apache 2"
__version__ = "2.0.0"
__maintainer__ = "Arun Manoharan"
__email__ = "arun@primediscoveries.com"

import sys

sys.path.append('../')

from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation
from keras.layers import Convolution2D, MaxPooling2D, Flatten, Dense
from keras.utils import np_utils
from sklearn.preprocessing import LabelEncoder
from utility.file_utility import FileUtility
import numpy as np
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import f1_score, precision_score, recall_score
import os
import matplotlib.pyplot as plt
import matplotlib


class DNNMutliclass16S(object):
    '''
    Deep MLP Neural Network
    '''

    def __init__(self, X, Y, model_arch=[500]):
        # rep. X
        self.X = X
        # encoding Y
        self.Y = Y
        self.C = len(set(Y))
        encoder = LabelEncoder()
        encoder.fit(Y)
        self.encoded_Y = encoder.transform(Y)
        self.onehot_y = np_utils.to_categorical(self.encoded_Y)
        # model's arch
        self.model_arch = model_arch

    def get_MLP_model(self):
        '''
        Create the model
        :return:
        '''
        # creating the model
        model = Sequential()
        for layer_idx, h_layer_size in enumerate(self.model_arch):
            if layer_idx == 0:
                model.add(Dense(h_layer_size, input_dim=self.X.shape[1], activation='relu'))
            else:
                if h_layer_size < 1:
                    model.add(Dropout(h_layer_size))
                else:
                    model.add(Dense(h_layer_size, activation='relu'))
        model.add(Dense(self.C, activation='softmax'))
        # Compile model
        model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
        return model

    def get_pretrained_model(self, file_name, trainable):
        pretrained_weights = FileUtility.load_obj(file_name)

        h_sizes = [float(x) for x in file_name.split('/')[-1].split('_')[3].split('-')]
        model = Sequential()
        for layer_idx, h_layer_size in enumerate(h_sizes):
            if layer_idx == 0:
                model.add(Dense(int(h_layer_size), input_dim=self.X.shape[1], weights=pretrained_weights[0],
                                activation='relu', trainable=trainable))
            else:
                if h_layer_size < 1:
                    model.add(Dropout(h_layer_size, weights=pretrained_weights[layer_idx], trainable=trainable))
                else:
                    model.add(Dense(int(h_layer_size), weights=pretrained_weights[layer_idx], activation='relu',
                                    trainable=trainable))
        if self.model_arch:
            for layer_idx, h_layer_size in enumerate(self.model_arch):
                if h_layer_size < 1:
                    model.add(Dropout(h_layer_size))
                else:
                    model.add(Dense(h_layer_size, activation='relu'))
        model.add(Dense(self.C, activation='softmax'))
        model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
        return model

    def cross_validation(self, result_filename, gpu_dev='2', n_fold=5, epochs=50, batch_size=100, model_strct='mlp',
                         pretrained_model=False, trainable=False):
        '''
        :param result_filename:
        :param gpu_dev:
        :param n_fold:
        :param epochs:
        :param batch_size:
        :param model_strct:
        :param pretrained_model:
        :param trainable:
        :return:
        '''
        os.environ["CUDA_VISIBLE_DEVICES"] = gpu_dev

        skf = StratifiedKFold(n_splits=n_fold, shuffle=True)

        p_micro = []
        p_macro = []
        r_micro = []
        r_macro = []
        f1_micro = []
        f1_macro = []

        for train_index, valid_index in skf.split(self.X, self.Y):
            print ('\n Evaluation on a new fold is now get started ..')
            X_train = self.X[train_index, :]
            y_train = self.onehot_y[train_index, :]
            y_class_train = self.encoded_Y[train_index]
            X_valid = self.X[valid_index, :]
            y_valid = self.onehot_y[valid_index, :]
            y_class_valid = self.encoded_Y[valid_index]

            if pretrained_model:
                model = self.get_pretrained_model(model_strct, trainable)
            else:
                if model_strct == 'mlp':
                    model = self.get_MLP_model()

            # fitting
            history = model.fit(X_train, y_train, epochs=epochs, batch_size=batch_size, shuffle=True,
                                validation_data=(X_valid, y_valid), verbose=0)
            pred = model.predict_classes(X_valid)
            # score-calculations
            f1_micro.append(f1_score(y_class_valid, pred, average='micro'))
            f1_macro.append(f1_score(y_class_valid, pred, average='macro'))
            p_micro.append(precision_score(y_class_valid, pred, average='micro'))
            p_macro.append(precision_score(y_class_valid, pred, average='macro'))
            r_micro.append(recall_score(y_class_valid, pred, average='micro'))
            r_macro.append(recall_score(y_class_valid, pred, average='macro'))

        # mean values
        f1mac = np.mean(f1_macro)
        f1mic = np.mean(f1_micro)
        prmac = np.mean(p_macro)
        prmic = np.mean(p_micro)
        remac = np.mean(r_macro)
        remic = np.mean(r_micro)
        # std values
        sf1mac = np.std(f1_macro)
        sf1mic = np.std(f1_micro)
        sprmac = np.std(p_macro)
        sprmic = np.std(p_micro)
        sremac = np.std(r_macro)
        sremic = np.std(r_micro)
        # table
        latex_line = ' & '.join([str(np.round(x, 2)) + ' $\\pm$ ' + str(np.round(y, 2)) for x, y in
                                 [[prmic, sprmic], [remic, sremic], [f1mic, sf1mic], [prmac, sprmac], [remac, sremac],
                                  [f1mac, sf1mac]]])

        print (latex_line)

        history_dict = history.history
        loss_values = history_dict['loss']
        val_loss_values = history_dict['val_loss']
        epochs = range(1, len(loss_values) + 1)

        '''
        Saving the results
        '''
        if pretrained_model:
            model_strct = 'pretrained'
            # print (model.summary())

        output_file_name = '_'.join([result_filename, model_strct, '-'.join([str(x) for x in self.model_arch]),str(np.round(f1mac, 2))])
        output_file_name += '_metrics.txt'
        print(output_file_name)
        metrics_file = open(output_file_name, 'w')

        metrics_file.write('loss_values\t' + str(loss_values)+'\n')
        metrics_file.write('validation_loss_values\t' + str(val_loss_values)+'\n')
        metrics_file.write('num_epochs\t' + str(epochs)+'\n')
        metrics_file.write('\n')
        metrics_file.write('p_micro\t' + str(p_micro)+'\n')
        metrics_file.write('p_macro\t' + str(p_macro)+'\n')
        metrics_file.write('r_micro\t' + str(r_micro)+'\n')
        metrics_file.write('r_macro\t' + str(r_macro)+'\n')
        metrics_file.write('f1_micro\t' + str(f1_micro)+'\n')
        metrics_file.write('f1_macro\t' + str(f1_macro)+'\n')


        weights = []
        for x in model.layers:
            weights.append(x.get_weights())

        '''
        Saving the parameters and weights
        '''

        FileUtility.save_obj('_'.join(
            [result_filename, 'layers', model_strct, '-'.join([str(x) for x in self.model_arch]),
             str(np.round(f1mac, 2))]), weights)

    def save_results(self, metrics, output_file):
        pass

    @staticmethod
    def load_history(filename, fileout):
        '''
        Plot the history
        :param filename:
        :param fileout:
        :return:
        '''
        [latex_line, p_micro, r_micro, f1_micro, p_macro, r_macro, f1_macro, history] = FileUtility.load_obj(filename)
        (loss_values, val_loss_values, epochs) = history
        matplotlib.rcParams['mathtext.fontset'] = 'stix'
        matplotlib.rcParams['font.family'] = 'STIXGeneral'
        matplotlib.rcParams['mathtext.fontset'] = 'custom'
        matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
        matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
        matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
        matplotlib.rcParams["axes.edgecolor"] = "black"
        matplotlib.rcParams["axes.linewidth"] = 0.6
        plt.rc('text', usetex=True)
        plt.plot(epochs, loss_values, 'ro', label='Loss for train set')
        plt.plot(epochs, val_loss_values, 'b+', label='Loss for test set')
        plt.xlabel('Epochs')
        plt.ylabel('Loss')
        plt.legend(loc=1, prop={'size': 8}, ncol=1, edgecolor='black', facecolor='white', frameon=True)
        plt.title('Loss with respect to the number of epochs for train and test sets')
        plt.savefig(fileout + '.pdf')
        plt.show()

    @staticmethod
    def make_activation_function(file_name, X, last_layer=None):
        pretrained_weights = FileUtility.load_obj(file_name)
        if last_layer:
            h_sizes = [float(x) for x in file_name.split('/')[-1].split('_')[-2].split('-')] + [last_layer]
        else:
            h_sizes = [float(x) for x in file_name.split('/')[-1].split('_')[-2].split('-')]
        model = Sequential()
        for layer_idx, h_layer_size in enumerate(h_sizes):
            if layer_idx == 0:
                model.add(
                    Dense(int(h_layer_size), input_dim=X.shape[1], weights=pretrained_weights[0], activation='relu'))
            else:
                if h_layer_size < 1:
                    model.add(Dropout(h_layer_size, weights=pretrained_weights[layer_idx]))
                else:
                    if layer_idx == len(h_sizes) - 1 and last_layer:
                        model.add(Dense(int(h_layer_size), weights=pretrained_weights[layer_idx], activation='softmax'))
                    else:
                        model.add(Dense(int(h_layer_size), weights=pretrained_weights[layer_idx], activation='relu'))
        activations = model.predict(X)
        np.savetxt(file_name.replace(file_name.split('/')[-1].split('_')[0], 'activationlayer'), activations)
        return activations

    @staticmethod
    def result_visualization(filename):
        [latex_line, p_micro, r_micro, f1_micro, p_macro, r_macro, f1_macro,
         (loss_values, val_loss_values, epochs)] = FileUtility.load_obj(filename)
        print(latex_line)
