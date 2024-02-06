import sys
from classifier.DNN import DNNMutliclass16S
from utility.file_utility import FileUtility
import matplotlib.pyplot as plt
import matplotlib
from utility.visualization_utility import create_mat_plot

#[latex_line, p_micro, r_micro, f1_micro, p_macro, r_macro, f1_macro, history]=FileUtility.load_obj('/Users/Prime/Documents/MicroPheno/449_kmer/output_dir/nn_classification_results_test_crohn_mlp_1000_0.61.pickle')
[
    latex_line, p_micro, r_micro, f1_micro, p_macro, r_macro, f1_macro, history
] = FileUtility.load_obj(
    '/Users/Prime/Documents/MicroPheno/449_kmer/output_dir/nn_classification_results_test_crohn_mlp_1024-0.2-512-0.1-128-64_0.63.pickle'
)

(loss_values, val_loss_values, epochs) = history

print('p_micro')
print(p_micro)
print('r_micro')
print(r_micro)
print('f1_micro')
print(f1_micro)

print('p_macro')
print(p_macro)
print('r_macro')
print(r_macro)
print('f1_macro')
print(f1_macro)

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
matplotlib.rcParams["axes.edgecolor"] = "black"
matplotlib.rcParams["axes.linewidth"] = 0.6
plt.plot(epochs, loss_values, 'ro', label='Loss for train set')
plt.plot(epochs, val_loss_values, 'b', label='Loss for test set')
plt.xlabel('Epochs')
plt.ylabel('Loss')
plt.legend(
    loc=1,
    prop={'size': 8},
    ncol=1,
    edgecolor='black',
    facecolor='white',
    frameon=True)
plt.title('Loss with respect to the number of epochs for train and test sets')
plt.savefig("449_1024.png")
print('done')
