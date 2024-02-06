# Written By Arun Manoharan
import numpy

from sklearn.preprocessing import LabelEncoder
from classifier.DNN import DNNMutliclass16S
from utility.file_utility import FileUtility

#file_name = sys.argv[1]
#prediction = sys.argv[2]
#X_file = sys.argv[3]
#Y_file = sys.argv[4]

file_name = '/home/arun/bioinformatics/deep_learning/449/body_habitat/nn_classification_results_fourfournine_layers_mlp_1024-0.2-512-0.1-128-64_0.85.pickle'
predict_X = '/home/arun/bioinformatics/eugene/kmers/eugene_6-mers_1000.npz'
train_X = '/home/arun/bioinformatics/449/kmers/449_6-mers_1000.npz'
train_Y = '/home/arun/bioinformatics/deep_learning/449/body_habitat/c.txt'

X = FileUtility.load_sparse_csr(train_X).toarray()

print(X)
print(X.shape)
# labels
Y = FileUtility.load_list(train_Y)
encoder = LabelEncoder()
encoder.fit(Y)
Y_encoded = encoder.transform(Y)
Y_onehot = np_utils.to_categorical(Y_encoded)

DNN = DNNMutliclass16S(X, Y, model_arch=[1024, 0.2, 512, 0.1, 128, 64])
#model = DNN.get_pretrained_model(file_name, trainable =False)

model = DNN.get_MLP_model()
model.fit(X, Y_onehot, epochs=20, batch_size=100, verbose=0)

predict_X = FileUtility.load_sparse_csr(predict_X).toarray()
pred = model.predict_classes(predict_X)

print(pred)
print(encoder.inverse_transform(pred[0]))

un_onehot_prediction = numpy.argmax(predict_X, axis=None, out=None)
un_encoded_prediction = encoder.inverse_transform(un_onehot_prediction)

print(un_encoded_prediction)
