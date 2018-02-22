"""
Utility used by the Network class to actually train.

Based on:
    https://github.com/fchollet/keras/blob/master/examples/mnist_mlp.py

"""
from keras.datasets import mnist, cifar10
from keras.models import Sequential
from keras.layers import Dense, Dropout
from keras.utils.np_utils import to_categorical
from keras.callbacks import EarlyStopping, ModelCheckpoint
from keras import optimizers
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
import numpy as np
import logging


# Helper: Early stopping.
early_stopper = EarlyStopping(patience=150)

logging.basicConfig(
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%m/%d/%Y %I:%M:%S %p',
    level=logging.DEBUG,
    filename='log9.txt'
)

optimizers_dict = {
        "rmsprop" : optimizers.rmsprop,
        "adam" : optimizers.adam,
        "sgd" : optimizers.SGD,
        "adagrad" : optimizers.adagrad,
        "adadelta" : optimizers.adadelta,
        "adamax" : optimizers.adamax,
        "nadam" : optimizers.nadam,
        }

def get_cifar10():
    """Retrieve the CIFAR dataset and process the data."""
    # Set defaults.
    nb_classes = 10
    batch_size = 64
    input_shape = (3072,)

    # Get the data.
    (x_train, y_train), (x_test, y_test) = cifar10.load_data()
    x_train = x_train.reshape(50000, 3072)
    x_test = x_test.reshape(10000, 3072)
    x_train = x_train.astype('float32')
    x_test = x_test.astype('float32')
    x_train /= 255
    x_test /= 255

    # convert class vectors to binary class matrices
    y_train = to_categorical(y_train, nb_classes)
    y_test = to_categorical(y_test, nb_classes)

    return (nb_classes, batch_size, input_shape, x_train, x_test, y_train, y_test)

def get_mnist():
    """Retrieve the MNIST dataset and process the data."""
    # Set defaults.
    nb_classes = 10
    batch_size = 128
    input_shape = (784,)

    # Get the data.
    (x_train, y_train), (x_test, y_test) = mnist.load_data()
    x_train = x_train.reshape(60000, 784)
    x_test = x_test.reshape(10000, 784)
    x_train = x_train.astype('float32')
    x_test = x_test.astype('float32')
    x_train /= 255
    x_test /= 255

    # convert class vectors to binary class matrices
    y_train = to_categorical(y_train, nb_classes)
    y_test = to_categorical(y_test, nb_classes)

    return (nb_classes, batch_size, input_shape, x_train, x_test, y_train, y_test)

def get_amps():
    nb_classes = 1
    
    inp = np.load("../hemolytik_array.npy", allow_pickle=False)

    X_train = np.array([i[2:-1] for i in inp if i[1] == 1])
    X_val = np.array([i[2:-1] for i in inp if i[1] == 0])
    X_fit = np.array([i[2:-1] for i in inp])

    Y_train = np.array([i[-1:] for i in inp if i[1] == 1])
    Y_val = np.array([i[-1:] for i in inp if i[1] == 0])

    
    ####Scaling#####
    scaler = StandardScaler()
    scaler.fit(X_fit)
    X_train = scaler.transform(X_train)
    X_val = scaler.transform(X_val)
    
    dims = X_train.shape[1]
    
    batch_size = dims
    
    input_shape = (dims,)

    return (nb_classes, batch_size, input_shape, X_train, X_val, Y_train, Y_val)

    

def compile_model(network, nb_classes, input_shape):
    """Compile a sequential model.

    Args:
        network (dict): the parameters of the network

    Returns:
        a compiled network.

    """
    # Get our network parameters.
    nb_layers = network['nb_layers']
    nb_neurons = network['nb_neurons']
    activation = network['activation']
    optimizer = network['optimizer']
    dropout = network['dropout']
    learning_rate = network['learning_rate']
    funnel = network['funnel']

    model = Sequential()
    
    # Add each layer.
    for i in range(nb_layers):

        # Need input shape for first layer.
        if i == 0:
            model.add(Dense(nb_neurons, activation=activation, input_shape=input_shape))
        else:
            model.add(Dense(nb_neurons, activation=activation))

        model.add(Dropout(dropout))  # hard-coded dropout
        
        if funnel:
            if nb_neurons // 2 == 0:
                nb_neurons = nb_neurons / 2
            elif nb_neurons // 2 == 1:
                nb_neurons = (nb_neurons / 2) + 1
        
    # Output layer.
    model.add(Dense(nb_classes, activation='sigmoid'))
    
    model.compile(optimizer=optimizers_dict[optimizer](lr=learning_rate),
                  metrics=['accuracy'],
                  loss='binary_crossentropy', 
                  )

    return model

def train_and_score(network, dataset):
    """Train the model, return test loss.

    Args:
        network (dict): the parameters of the network
        dataset (str): Dataset to use for training/evaluating

    """
    if dataset == 'amp':
        nb_classes, batch_size, input_shape, x_train, \
            x_test, y_train, y_test = get_amps()
    elif dataset == 'cifar10':
        nb_classes, batch_size, input_shape, x_train, \
            x_test, y_train, y_test = get_cifar10()
    elif dataset == 'mnist':
        nb_classes, batch_size, input_shape, x_train, \
            x_test, y_train, y_test = get_mnist()

    model = compile_model(network, nb_classes, input_shape,)
    
    #fBestModel = 'best_model.h5' 
    #best_model = ModelCheckpoint(fBestModel, verbose=0, save_best_only=True)

    
    model.fit(x_train, y_train,
              batch_size=batch_size,
              epochs=5000,  # using early stopping, so no real limit
              verbose=1,
              shuffle=True,
              validation_data=(x_test, y_test),
              callbacks=[early_stopper]
              )
    print(network)
    score = model.evaluate(x_test, y_test, verbose=1)
    print("Evaluation Score: "+ str(score))
    print("__________________")
    
    log_string = str(network) + "\n" + str(score) + "\n\n\n"
    logging.info(log_string)
    
    
    return score[1]  # 1 is accuracy. 0 is loss.
