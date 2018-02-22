# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 15:50:03 2017
98.5% Accuracy Reached
@author: Patrick
"""

import numpy as np
from keras.models import Sequential
from keras.layers import Dense, Activation, Dropout, Conv1D, Embedding
from sklearn.model_selection import train_test_split
from keras.callbacks import EarlyStopping, ModelCheckpoint
from keras import optimizers
from sklearn.preprocessing import StandardScaler, normalize

inp = np.load("hemolytik_array_doubled.npy", allow_pickle=False)

X_train = np.array([i[2:] for i in inp if i[1] == 1])
X_val = np.array([i[2:] for i in inp if i[1] == 0])
X_fit = np.array([i[2:] for i in inp])

Y_train = np.array([i[-1:] for i in inp if i[1] == 1])
Y_val = np.array([i[-1:] for i in inp if i[1] == 0])


####Scaling#####
scaler = StandardScaler()
scaler.fit(X_fit)
X_train = scaler.transform(X_train)
X_val = scaler.transform(X_val)
#######

###Normalise#####
#normalize(X_train_org, norm = 'l1', axis = 1, copy = False, )

dims = X_train.shape[1]

print(dims, 'dims')
print("Building model...")

nb_classes = Y_train.shape[1]
print(nb_classes, 'classes')

rate = 0.1

model = Sequential()
model.add(Dense(88, input_shape=(dims,), init='uniform', activation='tanh'))
model.add(Dropout(0.87, noise_shape=None, seed=42))
model.add(Dense(88, init='uniform', activation='tanh'))
model.add(Dropout(0.87, noise_shape=None, seed=42))
model.add(Dense(88, init='uniform', activation='tanh'))
model.add(Dropout(0.87, noise_shape=None, seed=42))

model.add(Dense(1, init='uniform', activation='sigmoid'))


sgd = optimizers.SGD(lr=0.01, decay=0.005, momentum=0.5, nesterov=True)
rmsprop = optimizers.RMSprop(lr=0.01, rho=0.9, epsilon=1e-08, decay=0.001)
rmsprop = optimizers.RMSprop(lr=0.01)
adagrad = optimizers.Adagrad(lr=0.005)
adadelta = optimizers.Adadelta(lr=0.1, rho=0.95, epsilon=1e-08, decay=0.005)
adamax = optimizers.Adamax(lr=0.002, beta_1=0.9, beta_2=0.999, epsilon=1e-08, decay=0.001)
adam = optimizers.Adam(lr=0.01, beta_1=0.9, beta_2=0.999, epsilon=1e-08, decay=0.001)

model.compile(optimizer=rmsprop, loss='binary_crossentropy', metrics=["accuracy"])
model.summary()

print("Model fitting...")

fBestModel = 'best_modelxd.h5' 
#early_stop = EarlyStopping(monitor='val_loss', patience=500, verbose=1) 
best_model = ModelCheckpoint(fBestModel, verbose=0, save_best_only=True)

model.fit(X_train, Y_train, validation_data = (X_val, Y_val), 
          epochs=20000, batch_size=dims, verbose=True,
          shuffle=True, 
          callbacks=[
                  best_model, 
                  #early_stop
                  ])
print("Model predicting...")
predictions = model.predict(X_train_org, batch_size=32, verbose=True)
score = model.evaluate(X_val, Y_val, batch_size=32, verbose=1)
#for a in range(len(predictions)):
#    print(str(predictions[a]) + "            "  + str(Y_train_org[a]))
print("Evaluation Score: "+ str(score))