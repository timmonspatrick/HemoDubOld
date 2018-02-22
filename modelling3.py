# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 15:50:03 2017

@author: Patrick
"""

import numpy as np
from keras.models import Sequential
from keras.layers import Dense, Activation, Dropout
from sklearn.model_selection import train_test_split
from keras.callbacks import EarlyStopping, ModelCheckpoint
from keras import optimizers

inp = np.load("hemolytik_array.npy", allow_pickle=False)

X_train_org = inp[0:,1:-1]
Y_train_org = inp[0:,-1:]

X_train, X_val, Y_train, Y_val = train_test_split(X_train_org, Y_train_org, test_size=0.15, random_state=3)

dims = X_train.shape[1]

print(dims, 'dims')
print("Building model...")

nb_classes = Y_train.shape[1]
print(nb_classes, 'classes')

rate = 0.1

model = Sequential()
model.add(Dense(64, input_shape=(dims,), activation='relu'))
model.add(Dropout(0.5, noise_shape=None, seed=42))
model.add(Dense(64, init='uniform', activation='relu'))
model.add(Dropout(0.5, noise_shape=None, seed=42))
model.add(Dense(1, init='uniform', activation='sigmoid'))


sgd = optimizers.SGD(lr=0.01, decay=0.005, momentum=0.5, nesterov=True)
rmsprop = optimizers.RMSprop(lr=0.01, rho=0.9, epsilon=1e-08, decay=0.001)
adagrad = optimizers.Adagrad(lr=0.01, epsilon=1e-09, decay=0.0001)
adadelta = optimizers.Adadelta(lr=0.1, rho=0.95, epsilon=1e-08, decay=0.005)
adamax = optimizers.Adamax(lr=0.002, beta_1=0.9, beta_2=0.999, epsilon=1e-08, decay=0.001)
adam = optimizers.Adam(lr=0.01, beta_1=0.9, beta_2=0.999, epsilon=1e-08, decay=0.001)

model.compile(optimizer=adagrad, loss='binary_crossentropy', metrics=["accuracy"])
model.summary()

print("Model fitting...")

fBestModel = 'best_model2.h5' 
#early_stop = EarlyStopping(monitor='val_loss', patience=500, verbose=1) 
best_model = ModelCheckpoint(fBestModel, verbose=0, save_best_only=True)

model.fit(X_train, Y_train, validation_data = (X_val, Y_val), 
          epochs=30000, batch_size=dims, verbose=True,
          shuffle=True, 
          callbacks=[
                  best_model, 
                  #early_stop
                  ])
print("Model predicting...")
predictions = model.predict(X_train_org, batch_size=32, verbose=True)
score = model.evaluate(X_train_org, Y_train_org, batch_size=32, verbose=1)
#for a in range(len(predictions)):
#    print(str(predictions[a]) + "            "  + str(Y_train_org[a]))
print("Evaluation Score: "+ str(score))