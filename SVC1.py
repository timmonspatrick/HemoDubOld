# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 17:32:38 2018

@author: Patrick
"""

from sklearn.svm import SVC, LinearSVC # "Support vector classifier"
from sklearn.model_selection import train_test_split
import numpy as np
from sklearn.preprocessing import StandardScaler

inp = np.load("hemolytik_array.npy", allow_pickle=False)

X_train_org = inp[0:,1:-1]
Y_train_org = inp[0:,-1:]

###Scaling#####
scaler = StandardScaler()
scaler.fit(X_train_org)
X_train_org = scaler.transform(X_train_org)
#######

X_train, X_val, Y_train, Y_val = train_test_split(X_train_org, Y_train_org, test_size=0.10, random_state=31)

Y_train = [i[0] for i in Y_train]
Y_val = [i[0] for i in Y_val]

model = SVC(kernel='poly', C=1E10)
#model = LinearSVC(random_state=0, C=0.001)
model.fit(X_train, Y_train)

predictions = model.predict(X_val)
predictions = list(predictions)

accurate_count = 0
for i in range(len(predictions)):
    print(str(predictions[i]) + "     " + str(Y_val[i]) )
    if predictions[i] == Y_val[i]:
        accurate_count = accurate_count + 1
        
print("Accuracy: " + str(accurate_count/len(predictions)))