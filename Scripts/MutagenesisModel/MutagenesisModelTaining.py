##################################################
# description
##################################################

# this script is used to train the mutagenesis model

##############################################################
# load packages
##############################################################

import rpy2.robjects as robjects
import numpy as np  
import pandas as pd
import sklearn
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import RandomizedSearchCV
import csv
import joblib

##############################################################
# read input
##############################################################
training_data = robjects.r['readRDS'](snakemake.input[0])

def read_R(raw_frame):
    
    colnames = raw_frame.colnames
    pandas_frame = pd.DataFrame(np.transpose(np.asarray(raw_frame)))
    pandas_frame.columns = colnames
    
    return pandas_frame

training_data = read_R(training_data)

trainy = training_data.iloc[:,0]
trainx= training_data.iloc[:,1:]

##############################################################
# hyperparameter optimization (unused)
##############################################################

#random_grid = {'bootstrap': [True, False],
#              'max_depth': [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, None],
#              'max_features': ['auto', 'sqrt'],
#              'min_samples_leaf': [1, 2, 4],
#              'min_samples_split': [2, 5, 10],
#              'n_estimators': [5,10,20,40,80,100,200,400,800]
#}

#rf = RandomForestClassifier(n_estimators = 10, random_state = 42)
#rf = RandomizedSearchCV(estimator = rf, param_distributions = random_grid, n_iter = 100, cv = 5, verbose=2, random_state=42, n_jobs = -1)

##############################################################
# training the model
##############################################################
rf = RandomForestClassifier(bootstrap=False, max_depth=70, min_samples_leaf=4,
                       min_samples_split=10, n_estimators=100, random_state=42)

rf.fit(trainx, trainy)

##############################################################
# save the model
############################################################## 
#pickle.dump(rf, open(snakemake.output[0], 'wb'))
joblib.dump(rf, snakemake.output[0], 9)

