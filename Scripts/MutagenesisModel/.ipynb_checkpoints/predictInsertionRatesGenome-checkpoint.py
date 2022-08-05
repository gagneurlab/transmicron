##############################################################################
# load packages
##############################################################################

import rpy2.robjects as robjects
import numpy as np
import joblib
import pandas as pd
import sklearn
from sklearn.ensemble import RandomForestClassifier
import csv
import h5py
import tables 

##############################################################################
# prepare input
##############################################################################

#read model 
read_path= snakemake.input["predict_insertion_model"]
trained_model = joblib.load(open(read_path, 'rb'))

#read table with all TTAA / TA sites in the genome scored on all features
path_to_input = snakemake.input["step1_input_negatives_whole_genome"]
h5file = tables.open_file(path_to_input)

#read the table of positives used for training the mutagenesis model
positives_whole_table = robjects.r['readRDS'](snakemake.input["step1_input_positives"])


def read_R(raw_frame):
    
    colnames = raw_frame.colnames
    pandas_frame = pd.DataFrame(np.transpose(np.asarray(raw_frame)))
    pandas_frame.columns = colnames
    
    return pandas_frame
positives_whole_table = read_R(positives_whole_table)

##############################################################################
# predict insertion rates
##############################################################################

#### we need the sharae of TTAA / TA sites used for training the mutagenesis modeöamong all TTAA / TA sites in the genome 

#calculate total number of TA / TTAA sites in the genome
numberNegatives = 0
for i in h5file.root:
    numberNegatives = numberNegatives+h5file.get_node(i).shape[0]

#calculate fraction of training set
numberPositives = positives_whole_table.shape[0] * 0.9 # number of negatives equals number of positives; 10% are used for testing.
shareUsedForTraining = numberPositives / numberNegatives

#### use the trained model to predict insertion rates (iteratively over chromosomes)

predictions = []
for i in h5file.root:
        chromosome = pd.DataFrame(h5file.get_node(i)[0:h5file.get_node(i).shape[0]]) #subset hdf5 file for chromosomes
        predI=(np.log(1+(trained_model.predict_proba(chromosome)[:, 1]/trained_model.predict_proba(chromosome)[:, 0])*shareUsedForTraining)) #caclulate insertion rates
        
        #in borderline cases, predicitions will result in INF values; replace these with 1
        flat=np.ravel(predI)
        if(max(flat)>1):
            predI[np.ravel(predI)>1]=1
        predictions.append(predI)
        
#format output        
predictions = np.transpose(np.matrix(np.concatenate(predictions, axis=0)))

##############################################################################
# save output
##############################################################################

np.savetxt(snakemake.output[0], predictions, delimiter=",")

