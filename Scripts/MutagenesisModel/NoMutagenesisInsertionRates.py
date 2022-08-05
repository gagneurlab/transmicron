# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.10.2
#   kernelspec:
#     display_name: Python [conda env:anaconda-bredthauer_transposon]
#     language: python
#     name: conda-env-anaconda-bredthauer_transposon-py
# ---

# +
import rpy2.robjects as robjects
import numpy as np
import matplotlib.pyplot as plt
import sklearn

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report, confusion_matrix
  
import pickle
from sklearn.metrics import roc_curve, roc_auc_score

import pandas as pd
import matplotlib.pyplot as plt
import sklearn
from sklearn.datasets import make_classification
from sklearn.metrics import PrecisionRecallDisplay
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.datasets import make_classification
from sklearn.model_selection import cross_val_score
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.preprocessing import OneHotEncoder

import csv
import h5py
import tables   # but in this tutorial we use "from tables import \*"
import numpy as np

path_to_input = snakemake.input["step1_input_negatives_whole_genome"]


#h5file = tables.open_file("/s/project/transposon/Output/Snakemake/DLBCL_PB/step1_input_negatives_whole_genome.h5")
h5file = tables.open_file(path_to_input)


sum_rows = 0
for i in h5file.root:
        length_chrom_i = h5file.get_node(i).shape[0]
        sum_rows= sum_rows +length_chrom_i
ones = np.array([1])
predictions=np.transpose(np.matrix((np.repeat(ones, sum_rows, axis=0))))

np.savetxt(snakemake.output[0], predictions, delimiter=",")

