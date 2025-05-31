import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
from matplotlib import pylab
import seaborn as sns 
from sklearn.metrics.pairwise import cosine_similarity
from statannot import add_stat_annotation
import matplotlib as mpl
from scipy import stats, cluster
import glob
import re
from sklearn.metrics import mean_squared_error, r2_score
from statsmodels.stats import multitest

import warnings
warnings.simplefilter('ignore')

## load_large_dataFrame
def load_large_dataFrame(input_file, sep=",", header=0, index_col=0, chunksize=100000, compressed=False):
    if compressed:
        TextFileReader = pd.read_csv(input_file, chunksize=chunksize, sep=sep, header=header,index_col=index_col, compression='gzip')
    else:
        TextFileReader = pd.read_csv(input_file, chunksize=chunksize, sep=sep, header=header,index_col=index_col)
    dfList=[]
    for df in TextFileReader:
        dfList.append(df)
    final_df = pd.concat(dfList,sort=False)
    return final_df

import numpy as np
import os
import shutil
from sklearn.preprocessing import QuantileTransformer, quantile_transform

if os.path.exists("all_raw_mean_qn"):
    pass
else:
    os.makedirs("all_raw_mean_qn")


files = glob.glob("all_raw_mean/*")
for file in files:
    print (file)
    s = file.split("/")[1]
    print (s)
    tpm_df = load_large_dataFrame(file, sep="\t")
    qt_tpm = quantile_transform(tpm_df, n_quantiles=10, random_state=0, axis=0, copy=True)
    qt_tpm_df = pd.DataFrame(qt_tpm, index=tpm_df.index, columns=tpm_df.columns)
    qt_tpm_df.to_csv("all_raw_mean_qn/"+s.strip("txt")+"quantile_transform.csv") 
