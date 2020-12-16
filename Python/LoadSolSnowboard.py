# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 16:19:30 2020

@author: Daniel.Feeney
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# Read in files
# only read .asc files for this work
fPath = 'C:/Users/Daniel.Feeney/Dropbox (Boa)/Snow Protocol/SnowboardLoadSol/'
entries = os.listdir(fPath)
fName = entries[1]
dat = pd.read_csv(fPath+fName,sep='         ', skiprows = 3, header = 0)
dat.columns = ['Time', 'LeftHeel', 'LeftMedial','LeftLateral','Total']


plt.plot(dat.Total[4000:7500])
plt.plot(dat.LeftHeel[4000:7500])
plt.plot(dat.LeftMedial[4000:7500])
plt.plot(dat.LeftLateral[4000:7500])
plt.legend()