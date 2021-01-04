# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 14:30:46 2020

@author: Daniel.Feeney
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 16:19:30 2020

@author: Daniel.Feeney
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import ctypes  

# Read in files
# only read .asc files for this work
fPath = 'C:/Users/Daniel.Feeney/Dropbox (Boa)/Snow Protocol/SnowboardLoadSol/'
entries = os.listdir(fPath)
fName = entries[1]
dat = pd.read_csv(fPath+fName,sep='         ', skiprows = 3, header = 0)
dat.columns = ['Time', 'LeftHeel', 'LeftMedial','LeftLateral','Total']
dat['Toes'] = dat.LeftMedial + dat.LeftLateral

fig1, ax1 = plt.subplots()
ax1.plot(dat.Total[4000:7500])
ax1.plot(dat.LeftHeel[4000:7500])
ax1.plot(dat.Toes[4000:7500])
plt.legend()
print('Select 2 points: start and end of analysis trial')
pts = np.asarray(plt.ginput(2, timeout=-1))

# downselect the region of the dataframe you'd like
dat = dat.iloc[int(np.floor(pts[0,0])) : int(np.floor(pts[1,0])),:]

#### Initial analysis plans: segment heel and toe turns, calcualte smoothness
## as the CV or SD of force during the turn, calcualte turn time, calculate symmetry

# Find indices of toe and heel turns but leave original data untouched 
heelThresh = 150 #below this value will be set to 0 temporarily to find indices to start/end turns
toeThresh = 75

tmpToes = np.array(dat.Toes)
tmpHeel = np.array(dat.LeftHeel)

tmpHeel[tmpHeel < heelThresh] = 0
tmpToes[tmpToes < toeThresh] = 0

# Select toe turn starts
plt.plot(tmpToes)
toestart = np.asarray(plt.ginput(5, timeout=-1))
realToeStart = [int(toestart[0,0]), int(toestart[1,0]), (toestart[2,0]), int(toestart[3,0]),int(toestart[4,0])]

# Select heel turn starts
plt.plot(tmpHeel)
heelstart = np.asarray(plt.ginput(5, timeout=-1))
realHeelStart = [int(heelstart[0,0]), int(heelstart[1,0]), (heelstart[2,0]), int(heelstart[3,0]),int(heelstart[4,0])]
   
    

