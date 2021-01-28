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
fPath = 'C:/Users/Daniel.Feeney/Dropbox (Boa)/Snow Protocol/SnowboardProtocol/'
entries = os.listdir(fPath)
fName = entries[1]
dat = pd.read_csv(fPath+fName,sep='         ', skiprows = 3, header = 0, index_col = False)
#dat.columns = ['Time', 'LeftHeel', 'LeftMedial','LeftLateral','Total']
dat.columns = ['Time', 'LHeel', 'LMedial','LLateral','LTotal', 'Time2', 'RLateral','RMedial','RHeel','RTotal']

dat['LToes'] = dat.LMedial + dat.LLateral
dat['RToes'] = dat.RMedial + dat.RLateral

fig, ax = plt.subplots()
ax.plot(dat.LTotal, label = 'Left Total Force')
ax.plot(dat.RTotal, label = 'Right Total Force')
fig.legend()
print('Select start and end of analysis trial')
pts = np.asarray(plt.ginput(2, timeout=15))
plt.close()
# downselect the region of the dataframe you'd like
dat = dat.iloc[int(np.floor(pts[0,0])) : int(np.floor(pts[1,0])),:]

#fig, ax = plt.subplots()
#ax.plot(dat.LToes, label='Left Toes')
#ax.plot(dat.LHeel, label='Left Heel')
#fig.legend()

#### Initial analysis plans: segment heel and toe turns, calcualte smoothness
## as the CV or SD of force during the turn, calcualte turn time, calculate symmetry

# Find indices of toe and heel turns but leave original data untouched 
heelThresh = 150 #below this value will be set to 0 temporarily to find indices to start/end turns
toeThresh = 75

tmpToes = np.array(dat.LToes)
tmpHeel = np.array(dat.LHeel)

tmpHeel[tmpHeel < heelThresh] = 0
tmpToes[tmpToes < toeThresh] = 0

# Select toe turn starts
plt.plot(tmpToes)
toestart = np.asarray(plt.ginput(5, timeout=-1))
realToeStart = [int(toestart[0,0]), int(toestart[1,0]), (toestart[2,0]), int(toestart[3,0]),int(toestart[4,0])]
plt.close()

# Select heel turn starts
plt.plot(tmpHeel)
heelstart = np.asarray(plt.ginput(5, timeout=-1))
realHeelStart = [int(heelstart[0,0]), int(heelstart[1,0]), (heelstart[2,0]), int(heelstart[3,0]),int(heelstart[4,0])]
plt.close()   

###### After finding starts of turns, find avg, SD, CV, etc. for each turn ####

plt.plot(dat.LToes[realToeStart[0]:realToeStart[1]], label='L Toes')
plt.plot(dat.RToes[realToeStart[0]:realToeStart[1]], label = 'R Toes')
plt.legend()
