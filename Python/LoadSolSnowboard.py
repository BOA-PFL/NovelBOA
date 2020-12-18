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
dat['Toes'] = dat.LeftMedial + dat.LeftLateral


plt.plot(dat.Total[4000:7500])
plt.plot(dat.LeftHeel[4000:7500])
plt.plot(dat.Toes[4000:7500])
plt.legend()
print('Select start and end of analysis trial')
pts = np.asarray(plt.ginput(2, timeout=-1))

# downselect the region of the dataframe you'd like
dat = dat.iloc[int(np.floor(pts[0,0])) : int(np.floor(pts[1,0])),:]

#### Initial analysis plans: segment heel and toe turns, calcualte smoothness
## as the CV or SD of force during the turn, calcualte turn time, calculate symmetry

# Find indices of toe and heel turns but leave original data untouched 
fThresh = 150 #below this value will be set to 0.
stepLen = 45 #Set value to look forward 
tmpToes = np.array(dat.Toes)
tmpHeel = np.array(dat.LeftHeel)

tmpHeel[tmpHeel < fThresh] = 0

def findTurnStart(force):
    lic = []
    for step in range(len(force)-1):
        if force[step] <= fThresh and np.mean(force[step:step+10]) >= fThresh:
            lic.append(step)
    return lic

heelStart = findTurnStart(tmpHeel) #gets too many entries, trim below
realHeelStart = []
# crass way to crim the indices too close together
for numel, entry in enumerate(heelStart):
    if heelStart[numel + 1] - entry > 50:
        realHeelStart.append(entry)

#Find takeoff from FP when force goes from above thresh to 0
def findTakeoffs(force):
    lto = []
    for step in range(len(force)-1):
        if force[step] >= fThresh and force[step + 1] == 0:
            lto.append(step + 1)
    return lto