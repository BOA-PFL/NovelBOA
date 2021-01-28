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
fName = entries[4]
dat = pd.read_csv(fPath+fName,sep='         ', skiprows = 3, header = 0, index_col = False)
#dat.columns = ['Time', 'LeftHeel', 'LeftMedial','LeftLateral','Total']
dat.columns = ['Time', 'LHeel', 'LMedial','LLateral','LTotal', 'Time2', 'RLateral','RMedial','RHeel','RTotal']
dat['LToes'] = dat.LMedial + dat.LLateral
dat['RToes'] = dat.RMedial + dat.RLateral


plt.plot(dat.LTotal)
plt.plot(dat.LHeel)
plt.plot(dat.LToes)
plt.legend()
print('Select start and end of analysis trial')
pts = np.asarray(plt.ginput(2, timeout=-1))

# downselect the region of the dataframe you'd like
dat = dat.iloc[int(np.floor(pts[0,0])) : int(np.floor(pts[1,0])),:]

#### Initial analysis plans: segment heel and toe turns, calcualte smoothness
## as the CV or SD of force during the turn, calcualte turn time, calculate symmetry

# Find indices of toe and heel turns but leave original data untouched 
heelThresh = 150 #below this value will be set to 0 temporarily to find indices to start/end turns
stepLen = 45 #Set value to look forward 
tmpToes = np.array(dat.LToes)
tmpHeel = np.array(dat.LHeel)

tmpHeel[tmpHeel < heelThresh] = 0

def findTurnStart(force, thresh):
    ## enter np array and a threshold value to detect 0s ##
    lic = []
    for step in range(len(force)-1):
        if force[step] <= thresh and np.mean(force[step:step+10]) >= thresh:
            lic.append(step)
    return lic

def trimEntries(zeroEntries, lenForward):
    realTurnStart = []
    # crass way to crim the indices too close together
    for numel, entry in enumerate(zeroEntries):
        if numel < (len(zeroEntries) - 1):
            if zeroEntries[numel + 1] - entry > lenForward:
                realTurnStart.append(entry)
    return realTurnStart

#### Use functions above to find the starts of turns for heel and toes
heelStart = findTurnStart(tmpHeel, 150) #gets too many entries, trim below
realHeelStart = trimEntries(heelStart, 150)

# Same will be done with toes, set below thresh to 0, find indices
toeThresh = 150
tmpToes[tmpToes < toeThresh] = 0

toeStart = findTurnStart(tmpToes, 100)
realToeStart = trimEntries(toeStart, 100)
plt.plot(tmpToes)
for xc in realToeStart:
    plt.axvline(x = xc)

correctStarts = ctypes.windll.user32.MessageBoxW(0, "Are the starts correct?", "Heel Start", 3)
if correctStarts == 7:
    print('Select new points')
    realToeStart = []
    plt.plot(tmpToes)
    toestart2 = np.asarray(plt.ginput(5, timeout=-1))


    
plt.plot(tmpHeel)
for xc in realHeelStart:
    plt.axvline(x = xc, color = 'r')
    
   

    

