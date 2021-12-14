# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 10:55:52 2020
# This relies on plotting based on heel strike approximations. For run data,
# this may not be a good enough approximation, so using the plotting from peaks 
# is recommended. 
@author: Daniel.Feeney
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns

#%matplotlib qt
# Read in files
# only read .asc files for this work
fPath = 'C:/Users/Daniel.Feeney/Dropbox (Boa)/EndurancePerformance/Altra_MontBlanc_Jan2021/PedarPressures/'
fileExt = r".asc"
entries = [fName for fName in os.listdir(fPath) if fName.endswith(fileExt)]

### in order to plot same time after landing, need to specify HS, Mid Stance
### and toe off times here. For running, 4, 8, and 12 at 50 Hz suggested. For walking
### 5, 20, and 30 are a good start at 50 Hz
hs = 4
ms = 8
to = 12

# Define constants and options
fThresh = 700 #below this value will be set to 0.
# list of functions 
# finding landings on the force plate once the filtered force exceeds the force threshold
def findLandings(force):
    lic = []
    for step in range(len(force)-1):
        if force[step] == 0 and force[step + 1] >= fThresh:
            lic.append(step)
    return lic

#Find takeoff from FP when force goes from above thresh to 0
def findTakeoffs(force):
    lto = []
    for step in range(len(force)-1):
        if force[step] >= fThresh and force[step + 1] == 0:
            lto.append(step + 1)
    return lto

def reshapeArray(listIn):

    listIn.insert(0,0)
    listIn.insert(6,0)
    listIn.insert(97,0)
    listIn.insert(98,0)
    listIn.insert(103,0)
    listIn.insert(104,0)    
    
    outDat = np.fliplr(np.array(listIn).reshape(15,7))
    outDat = np.flipud(outDat)
    return outDat

#fakeArray = list(np.arange(0,99))
#reshapeArray(fakeArray) #Testing to see if this works


fName = entries[1] #Load one file at a time
        
dat = pd.read_csv(fPath+fName,sep='\t', skiprows = 9, header = 0)

# Create Force signal
dat['forceTot'] = dat.iloc[:,100:198].sum(axis=1)
forceTot = dat['forceTot']
forceTot[forceTot<fThresh] = 0
plt.plot(forceTot[0:150])

#find the landings and offs of the FP as vectors
landings = findLandings(forceTot)
takeoffs = findTakeoffs(forceTot)

HSarray = []
HSdorsal = []
MSarray = []
MSdorsal = []
TOarray = []
TOdorsal = []

# loop through landings, extract rows defined above at landing
bufferLen = len(dat)
for landing in landings:
    hsindex = landing + hs
    mdindex = landing + ms
    toindex = landing + to
    if landing < 15:
        print(landing)
    elif (landing + 20 > bufferLen): 
        print(landing)
    else:
        HSarray.append(list(dat.iloc[hsindex,99:198]))
        HSdorsal.append(list(dat.iloc[hsindex,1:100]))
        MSarray.append(list(dat.iloc[mdindex,99:198]))
        MSdorsal.append(list(dat.iloc[mdindex,1:100]))
        TOarray.append(list(dat.iloc[toindex,99:198]))
        TOdorsal.append(list(dat.iloc[toindex,1:100]))


noRows = len(HSarray)
noCols = 99    
    
# Plantar pressure aggregation and mean
HSlist = list(np.mean(np.array(HSarray).reshape((noRows, noCols)), axis = 0))
MSlist = list(np.mean(np.array(MSarray).reshape((noRows, noCols)), axis = 0))
TOlist = list(np.mean(np.array(TOarray).reshape((noRows, noCols)), axis = 0))

# Dorsal pressure aggregation and mean
HSlistd = list(np.mean(np.array(HSdorsal).reshape((noRows, noCols)), axis = 0))
MSlistd = list(np.mean(np.array(MSdorsal).reshape((noRows, noCols)), axis = 0))
TOlistd = list(np.mean(np.array(TOdorsal).reshape((noRows, noCols)), axis = 0))

# Reshape to be correct for mapping
hsAvg = reshapeArray(HSlist)
msAvg = reshapeArray(MSlist)
toAvg = reshapeArray(TOlist)

# Reshape to be correct for mapping
hsAvgd = reshapeArray(HSlistd)
msAvgd = reshapeArray(MSlistd)
toAvgd = reshapeArray(TOlistd)


fig, ( ax1, ax2, ax3 ) = plt.subplots(1,3)
g1 = sns.heatmap(hsAvg, cmap="jet", ax = ax1, vmin = 0, vmax = 250)
g1.set(xticklabels=[])
g1.set_title('Plantar Initial Contact')
#mid stance
fig.suptitle('Average Pressures across the gait cycle')
g2 = sns.heatmap(msAvg, cmap="jet", ax = ax2, vmin = 0, vmax = 250)
g2.set(xticklabels=[])
g2.set_title('Plantar Mid Stance')
# Toe off
fig.suptitle('Average Pressures across the gait cycle')
g3 = sns.heatmap(toAvg, cmap="jet", ax = ax3, vmin = 0, vmax = 250)
g3.set(xticklabels=[])
g3.set_title('Plantar Toe Off')
fig.tight_layout()
### options for colors include winter, autumn, blue, and more ###
### need to find a good option for our purpose still 
