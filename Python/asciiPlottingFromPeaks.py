# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 10:55:52 2020
# use this to average ascii files from Novel around the peaks rather than 
# from heel strike as in asciiPlottingOneSide and asciiPlotting
# This relies on the peaks being over a threshold defined in the script. 
# Called fThresh. 800 works
@author: Daniel.Feeney
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns
from scipy.signal import find_peaks


fPath = 'C:/Users/Daniel.Feeney/Dropbox (Boa)/EndurancePerformance/Altra_MontBlanc_Jan2021/PedarPressures/'
fPath = 'C:\\Users\\daniel.feeney\\Boa Technology Inc\\PFL - General\\Cycling2021\\EH_CyclingPilot_2021\\Pressures\\'
fileExt = r".asc"
entries = [fName for fName in os.listdir(fPath) if fName.endswith(fileExt)]

### in order to plot same time after landing, need to specify HS, Mid Stance
### and toe off times here. For running, 4, 8, and 12 at 50 Hz suggested. For walking
### 5, 20, and 30 are a good start at 50 Hz
hs = -10
ms = 0
to = 5

#fThresh = 800
# Define constants and options

# list of functions 
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

fName = entries[5] #Load one file at a time
        
dat = pd.read_csv(fPath+fName,sep='\t', skiprows = 16, header = 0)

# define trial and force threshold
# Create Force signal
dat['forceTot'] = dat.iloc[:,100:198].sum(axis=1)
fig, ax = plt.subplots()
ax.plot(dat.forceTot, label = 'Right Total Force')
fig.legend()
print('Select start and end of analysis trial')
pts = np.asarray(plt.ginput(2, timeout=-1))
plt.close()
# downselect the region of the dataframe you'd like
dat = dat.iloc[int(np.floor(pts[0,0])) : int(np.floor(pts[1,0])),:]
dat = dat.reset_index()

# find threshold force
fig, ax = plt.subplots()
ax.plot(dat.forceTot, label = 'Right Foot Force')
print('Select a point to threshold peak force (should be above all secondary peaks')
pts = np.asarray(plt.ginput(1, timeout=-1))
plt.close()
fThresh = pts[0][1]
# make np array for below operation
forceTot = dat['forceTot']
#find peaks
#find the peaks and offs of the FP as vectors
peaks, _ = find_peaks(forceTot, height=fThresh)

HSarray = []
HSdorsal = []
MSarray = []
MSdorsal = []
TOarray = []
TOdorsal = []

# loop through peaks, extract rows defined above at landing
bufferLen = len(dat)
for peak in peaks:
    hsindex = peak + hs
    mdindex = peak + ms
    toindex = peak + to
    if peak < 15:
        print(peak)
    elif (peak + 20 > bufferLen): 
        print(peak)
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

maxP = 150
fig, ( ax1, ax2, ax3 ) = plt.subplots(1,3)
g1 = sns.heatmap(hsAvg, cmap="jet", ax = ax1, vmin = 0, vmax = maxP)
g1.set(xticklabels=[])
g1.set_title('Plantar Initial Contact')
#mid stance
fig.suptitle('Average Pressures across the gait cycle')
g2 = sns.heatmap(msAvg, cmap="jet", ax = ax2, vmin = 0, vmax = maxP)
g2.set(xticklabels=[])
g2.set_title('Plantar Mid Stance')
# Toe off
fig.suptitle('Average Pressures across the gait cycle')
g3 = sns.heatmap(toAvg, cmap="jet", ax = ax3, vmin = 0, vmax = maxP)
g3.set(xticklabels=[])
g3.set_title('Plantar Toe Off')
fig.tight_layout()
### options for colors include winter, autumn, blue, and more ###
### need to find a good option for our purpose still 
