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


fPath = 'C:/Users/kate.harrison/Dropbox (Boa)/AgilityPerformance/BOA_mechanisticStudy/CMJ_ASCII/Max/'
fileExt = r".asc"
entries = [fName for fName in os.listdir(fPath) if fName.endswith(fileExt)]

### in order to plot same time after landing, need to specify HS, Mid Stance
### and toe off times here. For running, 4, 8, and 12 at 50 Hz suggested. For walking
### 5, 20, and 30 are a good start at 50 Hz
hs = -12
ms = 0
to = 12

fThresh = 800
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

### Need to create a for loop here to include all files ###

HSarray = []
HSdorsal = []
MSarray = []
MSdorsal = []
TOarray = []
TOdorsal = []

for fName in entries:
    
    #fName = entries[0] #Load one file at a time. NEED TO MAKE A LOOP
        
    dat = pd.read_csv(fPath+fName,sep='\t', skiprows = 9, header = 0)

    # Create Force signal. Sum force of 99 sensors under the foot.
    dat['forceTot'] = dat.iloc[:,100:198].sum(axis=1)
    forceTot = dat['forceTot']

    #find the peaks and offs of the FP as vectors
    peaks, _ = find_peaks(forceTot, height=fThresh)

    plt.plot(dat.forceTot)

    print('Click a zero point at the start time of usable data on the plot')
    minClick = plt.ginput(1)
    print('Click end time point of usable data on the plot')
    endClick = plt.ginput(1)
    plt.close()
    minExtract = minClick[0]
    endExtract = endClick[0]
    (begin, fThresh) = minExtract
    (finish, y) = endExtract
    begin = round(begin)
    finish = round(finish)

    dat.forceTot[dat.forceTot<fThresh] = 0

    landings = findLandings(forceTot)
    landings[:] = [x for x in landings if x > begin] # remove landings before the start point specified on graph
    takeoffs = findTakeoffs(forceTot)
    takeoffs[:] = [x for x in takeoffs if x > landings[0]] # remove takeoffs before first landing
    landings[:] = [x for x in landings if x < finish] # remove landings after end of usable data
     
    

    # loop through peaks, extract rows defined above at landing
    bufferLen = len(dat)
    
    for step in range(len(landings)):
        
        hsindex = landings[step] + round((takeoffs[step] - landings[step])*0.1)
        toindex = landings[step] + round((takeoffs[step] - landings[step])*0.9)
        mdindex = round(hsindex + (toindex - hsindex)/2)
        HSarray.append(list(dat.iloc[hsindex,99:198]))
        HSdorsal.append(list(dat.iloc[hsindex,1:100]))
        MSarray.append(list(dat.iloc[mdindex,99:198]))
        MSdorsal.append(list(dat.iloc[mdindex,1:100]))
        TOarray.append(list(dat.iloc[toindex,99:198]))
        TOdorsal.append(list(dat.iloc[toindex,1:100]))
    
    HSarray.append(list(dat.iloc[hsindex,99:198]))


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
g1 = sns.heatmap(hsAvg, cmap="jet", ax = ax1, vmin = 0, vmax = 200)
g1.set(xticklabels=[])
g1.set_title('Plantar Initial Contact')
#mid stance
fig.suptitle('Average Pressures across the gait cycle')
g2 = sns.heatmap(msAvg, cmap="jet", ax = ax2, vmin = 0, vmax = 200)
g2.set(xticklabels=[])
g2.set_title('Plantar Mid Stance')
# Toe off
fig.suptitle('Average Pressures across the gait cycle')
g3 = sns.heatmap(toAvg, cmap="jet", ax = ax3, vmin = 0, vmax = 200)
g3.set(xticklabels=[])
g3.set_title('Plantar Toe Off')
fig.tight_layout()
### options for colors include winter, autumn, blue, and more ###
### need to find a good option for our purpose still 
