# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 10:55:52 2020

@author: Daniel.Feeney
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

#%matplotlib qt
# Read in files
# only read .asc files for this work
fPath = 'C:/Users/Daniel.Feeney/Dropbox (Boa)/EnduranceProtocolWork/PressureASCII/'
fileExt = r".asc"
entries = [fName for fName in os.listdir(fPath) if fName.endswith(fileExt)]

### in order to plot same time after landing, need to specify HS, Mid Stance
### and toe off times here. For running, 5, 12, and 18 work. For walking
### 5, 20, and 30 are a good start
hs = 5
ms = 20
to = 30

# Define constants and options
fThresh = 135 #below this value will be set to 0.
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

fakeArray = list(np.arange(0,99))
reshapeArray(fakeArray) #Testing to see if this works


fName = entries[3] #Load one file at a time
        
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
MSarray = []
TOarray = []
# loop through landings, extract row at 5, 12, and 18 after landing 
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
        MSarray.append(list(dat.iloc[mdindex,99:198]))
        TOarray.append(list(dat.iloc[toindex,99:198]))

noRows = len(HSarray)
noCols = 99    
    
HSlist = list(np.mean(np.array(HSarray).reshape((noRows, noCols)), axis = 0))
MSlist = list(np.mean(np.array(MSarray).reshape((noRows, noCols)), axis = 0))
TOlist = list(np.mean(np.array(TOarray).reshape((noRows, noCols)), axis = 0))


hsAvg = reshapeArray(HSlist)
plt.imshow(hsAvg, cmap='Blues')

msAvg = reshapeArray(MSlist)
plt.imshow(msAvg, cmap='Blues')

toAvg = reshapeArray(TOlist)
plt.imshow(toAvg, cmap='Blues')
### options for colors include winter, autumn, blue, and more ###
### need to find a good option for our purpose still 

startVal = 100
for i in np.arange(1,6,1):
    tmpArray = list(dat.iloc[i+startVal,99:198])
    newDat = reshapeArray(tmpArray)
    plt.imshow(newDat, cmap='Blues')
    plt.show()
    plt.pause(1)
    plt.close()