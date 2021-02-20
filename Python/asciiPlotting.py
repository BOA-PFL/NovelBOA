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


fName = entries[0] #Load one file at a time
        
dat = pd.read_csv(fPath+fName,sep='\t', skiprows = 9, header = 0)

# Create Force signal
dat['forceTot'] = dat.iloc[:,100:198].sum(axis=1)
forceTot = dat['forceTot']
forceTot[forceTot<fThresh] = 0
plt.plot(forceTot[0:150])

#find the landings and offs of the FP as vectors
landings = findLandings(forceTot)
takeoffs = findTakeoffs(forceTot)

tmpArray = list(dat.iloc[100,99:198])
newDat = reshapeArray(tmpArray)
plt.imshow(newDat, cmap='hot')


startVal = 100
for i in np.arange(1,6,1):
    tmpArray = list(dat.iloc[i+startVal,99:198])
    newDat = reshapeArray(tmpArray)
    plt.imshow(newDat, cmap='autumn')
    plt.show()
    plt.pause(1)
    plt.close()