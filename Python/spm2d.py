# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 10:28:46 2021

@author: Daniel.Feeney
"""

import spm1d
import numpy as np
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns

## example
#fname = 'C:/Users/Daniel.Feeney/OneDrive - Boa Technology Inc/Documents/data2d.npy'
#Y     = np.load(fname)

#print(Y.shape)

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
to = 35

# Define constants and options
fThresh = 215 #below this value will be set to 0.
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
#fakeReshaped = reshapeArray(fakeArray) #Testing to see if this works


fName = entries[3] #Load one file at a time
fName2 = entries[4]
        
dat = pd.read_csv(fPath+fName,sep='\t', skiprows = 9, header = 0)
dat2 = pd.read_csv(fPath+fName2,sep='\t', skiprows = 9, header = 0)

# Create Force signal
dat['forceTot'] = dat.iloc[:,100:198].sum(axis=1)
forceTot = dat['forceTot']
forceTot[forceTot<fThresh] = 0

dat2['forceTot'] = dat2.iloc[:,100:198].sum(axis=1)
forceTot2 = dat2['forceTot']
forceTot2[forceTot2<fThresh] = 0

#find the landings and offs of the FP as vectors
landings = findLandings(forceTot)
takeoffs = findTakeoffs(forceTot)


# Create arrays for heel strike from dataset 1
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
    
# Create arrays for heel strike from dataset 2
landings2 = findLandings(forceTot2)

HSarray2 = []
HSdorsal2 = []
MSarray2 = []
MSdorsal2 = []
TOarray2 = []
TOdorsal2 = []

# loop through landings, extract rows defined above at landing
bufferLen = len(dat2)
for landing in landings2:
    hsindex = landing + hs
    mdindex = landing + ms
    toindex = landing + to
    if landing < 15:
        print(landing)
    elif (landing + 20 > bufferLen): 
        print(landing)
    else:
        HSarray2.append(list(dat2.iloc[hsindex,99:198]))
        HSdorsal2.append(list(dat2.iloc[hsindex,1:100]))
        MSarray2.append(list(dat2.iloc[mdindex,99:198]))
        MSdorsal2.append(list(dat2.iloc[mdindex,1:100]))
        TOarray2.append(list(dat2.iloc[toindex,99:198]))
        TOdorsal2.append(list(dat2.iloc[toindex,1:100]))

### need to reshape these to no steps x width x height matrices ###
### Formatting data into spm requirements

# These need to be the same dimensions #
if len(HSarray) < len(HSarray2):
    correctLengthHS = len(HSarray)
else:
    correctLengthHS = len(HSarray2)

stackedHSData = np.zeros((correctLengthHS,15,7))
rowIndex = 0
for row in HSarray:
    if rowIndex < correctLengthHS:
        tmpArray = row
        stackedHSData[rowIndex,:,:] = np.flipud(reshapeArray(row))
        rowIndex = rowIndex + 1
    else:
        print('skipped row')
    
# Dataset 2
stackedHSData2 = np.zeros((correctLengthHS,15,7))
rowIndex = 0
for row in HSarray2:
    if rowIndex < correctLengthHS:
        tmpArray = row
        stackedHSData2[rowIndex,:,:] = np.flipud(reshapeArray(row))
        rowIndex = rowIndex + 1
    else:
        print('skipped row')
    
### plotting average value
### plotting subsets of data

## compute mean of two subsets of data
mA  = stackedHSData.mean(axis=0)
mB  = stackedHSData2.mean(axis=0)
m   = np.hstack( [mA, mB] )

plt.figure()
ax = plt.axes()
ax.imshow(np.ma.masked_array(m, m==0), origin='lower', cmap='jet')
ax.set_title('Mean pressure distributions')
plt.axvline(x = 6.5, color = 'r')
cb = plt.colorbar(mappable=ax.images[0])
cb.set_label('Average pressure (kPa)')
plt.show()

### Statistical testing ###
y       = np.array([yy.flatten() for yy in stackedHSData])
J,Q     = y.shape
y2       = np.array([yy.flatten() for yy in stackedHSData2])
J2,Q2     = y2.shape
#separate the observations:
yA      = y #flattend dataset 1
yB      = y2  #flattend observations for Task B

#find zero-variance nodes:
eps     = np.finfo(float).eps     #smallest float number
iA      = yA.std(axis=0) > eps    #Task A indices where variance > 0
iB      = yB.std(axis=0) > eps    #Task B indices where variance > 0
i       = np.logical_and(iA, iB)  #indices where variance > 0 for both tasks

ynz     = y[:,i]    #all observations with non-zero variance nodes removed
ynzA    = ynz[:10]  #Task A observations with non-zero variance nodes removed
ynzB    = ynz[10:]  #Task B observations with non-zero variance nodes removed

# Run the tests
snpm    = spm1d.stats.nonparam.ttest2(ynzA, ynzB)
snpmi   = snpm.inference(0.05, two_tailed=True, iterations=10000)

# Reshaping for visualization
znz     = snpmi.z      #flattened test statistic  (i.e., t value) over only non-zero-variance nodes
zstar   = snpmi.zstar  #critical test statistic  (i.e., critical t value)

z       = np.zeros(Q)  #flattened test statistic including all nodes
z[i]    = znz          #add test statistic values for the non-zero-variance nodes

Z       = np.reshape(z, stackedHSData.shape[1:])   #2D test statistic image
Zi      = Z.copy()                     #2D inference image (temporary)
Zi[np.abs(Z)<zstar] = 0                #thresholded test statistic image

ZZ      = np.hstack( [Z, Zi] )         #both raw test statistic image and inference image


#plot:
plt.figure()
ax = plt.axes()
ax.imshow(np.ma.masked_array(ZZ, ZZ==0), origin='lower', cmap='jet', vmin=-15, vmax=15)
ax.set_title('SPM results at Heel Strike')
plt.axvline(x=6.5)
cb = plt.colorbar(mappable=ax.images[0])
cb.set_label('t value')
plt.show()