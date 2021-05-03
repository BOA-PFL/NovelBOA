# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 10:28:46 2021
This script takes in two ascii files from Pedar and runs 2D SPM on them
They are assumed to be 201 columns and the standard export from Novel.
The output is plots and SPM at 3 phases of gait (dorsal and plantar) or 
Left and right. This script is based off the peaks in the data not landings!
You need to specify how many indices after heel strike you want to 
normalize the parts of the stride to. Hike/walk is usually around 5, 25, and 35
where run should be closer to 3, 12, and 18

@author: Daniel.Feeney
"""

import spm1d
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from scipy.signal import find_peaks



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


def findRightLength(array1, array2):
    # input the arrays 1 and 2 from the two conditions an find the shorter ##
    # Array to make them the same shape (required for later in spm)
    if len(array1) < len(array2):
        correctLength = len(array1)
    else:
        correctLength = len(array2)
    return(correctLength)
    
def calcStackedArray(inputArray, correctArrayLen):
    ### Used to input array data in list format and output matrix ###
    ### in a 3D format for SPM ###
    preallocArray = np.zeros((correctLengthHS,15,7))
    rowIndex = 0
    for row in inputArray:
        if rowIndex < correctArrayLen:
            preallocArray[rowIndex,:,:] = np.flipud(reshapeArray(row))
            rowIndex = rowIndex + 1
        else:
            print('skipped row')
    return preallocArray

#%matplotlib qt
# Read in files
# only read .asc files for this work
fPath = 'C:/Users/kate.harrison/Dropbox (Boa)/EndurancePerformance/SauconyEndorphinPro_Feb2021/PedarData/'
fileExt = r".asc"
entries = [fName for fName in os.listdir(fPath) if fName.endswith(fileExt)]

fName = entries[0] # Pick 2 files
fName2 = entries[1]
        
dat = pd.read_csv(fPath+fName,sep='\t', skiprows = 9, header = 0)
dat2 = pd.read_csv(fPath+fName2,sep='\t', skiprows = 9, header = 0)

### in order to plot same time around the peak force, need to specify HS, Mid Stance
### and toe off times here. For running, -4, 0, and +5 seem to work
hs = -5
ms = 0
to = 7
# Define constants and options

# Create Force signal
dat['forceTot'] = dat.iloc[:,100:198].sum(axis=1)
forceTot = dat['forceTot']
forceTot[forceTot<700] = 0

dat2['forceTot'] = dat2.iloc[:,100:198].sum(axis=1)
forceTot2 = dat2['forceTot']

#find the peaks and offs of the FP as vectors
peaks, _ = find_peaks(forceTot, height=800)



# Create arrays for heel strike from dataset 1
HSarray = []
HSdorsal = []
MSarray = []
MSdorsal = []
TOarray = []
TOdorsal = []

# loop through landings, extract rows defined above at landing
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
    
# Create arrays for heel strike from dataset 2
peaks2, _ = find_peaks(forceTot2, height=800)


HSarray2 = []
HSdorsal2 = []
MSarray2 = []
MSdorsal2 = []
TOarray2 = []
TOdorsal2 = []

# loop through landings, extract rows defined above at landing
bufferLen = len(dat2)
for peak in peaks2:
    hsindex = peak + hs
    mdindex = peak + ms
    toindex = peak + to
    if peak < 15:
        print(peak)
    elif (peak + 20 > bufferLen): 
        print(peak)
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
# Heel Strike Data
correctLengthHS = findRightLength(HSarray, HSarray2)
stackedHSData = calcStackedArray(HSarray, correctLengthHS)
stackedHSData2 = calcStackedArray(HSarray2, correctLengthHS)
stackedHSDorsal = calcStackedArray(HSdorsal, correctLengthHS)
stackedHSDorsal2 = calcStackedArray(HSdorsal2, correctLengthHS)
# MS Data
correctLengthMS = findRightLength(MSarray, MSarray2)
stackedMSData = calcStackedArray(MSarray, correctLengthMS)
stackedMSData2 = calcStackedArray(MSarray2, correctLengthMS)
stackedMSDorsal = calcStackedArray(MSdorsal, correctLengthMS)
stackedMSDorsal2 = calcStackedArray(MSdorsal2, correctLengthMS)
# TO Data
correctLengthTO = findRightLength(TOarray, TOarray2)
stackedTOData = calcStackedArray(TOarray, correctLengthTO)
stackedTOData2 = calcStackedArray(TOarray2, correctLengthTO)
stackedTODorsal = calcStackedArray(TOdorsal, correctLengthTO)
stackedTODorsal2 = calcStackedArray(TOdorsal2, correctLengthTO)

######################################################################## 
# Workhorse function to calculate SPM and make plots
def calcSPM2d(stackedDat1, stackedDat2, gaitPhase):
    ### Run all parts of SPM 2D and return plots ###
    mA  = stackedDat1.mean(axis=0)
    mB  = stackedDat2.mean(axis=0)
    m   = np.hstack( [mA, mB] )
    
    plt.figure()
    ax = plt.axes()
    ax.imshow(np.ma.masked_array(m, m==0), origin='lower', cmap='jet')
    ax.set_title(gaitPhase)
    plt.axvline(x = 6.5, color = 'r')
    cb = plt.colorbar(mappable=ax.images[0])
    cb.set_label('Average pressure (kPa)')
    plt.show()
    
    ### Statistical testing ###
    y       = np.array([yy.flatten() for yy in stackedDat1])
    J,Q     = y.shape
    y2       = np.array([yy.flatten() for yy in stackedDat2])
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
    ynz2    = y2[:,i]
    ynzA    = ynz  #Task A observations with non-zero variance nodes removed
    ynzB    = ynz2  #Task B observations with non-zero variance nodes removed
    
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
    ax.set_title(gaitPhase)
    plt.axvline(x=6.5)
    cb = plt.colorbar(mappable=ax.images[0])
    cb.set_label('t value')
    plt.show()
    
## Calculate SPM with functions below for 6 phases
    
calcSPM2d(stackedHSData, stackedHSData2, 'Heel Strike')
#calcSPM2d(stackedHSDorsal, stackedHSDorsal2, 'Heel Strike Dorsal')

calcSPM2d(stackedMSData, stackedMSData2, 'Mid Stance')
#calcSPM2d(stackedMSDorsal, stackedMSDorsal2, 'Mid Stance Dorsal')

calcSPM2d(stackedTOData, stackedTOData2, 'Toe Off')
#calcSPM2d(stackedTODorsal, stackedTODorsal2, 'Toe Off Dorsal')
