# -*- coding: utf-8 -*-
"""
Created on Wed May 12 08:18:38 2021

@author: Daniel.Feeney
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os


# Delineate the climb and downhill portion of trial
def delimitTrial(inputDF):
    # generic function to plot and start/end trial #
    fig, ax = plt.subplots()
    ax.plot(inputDF.Left, label = 'Left Force')
    fig.legend()
    pts = np.asarray(plt.ginput(2, timeout=-1))
    plt.close()
    outputDat = dat.iloc[int(np.floor(pts[0,0])) : int(np.floor(pts[1,0])),:]
    outputDat = outputDat.reset_index()
    return(outputDat)

def defThreshold(inputDF):
        # find threshold force
        fig, ax = plt.subplots()
        ax.plot(inputDF.Left, label = 'Right Foot Force')
        print('Select a point to represent 0 in the trial')
        pts = np.asarray(plt.ginput(1, timeout=-1))
        plt.close()
        fThresh = pts[0][1]
        return(fThresh)

def trimForce(inputDF, threshForce):
    forceTot = inputDF.Left
    forceTot[forceTot<threshForce] = 0
    forceTot = np.array(forceTot)
    return(forceTot)

def findLandings(force):
    lic = []
    for step in range(len(force)-1):
        if force[step] == 0 and force[step + 1] >= fThresh:
            lic.append(step)
    return lic

# preallocate matrix for force and fill in with force data
def forceMatrix(inputForce, landings, noSteps, stepLength):
    #input a force signal, return matrix with n rows (for each landing) by m col
    #for each point in stepLen
    preForce = np.zeros((noSteps,stepLength))
    
    for iterVar, landing in enumerate(landings):
        try:
            preForce[iterVar,] = inputForce[landing:landing+stepLength]
        except:
            print(landing)
            
    return preForce

# Read in file(s) and add names
fPath = 'C:\\Users\\Daniel.Feeney\\Dropbox (Boa)\\EndurancePerformance\\SalomonQuicklace_Aug2020\\'
fileExt = r".txt"
entries = [fName for fName in os.listdir(fPath) if fName.endswith(fileExt)]

fName = entries[0] #temporarily hard coding one file

dat = pd.read_csv(fPath+fName,sep='\s+', skiprows = 3, header = 0)
dat.columns = ['Time', 'LeftHeel', 'LeftMedial','LeftLateral','Left','Time2','RightLateral','RightMedial','RightHeel','Right']

print('Select start and end of climbing part of trial')
climbDat = delimitTrial(dat)
climbThresh = defThreshold(climbDat)
climbForce = trimForce(climbDat, climbThresh)

print('Select start and end of descending part of trial')
downDat = delimitTrial(dat)
downThresh = defThreshold(downDat)
downForce = trimForce(downDat, downThresh)

# delimit steps on left side
climbLandings = findLandings(climbForce)
downLandings = findLandings(downForce)



stepLen = 25
x = np.linspace(0,stepLen,stepLen)
stackedClimbF = forceMatrix(climbForce, climbLandings, 10, 25)
stackedDownF = forceMatrix(downForce, downLandings, 10, 25)


avgClimbF = np.mean(stackedClimbF, axis = 0)
sdClimbF = np.std(stackedClimbF, axis = 0)
avgDownF = np.mean(stackedDownF, axis = 0)
sdDownF = np.std(stackedDownF, axis = 0)

fig, (ax1, ax2) = plt.subplots(2)
ax1.plot(x,avgClimbF, 'k', color='#003D4C', label = 'climb force')
ax1.fill_between(x,avgClimbF-sdClimbF, avgClimbF+sdClimbF,
    alpha=0.5, edgecolor='#003D4C', facecolor='#003D4C')
ax1.set_ylim([0,2200])
ax1.set_title('Uphill Force')
ax1.set_ylabel('Force (N)')
ax2.plot(x,avgDownF, 'k', color='#00966C', label = 'climb force')
ax2.fill_between(x,avgDownF-sdDownF, avgDownF+sdDownF,
    alpha=0.5, edgecolor='#00966C', facecolor='#00966C')
ax2.set_ylim([0,2200])
ax2.set_title('Downhill Force')
ax2.set_ylabel('Force (N)')
ax2.set_xlabel('Time (cs)')
plt.tight_layout()