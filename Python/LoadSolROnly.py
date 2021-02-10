# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 13:14:07 2021

@author: Daniel.Feeney
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# Define constants and options
fThresh = 85; #below this value will be set to 0.
minStepLen = 10; #minimal step length
writeData = 0; #will write to spreadsheet if 1 entered
desiredStepLength = 25; #length to look forward after initial contact
apple = 0; #1 for apple 0 otherwise

# Read in file and add names
fPath = 'C:/Users/Daniel.Feeney/Dropbox (Boa)/EndurancePerformance/Altra_MontBlanc_Jan2021/'
entries = os.listdir(fPath)
fName = entries[3] #temporarily hard coding one file

dat = pd.read_csv(fPath+fName,sep='\t', skiprows = 3, header = 0)
dat = dat.drop(dat.columns[[5,6,7,8,9]], axis=1)  
dat.columns = ['Time', 'RightLateral','RightMedial','RightHeel','RightTotal']

subName = fName.split(sep = "_")[0]
Test = fName.split(sep = "_")[1]
Config = fName.split(sep = "_")[2]
Condition = fName.split(sep = "_")[3].split(sep=".")[0]

#### Do the same thing for the right side #### 
##### Filter force below threshold to 0 #####
RForce = dat.RightTotal
RForce[RForce<fThresh] = 0

# delimit steps on left side
ric = []
count = 1;
for step in range(len(RForce)-1):
    if RForce[step] == 0 and RForce[step + 1] >= fThresh:
        ric.append(step)
        count = count + 1
#left toe off
rto = []
count = 1;
for step in range(len(RForce)-1):
    if RForce[step] >= fThresh and RForce[step + 1] == 0:
        rto.append(step + 1)
        count = count + 1

# Remove 0s from lists
#ric = [s for s in ric if s != 0]
#rto = [s for s in rto if s != 0]
#%%
# trim first contact/toe off if not a full step
if ric[0] > rto[0]:
    rto = np.delete(rto, [0])

    
if ric[-1] > rto[-1]:
    ric = ric[:-1]

### Calculate step lengths, remove steps that were too short and first three and last 2 steps
RStepLens = np.array(rto) - np.array(ric)
# remove first 3 and last two steps #
RStepLens = np.delete(RStepLens, [0,1,2])   #First 3
RStepLens = RStepLens[:-2]                  #Last 2
ric = np.delete(ric, [0,1,2])
ric = ric[:-2]
rto = np.delete(rto, [0,1,2])
rto = rto[:-2]

# Removing false steps
falseSteps = np.where(RStepLens < minStepLen) #Find indices to remove
# Remove all occurrences of elements below below minStepLen
RStepLens = RStepLens[ RStepLens > minStepLen ]
#Trim IC and TO for false steps
ric = np.delete(ric, falseSteps)
rto = np.delete(rto, falseSteps)

#%%
###### extract relevent features ###### 
RTot = []
RHeel = []
RLat = []
RMed = []
    
for landing in ric:
    RTot.append(RForce[landing:landing+desiredStepLength])
    RHeel.append(dat.RightHeel[landing:landing+desiredStepLength])
    RLat.append(dat.RightLateral[landing:landing+desiredStepLength])
    RMed.append(dat.RightMedial[landing:landing+desiredStepLength])
noSteps = len(ric)

####
RightMat = np.reshape(RTot, (noSteps,desiredStepLength))
RHeelMat = np.reshape(RHeel, (noSteps, desiredStepLength))
RLatMat = np.reshape(RLat, (noSteps, desiredStepLength))
RMedMat = np.reshape(RMed, (noSteps,desiredStepLength))

## plotting left side
fig, ax = plt.subplots(4)
for i in range(noSteps):
    ax[0].plot(RightMat[i,:])
    ax[1].plot(RHeelMat[i,:])
    ax[2].plot(RLatMat[i,:])
    ax[3].plot(RMedMat[i,:])
    
ax[0].set_title('Total Right Force')
ax[1].set_title('Right Heel Force')
ax[2].set_title('Right Lateral Force')
ax[3].set_title('Right Medial Force')

#%%

# Right side
MaxR = []
totImpulseR = []
pkHeelR = []
heelImpulseR = []
pkLatR = []
latImpulseR = []
pkMedR = []
medImpulseR = []
stanceTimeR = []
rateTotR = []  
nameR = []
configR = []
tmpR = []
 
for step in range(noSteps):
    #left
    MaxR.append(max(RightMat[step,:]))
    totImpulseR.append(sum(RightMat[step,:]))
    pkHeelR.append(max(RHeelMat[step,:]))
    heelImpulseR.append(sum(RHeelMat[step,:]))
    pkLatR.append(max(RLatMat[step,:]))
    latImpulseR.append(sum(RLatMat[step,:]))
    pkMedR.append(max(RMedMat[step,:]))
    medImpulseR.append(sum(RMedMat[step,:]))
    tmpF = RightMat[step,:]
    tmpMaxLoc = np.argmax(tmpF)
    rateTotR.append(max(tmpF) / (tmpMaxLoc/100))
    stanceTimeR.append(len(tmpF[tmpF>0]))
    nameR.append(subName)
    configR.append(Config)
    tmpR.append('R')

rightDat = pd.DataFrame({'Sub':list(nameR), 'Config': list(configR), 'Side': list(tmpR),'StanceTime': list(stanceTimeR),
              'VLR':list(rateTotR),'MaxF': list(MaxR),'pkHeel': list(pkHeelR), 'HeelImpulse': list(heelImpulseR),
              'PkLat': list(pkLatR), 'LatImp': list(latImpulseR), 'PkMed': list(pkMedR),
              'MedImp': list(medImpulseR)})

    


