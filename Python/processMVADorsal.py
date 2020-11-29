# -*- coding: utf-8 -*-
"""
Created on Sun Nov 29 12:10:19 2020

@author: Daniel.Feeney
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# Read in files
# only read .asc files for this work
fPath = 'C:/Users/Daniel.Feeney/Dropbox (Boa)/EnduranceProtocolWork/PressureASCII/'
fileExt = r".mva"
entries = [fName for fName in os.listdir(fPath) if fName.endswith(fileExt)]

# Define constants and options
fThresh = 100 #below this value will be set to 0.
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

### File Structure: Dorsal Forefoot, Metatarsals (.1), Midfoot (.2), plantar mets  (.3), plantar toes (.4), plantar heel  (.5)
cvMF = []
cvFF = []
cvMets = []
maxMF = []
maxFF = []
maxMets = []

maxPlantMetP = []
maxToeP = []
maxHeelP = []

trial = []
for file in entries:
    try:
        
        fName = file #Load one file at a time
        
        dat = pd.read_csv(fPath+fName,sep='\t', skiprows = 14, header = 0)
        
        dat.columns = ['Time','FFForce', 'FFMaxP', 'FFMeanP', 'FFpct', 'MetsForce', 'MetsMaxP', 'MetsMeanP','Metspct', 
            'MFForce','MFMaxP', 'MFMeanP', 'MFpct', 'PlantMetsForce','PlantMetsMaxP', 'PlantMetsMeanP', 'PlantMetsPct',
            'ToesForce','ToesMaxP','ToesMeanP','ToesPct','HeelForce', 'HeelMaxP', 'HeelMeanP','HeelPct']
        dat['Force'] = dat.HeelForce + dat.ToesForce + dat.PlantMetsForce
         # filtering force to find landings/takeoffs        
        forceTot = dat.Force
        forceTot[forceTot<fThresh] = 0
        
       #find the landings and offs of the FP as vectors
        landings = findLandings(forceTot)
        takeoffs = findTakeoffs(forceTot)
        
        for landing in landings:
            try:
               # Appending dorsal pressure data

                cvFF.append(np.std(dat.FFMeanP[landing:landing+45]) / np.mean(dat.FFMeanP[landing:landing+45]))
                cvMets.append(np.std(dat.MetsMeanP[landing:landing+45]) / np.mean(dat.MetsMeanP[landing:landing+45]))
                maxFF.append(np.max(dat.MFMaxP[landing:landing+45]))
                maxMF.append(np.max(dat.MFMaxP[landing:landing+45]))
                maxMets.append(np.max(dat.MetsMaxP[landing:landing+45]))
                cvMF.append(np.std(dat.MFMeanP[landing:landing+45]) / np.mean(dat.MFMeanP[landing:landing+45]))
                trial.append(fName)
                # Appending Plantar pressure data
                maxPlantMetP.append(np.max(dat.PlantMetsMaxP[landing:landing+45]))
                maxToeP.append(np.max(dat.ToesMaxP[landing:landing+45]))
                maxHeelP.append(np.max(dat.HeelMaxP[landing:landing+45]))  

            except:
                print(landing)
    except:
        print(file)    
        
       
outcomes = pd.DataFrame({'Trial':list(trial), 'cvMF': list(cvMF), 'cvFF':list(cvFF), 'cvMets':list(cvMets),
                         'maxFF':list(maxFF), 'maxMF':list(maxMF), 'maxMets':list(maxMets), 'maxPlantMetP':list(maxPlantMetP),
                         'maxToeP':list(maxToeP),'maxHeelP':list(maxHeelP)})
#
#fig, ax1 = plt.subplots()
#
#color = 'tab:red'
#ax1.set_xlabel('time')
#ax1.set_ylabel('TotalForce(N)', color=color)
#ax1.plot(dat.Force[40:85], color=color)
#ax1.tick_params(axis='y', labelcolor=color)
#
#ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
#
#ax2.set_ylabel('Max Pressures')  # we already handled the x-label with ax1
#ax2.plot(dat.FFMaxP[40:85])
#ax2.plot(dat.MetsForce[40:85])
#ax2.plot(dat.MFMaxP[40:85])
#plt.legend()
#
#fig.tight_layout()  # otherwise the right y-label is slightly clipped
#plt.show()
#
