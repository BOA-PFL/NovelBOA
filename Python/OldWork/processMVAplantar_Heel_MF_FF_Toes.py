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
fPath = 'C:/Users/kate.harrison/Dropbox (Boa)/EndurancePerformance/NewBalanceRoadRacing_Jan2020/PressureData/'
fileExt = r".mva"
entries = [fName for fName in os.listdir(fPath) if fName.endswith(fileExt)]

# Define constants and options
fThresh = 30 #below this value will be set to 0.
stepLen = 45 #Set value to look forward 
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
sdFF = []
meanFF = []
sdMF = []
meanMF = []
sdMets = []
meanMets = []
maxMF = []
maxFF = []
maxMets = []

sdPlantMetP = []
meanPlantMetP = []
sdToeP = []
meanToeP = []
sdHeelP = []
meanHeelP = []
maxPlantMetP = []
maxToeP = []
maxHeelP = []

HeelPMidStance = []
HeelRateDecay = []
trial = []
Subject = []
Condition = []
Config = []

#first columns (FF, Mets, and MF) all relate to dorsal values. Once you get to PlantMetsForce it is plantar metatarsal force 
#and everything to the right of that column is plantar side. Each location (e.g. FF, MF, etc.) has force, max Pressure, Mean Pressure, and pct

for file in entries:
    try:
        
        fName = file #Load one file at a time
        
        subName = fName.split(sep = "_")[0]
        ConditionTmp = fName.split(sep="_")[1]
        ConfigTmp = fName.split(sep="_")[2]
        
        dat = pd.read_csv(fPath+fName,sep='\t', skiprows = 15, header = 0)
        
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
              
                HeelPMidStance.append(dat.HeelMaxP[landing+25])
                HeelRateDecay.append(dat.HeelMaxP[landing+10] - dat.HeelMaxP[landing+18])
                sdFF.append(np.std(dat.FFMeanP[landing:landing+stepLen])) 
                meanFF.append(np.mean(dat.FFMeanP[landing:landing+stepLen]))
                sdMets.append(np.std(dat.MetsMeanP[landing:landing+stepLen]))
                meanMets.append(np.mean(dat.MetsMeanP[landing:landing+stepLen]))
                maxFF.append(np.max(dat.FFMaxP[landing:landing+stepLen]))
                maxMF.append(np.max(dat.MFMaxP[landing:landing+stepLen]))
                maxMets.append(np.max(dat.MetsMaxP[landing:landing+stepLen]))
                sdMF.append(np.std(dat.MFMeanP[landing:landing+stepLen])) 
                meanMF.append(np.mean(dat.MFMeanP[landing:landing+stepLen]))
                trial.append(fName)
                # Appending Plantar pressure data
                sdPlantMetP.append(np.std(dat.PlantMetsMeanP[landing:landing+stepLen]))
                meanPlantMetP.append(np.mean(dat.PlantMetsMeanP[landing:landing+stepLen]))
                sdToeP.append(np.std(dat.ToesMeanP[landing:landing+stepLen]))
                meanToeP.append(np.mean(dat.ToesMeanP[landing:landing+stepLen]))
                sdHeelP.append(np.std(dat.HeelMeanP[landing:landing+stepLen]))
                meanHeelP.append(np.mean(dat.HeelMeanP[landing:landing+stepLen]))
                maxPlantMetP.append(np.max(dat.PlantMetsMaxP[landing:landing+stepLen]))
                maxToeP.append(np.max(dat.ToesMaxP[landing:landing+stepLen]))
                maxHeelP.append(np.max(dat.HeelMaxP[landing:landing+stepLen]))  
                Subject.append(subName)
                Condition.append(ConditionTmp)
                Config.append(ConfigTmp)


            except:
                print(landing)
    except:
        print(file)    
    
#the outcomes dataframe uses the same logic (first set are dorsal pressures until maxPlantMetP, 
#then everything is plantar unless noted otherwise), then midstanceHeelP is the heel pressure at 
#mid stance and the heel decay from peak (which is a bit finnicky at the moment). 
        
outcomes = pd.DataFrame({'Subject':list(Subject),'Condition':list(Condition), 'Config':list(Config),'sdMets': list(sdMets),
                         'meanMets':list(meanMets), 'sdFF':list(sdFF), 'meanFF':list(meanFF),'sdMF':list(sdMF),
                         'meanMF':list(meanMF), 'maxFF':list(maxFF), 'maxMF':list(maxMF), 'maxMets':list(maxMets), 
                         'sdPlantMetP':list(sdPlantMetP), 'meanPlantMetP':list(meanPlantMetP),
                         'sdToeP':list(sdToeP), 'meanToeP':list(meanToeP),
                         'sdHeelP':list(sdHeelP), 'meanHeelP':list(meanHeelP),
                         'maxPlantMetP':list(maxPlantMetP),
                         'maxToeP':list(maxToeP),'maxHeelP':list(maxHeelP), 'MidStanceHeelP':list(HeelPMidStance),
                         'HeelDecay':list(HeelRateDecay)})

outcomes.to_csv(path_or_buf="C:/Users/kate.harrison/Dropbox (Boa)/EndurancePerformance/NewBalanceRoadRacing_Jan2020/PressureData/CompiledPressureData.csv", sep =",")
    
## plotting
#landingToPlot = 3
#fig, ax1 = plt.subplots()
#
#color = 'tab:red'
#ax1.set_xlabel('time')
#ax1.set_ylabel('TotalForce(N)', color=color)
#ax1.plot(dat.Force[landings[landingToPlot]:landings[landingToPlot]+stepLen], color=color)
#ax1.tick_params(axis='y', labelcolor=color)
#
#ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
#
#ax2.set_ylabel('Max Pressures')  # we already handled the x-label with ax1
#ax2.plot(dat.FFMaxP[landings[landingToPlot]:landings[landingToPlot]+stepLen])
#ax2.plot(dat.MetsMaxP[landings[landingToPlot]:landings[landingToPlot]+stepLen])
#ax2.plot(dat.MFMaxP[landings[landingToPlot]:landings[landingToPlot]+stepLen])
#plt.legend()
#
#fig.tight_layout()  # otherwise the right y-label is slightly clipped
#plt.show()
#
## mean pressure
#fig, ax1 = plt.subplots()
#
#color = 'tab:red'
#ax1.set_xlabel('time')
#ax1.set_ylabel('TotalForce(N)', color=color)
#ax1.plot(dat.Force[landings[landingToPlot]:landings[landingToPlot]+stepLen], color=color)
#ax1.tick_params(axis='y', labelcolor=color)
#
#ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
#
#ax2.set_ylabel('Max Pressures')  # we already handled the x-label with ax1
#ax2.plot(dat.FFMeanP[landings[landingToPlot]:landings[landingToPlot]+stepLen])
#ax2.plot(dat.MetsMeanP[landings[landingToPlot]:landings[landingToPlot]+stepLen])
#ax2.plot(dat.MFMeanP[landings[landingToPlot]:landings[landingToPlot]+stepLen])
#plt.legend()
#
#fig.tight_layout()  # otherwise the right y-label is slightly clipped
#plt.show()
    
# mean pressure
#fig, ax1 = plt.subplots()
#
#color = 'tab:red'
#ax1.set_xlabel('time')
#ax1.set_ylabel('TotalForce(N)', color=color)
#ax1.plot(dat.Force[landings[landingToPlot]:landings[landingToPlot]+stepLen], color=color)
#ax1.tick_params(axis='y', labelcolor=color)
#
#ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
#
#ax2.set_ylabel('Max Pressures')  # we already handled the x-label with ax1
##ax2.plot(dat.FFMeanP[landings[landingToPlot]:landings[landingToPlot]+stepLen])
#ax2.plot(dat.PlantMetsMeanP[landings[landingToPlot]:landings[landingToPlot]+stepLen])
#ax2.plot(dat.HeelMaxP[landings[landingToPlot]:landings[landingToPlot]+stepLen])
#plt.legend()
#
#fig.tight_layout()  # otherwise the right y-label is slightly clipped
#plt.show()
