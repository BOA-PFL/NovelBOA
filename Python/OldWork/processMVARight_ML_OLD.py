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
fPath = 'C:/Users/kate.harrison/Dropbox (Boa)/AgilityPerformance/BOA_mechanisticStudy/Skater_MVA/'
fileExt = r".mva"
entries = [fName for fName in os.listdir(fPath) if fName.endswith(fileExt)]

# Define constants and options
fThresh = 0 #below this value will be set to 0.
#stepLen = 45 #Set value to look forward 
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
sdRHeel = []
meanRHeel = []
sdRLatFF = []
meanRLatFF = []
sdRMedFF = []
meanRMedFF = []
sdRLatMF = []
meanRLatMF = []
sdRMedMF = []
meanRMedMF = []
sdRToes = []
meanRToes = []

maxRHeel = []
maxRLatFF = []
maxRMedFF = []
maxRLatMF = []
maxRMedMF = []
maxRToes = []

CT = []

# sdRDorsalLatFF = []
# meanRDorsalLatFF = []
# sdRDorsalMedFF = []
# meanRDorsalMedFF = []
# sdRDorsalLatMF = []
# meanRDorsalLatMF = []
# sdRDorsalMedMF = []
# meanRDorsalMedMF = []

# maxRDorsalLatFF = []
# maxRDorsalMedFF = []
# maxRDorsalLatMF = []
# maxRDorsalMedMF = []

#HeelPMidStance = []
#HeelRateDecay = []
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
        ConditionTmp = fName.split(sep="_")[2]
        ConfigTmp = fName.split(sep="_")[1]
        
        dat = pd.read_csv(fPath+fName,sep='\t', skiprows = 16, header = 0)
        
        dat.columns = ['Time','RHeel_Force', 'RHeel_MaxP', 'RHeel_MeanP', 'RHeel_pct', 
                       'RMedMF_Force','RMedMF_MaxP', 'RMedMF_MeanP', 'RMedMF_pct', 
                       'RLatMF_Force', 'RLatMF_MaxP', 'RLatMF_MeanP','RLatMF_pct',
                       'RMedFF_Force','RMedFF_MaxP', 'RMedFF_MeanP', 'RMedFF_pct', 
                       'RLatFF_Force', 'RLatFF_MaxP', 'RLatFF_MeanP','RLatFF_pct',
                       'RToes_Force','RToes_MaxP', 'RToes_MeanP', 'RToes_Pct'       
            ]
        
        dat['Force'] = dat.RHeel_Force + dat.RMedMF_Force + dat.RLatMF_Force + dat.RMedFF_Force + dat.RLatFF_Force + dat.RToes_Force
          #filtering force to find landings/takeoffs    
          
        plt.plot(dat.Force)
        minClick = plt.ginput(1)
        minExtract = minClick[0]
        (start, fThresh) = minExtract
        start = round(start)
        endClick = plt.ginput(1)
        endExtract = endClick[0]
        (finish, y) = endExtract
        finish = round(finish)
        forceTot = dat.Force
        forceTot[forceTot<fThresh] = 0
        forceTot[0:start] = 0
        forceTot[finish:-1] = 0
        plt.clf()
        
    
       #find the landings and offs of the FP as vectors
        landings = findLandings(forceTot)
        takeoffs = findTakeoffs(forceTot)
        
        for landing in range(len(landings)):
            try:
               # Appending dorsal pressure data
              
                #HeelPMidStance.append(dat.HeelMaxP[landing+25])
                #HeelRateDecay.append(dat.HeelMaxP[landing+10] - dat.HeelMaxP[landing+18])
                sdRHeel.append(np.std(dat.RHeel_MeanP[landings(landing):takeoffs(landing)])) 
                meanRHeel.append(np.mean(dat.RHeel_MeanP[landings(landing):takeoffs(landing)])) 
                sdRLatMF.append(np.std(dat.RLatMF_MeanP[landings(landing):takeoffs(landing)])) 
                meanRLatMF.append(np.mean(dat.RLatMF_MeanP[landings(landing):takeoffs(landing)])) 
                sdRMedMF.append(np.std(dat.RMedMF_MeanP[landings(landing):takeoffs(landing)])) 
                meanRMedMF.append(np.mean(dat.RMedMF_MeanP[landings(landing):takeoffs(landing)])) 
                sdRLatFF.append(np.std(dat.RLatFF_MeanP[landings(landing):takeoffs(landing)])) 
                meanRLatFF.append(np.mean(dat.RLatFF_MeanP[landings(landing):takeoffs(landing)])) 
                sdRMedFF.append(np.std(dat.RMedFF_MeanP[landings(landing):takeoffs(landing)])) 
                meanRMedFF.append(np.mean(dat.RMedFF_MeanP[landings(landing):takeoffs(landing)])) 
                sdRToes.append(np.std(dat.RToes_MeanP[landings(landing):takeoffs(landing)])) 
                meanRToes.append(np.mean(dat.RToes_MeanP[landings(landing):takeoffs(landing)])) 
                maxRHeel.append(np.max(dat.RHeel_MaxP[landings(landing):takeoffs(landing)])) 
                maxRLatMF.append(np.max(dat.RLatMF_MaxP[landings(landing):takeoffs(landing)])) 
                maxRMedMF.append(np.max(dat.RMedMF_MaxP[landings(landing):takeoffs(landing)])) 
                maxRLatFF.append(np.max(dat.RLatFF_MaxP[landings(landing):takeoffs(landing)])) 
                maxRMedFF.append(np.max(dat.RMedFF_MaxP[landings(landing):takeoffs(landing)])) 
                maxRToes.append(np.max(dat.RToes_MaxP[landings(landing):takeoffs(landing)])) 
               
                trial.append(fName)
                # Appending Plantar pressure data
                # sdRDorsalMedFF.append(np.std(dat.RDorsalMedFF_MeanP[landing:landing+stepLen]))
                # meanRDorsalMedFF.append(np.mean(dat.RDorsalMedFF_MeanP[landing:landing+stepLen]))
                # sdRDorsalLatFF.append(np.std(dat.RDorsalLatFF_MeanP[landing:landing+stepLen]))
                # meanRDorsalLatFF.append(np.mean(dat.RDorsalLatFF_MeanP[landing:landing+stepLen]))
                # sdRDorsalMedMF.append(np.std(dat.RDorsalMedMF_MeanP[landing:landing+stepLen]))
                # meanRDorsalMedMF.append(np.mean(dat.RDorsalMedMF_MeanP[landing:landing+stepLen]))
                # sdRDorsalLatMF.append(np.std(dat.RDorsalLatMF_MeanP[landing:landing+stepLen]))
                # meanRDorsalLatMF.append(np.mean(dat.RDorsalLatMF_MeanP[landing:landing+stepLen]))
                # maxRDorsalMedFF.append(np.max(dat.RDorsalMedFF_MaxP[landing:landing+stepLen]))
                # maxRDorsalLatFF.append(np.max(dat.RDorsalLatFF_MaxP[landing:landing+stepLen]))
                # maxRDorsalMedMF.append(np.max(dat.RDorsalMedMF_MaxP[landing:landing+stepLen]))  
                # maxRDorsalLatMF.append(np.max(dat.RDorsalLatMF_MaxP[landing:landing+stepLen])) 
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
        
outcomes = pd.DataFrame({'Subject':list(Subject),'Condition':list(Condition), 'Config':list(Config),
                         'sdRHeel': list(sdRHeel),'meanRHeel':list(meanRHeel), 
                         'sdRLatMF':list(sdRLatMF), 'meanRLatMF':list(meanRLatMF),
                         'sdRMedMF':list(sdRMedMF),'meanRMedMF':list(meanRMedMF),
                         'sdRLatFF':list(sdRLatFF), 'meanRLatFF':list(meanRLatFF),
                         'sdRMedFF':list(sdRMedFF),'meanRMedFF':list(meanRMedFF), 
                         'sdRToes':list(sdRToes),'meanRToes':list(meanRToes), 
                         'maxRHeel':list(maxRHeel), 'maxRLatMF':list(maxRLatMF), 'maxRMedMF':list(maxRMedMF), 'maxRLatFF':list(maxRLatFF), 'maxRMedFF':list(maxRMedFF), 'maxRToes':list(maxRToes)
                         })

outcomes.to_csv(path_or_buf="C:/Users/kate.harrison/Dropbox (Boa)/AgilityPerformance/BOA_mechanisticStudy/Skater_MVA/CompiledPressureData.csv", sep =",")
    
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
