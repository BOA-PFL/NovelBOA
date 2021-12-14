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
fPath = 'C:/Users/kate.harrison/Dropbox (Boa)/EndurancePerformance/PFL_ISB2021_HeelHoldStudy/MVA Files/Run/'
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

#cvHeel = [] 
#meanHeel = []
icHeel = []
msHeel = []
toHeel = []

#cvLatFF = []
#meanLatFF = []
icLatFF = []
msLatFF = []
toLatFF = []

#cvMedFF = []
#meanMedFF = []
icMedFF = []
msMedFF = []
toMedFF = []

#cvToes = []
#meanToes = []
icToes = []
msToes = []
toToes = []
        
# cvDorsMedFF = []
# meanDorsMedFF = []
# cvDorsLatFF = []
# meanDorsLatFF = []
# cvDorsMedMF = []
# meanDorsMedMF = []
# cvDorsLatMF = []
# meanDorsLatMF = []

trial = []
Subject = []
Condition = []
Config = []
contactTime = []
#first columns (FF, Mets, and MF) all relate to dorsal values. Once you get to PlantMetsForce it is plantar metatarsal force 
#and everything to the right of that column is plantar side. Each location (e.g. FF, MF, etc.) has force, max Pressure, Mean Pressure, and pct

for fName in entries:
    try:
        
        
        dat = pd.read_csv(fPath+fName,sep='\t', skiprows = 18, header = 0)
              
        subName = fName.split(sep = "_")[0]
        ConditionTmp = fName.split(sep="_")[2]
        ConfigTmp = fName.split(sep="_")[1]
       
        
        dat.columns = ['Time','RHeel_Force', 'RHeel_MaxP', 'RHeel_MeanP', 'RHeel_pct', 
                       'RMedFF_Force','RMedFF_MaxP', 'RMedFF_MeanP', 'RMedFF_pct', 
                       'RLatFF_Force', 'RLatFF_MaxP', 'RLatFF_MeanP','RLatFF_pct',
                       'RToes_Force','RToes_MaxP', 'RToes_MeanP', 'RToes_Pct',
                       'DorsMedFF_Force', 'DorsMedFF_MaxP', 'DorsMedFF_MeanP', 'DorsMedFF_Pct',
                       'DorsLatFF_Force', 'DorsLatFF_MaxP', 'DorsLatFF_MeanP', 'DorsLatFF_Pct',
                       'DorsMedMF_Force', 'DorsMedMF_MaxP', 'DorsMedMF_MeanP', 'DorsMedMF_Pct',
                       'DorsLatMF_Force', 'DorsLatMF_MaxP', 'DorsLatMF_MeanP', 'DorsLatMF_Pct'
                       ]
        
        dat['Force'] = dat.RHeel_Force + dat.RLatFF_Force + dat.RMedFF_Force + dat.RToes_Force
        #filtering force to find landings/takeoffs    
          
        plt.plot(dat.Force)
        print('Click a zero point/starting time point on the plot')
        minClick = plt.ginput(1)
        endClick = plt.ginput(1)
        plt.close()
        minExtract = minClick[0]
        (begin, fThresh) = minExtract
        begin = round(begin)
        endExtract = endClick[0]
        (finish, y) = endExtract
        finish = round(finish)
        
        
        
        forceTot = dat.Force
        forceTot[forceTot<fThresh] = 0
        #dat.forceFilt = forceTot
        
        
    
       #find the landings and offs of the FP as vectors
        landings = findLandings(forceTot)
        landings[:] = [x for x in landings if x > begin] # remove landings before the start point specified on graph
        takeoffs = findTakeoffs(forceTot)
        takeoffs[:] = [x for x in takeoffs if x > landings[0]] # remove takeoffs before first landing
        landings[:] = [x for x in landings if x < finish] # remove landings after end of usable data

        for step in range(len(landings)):
            try:
                
                
                ic = landings[step]
                to = takeoffs[step]
                msStart = landings[step] + round((takeoffs[step] - landings[step]) *0.4)
                msEnd = landings[step] + round((takeoffs[step] - landings[step]) * 0.6)
                icStart = landings[step] + round((takeoffs[step] - landings[step]) *0.1)
                icEnd = landings[step] + round((takeoffs[step] - landings[step]) * 0.3)
                toStart = landings[step] + round((takeoffs[step] - landings[step]) *0.7)
                toEnd = landings[step] + round((takeoffs[step] - landings[step]) * 0.9)
                
                stanceForce = np.mean(forceTot[ic:to])
                icForce = np.mean(forceTot[icStart:icEnd])
                msForce = np.mean(forceTot[msStart:msEnd])
                toForce = np.mean(forceTot[toStart:toEnd])
                
                #cvHeel.append(np.std(dat.RHeel_MeanP[ic:to])/np.mean(dat.RHeel_MeanP[ic:to]))  
                #meanHeel.append(np.mean(dat.RHeel_MeanP[ic:to]))
                icHeel.append((np.mean(dat.RHeel_MeanP[icStart:icEnd]))/icForce)
                msHeel.append((np.mean(dat.RHeel_MeanP[msStart:msEnd]))/msForce)
                toHeel.append((np.mean(dat.RHeel_MeanP[toStart:toEnd]))/toForce)
                
                
                #cvLatFF.append(np.std(dat.RLatFF_MeanP[ic:to])/np.mean(dat.RLatFF_MeanP[ic:to]))
                #meanLatFF.append(np.mean(dat.RLatFF_MeanP[ic:to]))
                icLatFF.append((np.mean(dat.RLatFF_MeanP[icStart:icEnd]))/icForce)
                msLatFF.append((np.mean(dat.RLatFF_MeanP[msStart:msEnd]))/msForce)
                toLatFF.append((np.mean(dat.RLatFF_MeanP[toStart:toEnd]))/toForce)
                
                             
                #cvMedFF.append(np.std(dat.RMedFF_MeanP[ic:to])/np.mean(dat.RMedFF_MeanP[ic:to]))
                #meanMedFF.append(np.mean(dat.RMedFF_MeanP[ic:to]))
                icMedFF.append((np.mean(dat.RMedFF_MeanP[icStart:icEnd]))/icForce)
                msMedFF.append((np.mean(dat.RMedFF_MeanP[msStart:msEnd]))/msForce)
                toMedFF.append((np.mean(dat.RMedFF_MeanP[toStart:toEnd]))/toForce)
                
                #cvToes.append(np.std(dat.RToes_MeanP[ic:to])/np.mean(dat.RToes_MeanP[ic:to]))
                #meanToes.append(np.mean(dat.RToes_MeanP[ic:to]))
                icToes.append((np.mean(dat.RToes_MeanP[icStart:icEnd]))/icForce)
                msToes.append((np.mean(dat.RToes_MeanP[msStart:msEnd]))/msForce)
                toToes.append((np.mean(dat.RToes_MeanP[toStart:toEnd]))/toForce)
                
                # cvDorsLatFF.append(np.std(dat.DorsLatFF_MeanP[ic:to])/np.mean(dat.DorsLatFF_MeanP[ic:to]))
                # meanDorsLatFF.append(np.mean(dat.DorsLatFF_MeanP[ic:to]))
                # cvDorsMedFF.append(np.std(dat.DorsMedFF_MeanP[ic:to])/np.mean(dat.DorsMedFF_MeanP[ic:to]))
                # meanDorsMedFF.append(np.mean(dat.DorsMedFF_MeanP[ic:to]))
                # cvDorsLatMF.append(np.std(dat.DorsLatMF_MeanP[ic:to])/np.mean(dat.DorsLatMF_MeanP[ic:to]))
                # meanDorsLatMF.append(np.mean(dat.DorsLatMF_MeanP[ic:to]))
                # cvDorsMedMF.append(np.std(dat.DorsMedMF_MeanP[ic:to])/np.mean(dat.DorsMedMF_MeanP[ic:to]))
                # meanDorsMedMF.append(np.mean(dat.DorsMedMF_MeanP[ic:to]))
                
                contactTime.append(to - ic)
                trial.append(fName)
                Subject.append(subName)
                Condition.append(ConditionTmp)
                Config.append(ConfigTmp)

                #the outcomes dataframe uses the same logic (first set are dorsal pressures until maxPlantMetP, 
                #then everything is plantar unless noted otherwise), then midstanceHeelP is the heel pressure at 
                #mid stance and the heel decay from peak (which is a bit finnicky at the moment). 
                        
                
            except:
                print(step)
        
        

    except:
        print(fName)       

outcomes = pd.DataFrame({'SubjectName':list(Subject),'Condition':list(Condition), 'Shoe':list(Config), 'ContactTime':list(contactTime),
                                         'icHeel':list(icHeel), 'msHeel':list(msHeel), 'toHeel':list(toHeel),
                                         'icLatFF':list(icLatFF), 'msLatFF':list(msLatFF), 'toLatFF':list(toLatFF),
                                         'icMedFF':list(icMedFF), 'msMedFF':list(msMedFF), 'toMedFF':list(toMedFF),
                                         'icToes':list(icToes), 'msToes':list(msToes), 'toToes':list(toToes),
                                         
                                         })

outcomes.to_csv("C:/Users/kate.harrison/Dropbox (Boa)/EndurancePerformance/PFL_ISB2021_HeelHoldStudy/MVA Files/Run/CompiledPressureData.csv", index = False)

#Removed Dorsal Variables from output
# 'cvDorsLatFF':list(cvDorsLatFF), 'meanDorsLatFF':list(meanDorsLatFF),
#                                          'cvDorsMedFF':list(cvDorsMedFF),'meanDorsMedFF':list(meanDorsMedFF),
#                                          'cvDorsLatMF':list(cvDorsLatMF), 'meanDorsLatMF':list(meanDorsLatMF),
#                                          'cvDorsMedMF':list(cvDorsMedMF),'meanDorsMedMF':list(meanDorsMedMF)

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
