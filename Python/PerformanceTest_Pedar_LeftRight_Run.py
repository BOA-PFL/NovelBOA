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
# only read .mva files for this work
fPath = 'C:/Users/kate.harrison/Boa Technology Inc/PFL - Documents/General/AgilityPerformanceData/BOA_DualPanel_Sept2021/Novel/'
fileExt = r".mva"
entries = [fName for fName in os.listdir(fPath) if fName.endswith(fileExt)]

# list of functions 
# finding landings on the force plate once the filtered force exceeds the force threshold
def findLandings(force, fThresh):
    lic = []
    for step in range(len(force)-1):
        if force[step] == 0 and force[step + 1] >= fThresh:
            lic.append(step)
    return lic

#Find takeoff from FP when force goes from above thresh to 0
def findTakeoffs(force, fThresh):
    lto = []
    for step in range(len(force)-1):
        if force[step] >= fThresh and force[step + 1] == 0:
            lto.append(step + 1)
    return lto

#initialize arrays
sdHeel = []
meanToes = []
maxHeel = []
maxToes = []
cvHeel = []
meanTotalP = []
CT = []


Subject = []
Config = []

for fName in entries:
    try:

        #fName = entries[0]
        subName = fName.split(sep = "_")[0]
        ConfigTmp = fName.split(sep="_")[2]
        
        dat = pd.read_csv(fPath+fName,sep='\t', skiprows = 13, header = 0)
        
        dat.columns = ['Time', 'Heel_Force', 'Heel_MaxP', 'Heel_MeanP', 'Heel_Pct', 
                       'Midfoot_Force', 'Midfoot_MaxP', 'Midfoot_MeanP', 'Midfoot_Pct',
                       'Forefoot_Force', 'Forefoot_MaxP', 'Forefoot_MeanP', 'Forefoot_Pct',
                       'Toe_Force', 'Toe_MaxP', 'Toe_MeanP', 'Toe_Pct']
        
        dat['Force'] = dat.Heel_Force + dat.Midfoot_Force + dat.Forefoot_Force + dat.Toe_Force
        
        # delimit trial
        fig, ax = plt.subplots()
        ax.plot(dat.Force, label = 'Right Total Force')
        fig.legend()
        print('Select start and end of analysis trial. CLICK AT THE LEVEL OF THE APPROPRIATE FORCE THRESHOLD')
        pts = np.asarray(plt.ginput(2, timeout=-1))
        plt.close()
        fThresh = pts[0,1]
        dat = dat.iloc[int(np.floor(pts[0,0])) : int(np.floor(pts[1,0])),:]
        dat = dat.reset_index()
        
        
        dat.Force[dat.Force<fThresh] = 0
        
               
        
       #find the landings and offs of the FP as vectors
        landings = findLandings(dat.Force, fThresh)
        takeoffs = findTakeoffs(dat.Force, fThresh)

        landings[:] = [x for x in landings if x < takeoffs[-1]]
        takeoffs[:] = [x for x in takeoffs if x > landings[0]]
        
        for counterVar, landing in enumerate(landings):
            try:
                    
                    meanHeel = np.mean(dat.Heel_MeanP[landing:takeoffs[counterVar]]) 
                    meanMidfoot = np.mean(dat.Midfoot_MeanP[landing:takeoffs[counterVar]])    
                    meanForefoot = np.mean(dat.Forefoot_MeanP[landing:takeoffs[counterVar]])
                    meanToe = np.mean(dat.Toe_MeanP[landing:takeoffs[counterVar]])  
                    
                    stdevHeel = np.std(dat.Heel_MeanP[landing:takeoffs[counterVar]])
                    maximumHeel = np.max(dat.Heel_MeanP[landing:takeoffs[counterVar]])
                    maximumToe = np.max(dat.Toe_MeanP[landing:takeoffs[counterVar]])
                    meanFoot = (meanHeel + meanMidfoot + meanForefoot + meanToe)/4
                    
                    meanTotalP.append(meanFoot)
                    sdHeel.append(stdevHeel)
                    meanToes.append(meanToe/meanFoot)
                    
                    maxHeel.append(maximumHeel)
                    cvHeel.append(stdevHeel/meanFoot)
                    
                    maxToes.append(maximumToe/meanFoot)
                    CT.append(takeoffs[counterVar] - landing)
            
                
                    Subject.append(subName)
                    Config.append(ConfigTmp)

                
                        

            except:
                print(fName, landing)
        

    except:
        print(fName)       

outcomes = pd.DataFrame({'Subject':list(Subject),'Config':list(Config),'meanTotalP':list(meanTotalP),
                         'sdHeel': list(sdHeel),'cvHeel':list(cvHeel),
                         'meanToes':list(meanToes), 
                         'maxHeel':list(maxHeel), 'maxToes':list(maxToes),
                         'CT':list(CT)
                         })

outFileName = fPath + 'CompiledPressureData.csv'
outcomes.to_csv(outFileName, index = False)

