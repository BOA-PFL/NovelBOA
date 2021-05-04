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
fPath = 'C:\\Users\\Daniel.Feeney\\Dropbox (Boa)\\Hike Work Research\\Work Pilot 2021\\Pressures\\'
fileExt = r".mva"
entries = [fName for fName in os.listdir(fPath) if fName.endswith(fileExt)]

# Define constants and options
#fThresh = 0 #below this value will be set to 0.
stepLen = 45 #Set value to look forward 
autoDetectTakeoff = 1 #if this is 1, it will try to find landings and takeoffs. 0 means it will look forward the step length
plottingEnabled = 1

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


#first columns (FF, Mets, and MF) all relate to dorsal values. Once you get to PlantMetsForce it is plantar metatarsal force 
#and everything to the right of that column is plantar side. Each location (e.g. FF, MF, etc.) has force, max Pressure, Mean Pressure, and pct

for file in entries:
    try:
 
        fName = file 
        #print(file)
        
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
        
        sdDorsalLatMF = []
        sdDorsalMedMF = []
        
        CT = []
        
        trial = []
        Subject = []
        Condition = []
        Config = []
        
        
        subName = fName.split(sep = "_")[0]
        ConditionTmp = fName.split(sep="_")[2]
        ConfigTmp = fName.split(sep="_")[1]
        
        dat = pd.read_csv(fPath+fName,sep='\t', skiprows = 16, header = 0)
        
        dat.columns = ['Time','RHeel_Force', 'RHeel_MaxP', 'RHeel_MeanP', 'RHeel_pct', 
                       'RLatFF_Force','RLatFF_MaxP', 'RLatFF_MeanP', 'RLatFF_pct', 
                       'RMedFF_Force', 'RMedFF_MaxP', 'RMedFF_MeanP','RMedFF_pct',
                       'RToes_Force','RToes_MaxP', 'RToes_MeanP', 'RToes_pct', 
                       'RDMedFF_Force', 'RDMedFF_MaxP', 'RDMedFF_MeanP','RDMedFF_pct',
                       'RDLatFF_Force','RDLatFF_MaxP', 'RDLatFF_MeanP', 'RDLatFF_Pct',
                       'RDMedMF_Force','RDMedMF_MaxP','RDMedMF_MeanP','RDMedMF_Pct',
                       'RDLatMF_Force','RDLatMF_MaxP','RDLatMF_MeanP','RDLatMF_Pct']
        
        # dat.columns = ['Time','RHeel_Force', 'RHeel_MaxP', 'RHeel_MeanP', 'RHeel_pct', 
        #        'RMedMF_Force','RMedMF_MaxP', 'RMedMF_MeanP', 'RMedMF_pct', 
        #        'RLatMF_Force', 'RLatMF_MaxP', 'RLatMF_MeanP','RLatMF_pct',
        #        'RMedFF_Force','RMedFF_MaxP', 'RMedFF_MeanP', 'RMedFF_pct', 
        #        'RLatFF_Force', 'RLatFF_MaxP', 'RLatFF_MeanP','RLatFF_pct',
        #        'RToes_Force','RToes_MaxP', 'RToes_MeanP', 'RToes_Pct']
        
        dat['Force'] = dat.RHeel_Force + dat.RLatFF_Force + dat.RMedFF_Force + dat.RToes_Force
        #filtering force to find landings/takeoffs    
        # delimit trial
        fig, ax = plt.subplots()
        ax.plot(dat.Force, label = 'Right Total Force')
        fig.legend()
        print('Select start and end of analysis trial')
        pts = np.asarray(plt.ginput(2, timeout=-1))
        plt.close()
        # downselect the region of the dataframe you'd like
        dat = dat.iloc[int(np.floor(pts[0,0])) : int(np.floor(pts[1,0])),:]
        dat = dat.reset_index()
        
        # find threshold force
        fig, ax = plt.subplots()
        ax.plot(dat.Force, label = 'Right Foot Force')
        print('Select a point to represent 0 in the trial')
        pts = np.asarray(plt.ginput(1, timeout=-1))
        plt.close()
        fThresh = pts[0][1]
        # downselect the region of the dataframe you'd like

        forceTot = dat.Force
        forceTot[forceTot<fThresh] = 0
        forceTot = np.array(forceTot)
        # plt.plot(forceTot)
        # plt.close()
        #dat.forceFilt = forceTot
        
        
       #find the landings and offs of the FP as vectors
        landings = findLandings(forceTot)
        takeoffs = findTakeoffs(forceTot)


        for counterVar, landing in enumerate(landings):
            try:
                
                                # Without takeoffs enabled
                if autoDetectTakeoff:
                    sdRHeel.append(np.std(dat.RHeel_MeanP[landing:takeoffs[counterVar]])) 
                    meanRHeel.append(np.mean(dat.RHeel_MeanP[landing:takeoffs[counterVar]]))
                    
                    sdRLatFF.append(np.std(dat.RLatFF_MeanP[landing:takeoffs[counterVar]]))
                    meanRLatFF.append(np.mean(dat.RLatFF_MeanP[landing:takeoffs[counterVar]]))
                    
                    sdRMedFF.append(np.std(dat.RMedFF_MeanP[landing:takeoffs[counterVar]]))
                    meanRMedFF.append(np.mean(dat.RMedFF_MeanP[landing:takeoffs[counterVar]]))
                    
                    sdRToes.append(np.std(dat.RToes_MeanP[landing:takeoffs[counterVar]]))
                    meanRToes.append(np.mean(dat.RToes_MeanP[landing:takeoffs[counterVar]]))
                    
                    maxRHeel.append(np.max(dat.RHeel_MaxP[landing:takeoffs[counterVar]]))
                    maxRLatFF.append(np.max(dat.RLatFF_MaxP[landing:takeoffs[counterVar]]))
                    maxRMedFF.append(np.max(dat.RMedFF_MaxP[landing:takeoffs[counterVar]]))
                    maxRToes.append(np.max(dat.RToes_MaxP[landing:takeoffs[counterVar]]))
                    # Dorsal values
                    sdDorsalLatMF.append(np.std(dat.RDLatMF_MeanP[landing:takeoffs[counterVar]]))
                    sdDorsalMedMF.append(np.std(dat.RDMedMF_MeanP[landing:takeoffs[counterVar]]))
                    
                    CT.append(takeoffs[counterVar] - landing)
                else:
                
                    sdRHeel.append(np.std(dat.RHeel_MeanP[landing:landing+stepLen])) 
                    meanRHeel.append(np.mean(dat.RHeel_MeanP[landing:landing+stepLen]))
                    
                    sdRLatFF.append(np.std(dat.RLatFF_MeanP[landing:landing+stepLen]))
                    meanRLatFF.append(np.mean(dat.RLatFF_MeanP[landing:landing+stepLen]))
                    
                    sdRMedFF.append(np.std(dat.RMedFF_MeanP[landing:landing+stepLen]))
                    meanRMedFF.append(np.mean(dat.RMedFF_MeanP[landing:landing+stepLen]))
                    
                    sdRToes.append(np.std(dat.RToes_MeanP[landing:landing+stepLen]))
                    meanRToes.append(np.mean(dat.RToes_MeanP[landing:landing+stepLen]))
                    
                    maxRHeel.append(np.max(dat.RHeel_MaxP[landing:landing+stepLen]))
                    maxRLatFF.append(np.max(dat.RLatFF_MaxP[landing:landing+stepLen]))
                    maxRMedFF.append(np.max(dat.RMedFF_MaxP[landing:landing+stepLen]))
                    maxRToes.append(np.max(dat.RToes_MaxP[landing:landing+stepLen]))
               
                    # Add std and average dorsal midfoot values
                    sdDorsalLatMF.append(np.std(dat.RDLatMF_MeanP[landing:landing+stepLen]))
                    sdDorsalMedMF.append(np.std(dat.RDMedMF_MeanP[landing:landing+stepLen]))
                    
                    CT.append(stepLen)
               
                trial.append(fName)
                Subject.append(subName)
                Condition.append(ConditionTmp)
                Config.append(ConfigTmp)

                #the outcomes below are entirely on the plantar aspect of the foot
                        
                outcomes = pd.DataFrame({'Subject':list(Subject),'Condition':list(Condition), 'Config':list(Config),
                                         'sdRHeel': list(sdRHeel),'meanRHeel':list(meanRHeel), 
                                         'sdRLatFF':list(sdRLatFF), 'meanRLatFF':list(meanRLatFF),
                                         'sdRMedFF':list(sdRMedFF),'meanRMedFF':list(meanRMedFF),
                                         'sdRToes':list(sdRToes),'meanRToes':list(meanRToes), 
                                         'maxRHeel':list(maxRHeel), 'maxRLatMF':list(maxRLatFF), 
                                         'maxRMedMF':list(maxRMedFF), 'maxRLatFF':list(maxRLatFF), 
                                         'maxRMedFF':list(maxRMedFF), 'maxRToes':list(maxRToes),
                                         'sdDorsalLatMF':list(sdDorsalLatMF), 'sdDorsalMedMF':list(sdDorsalMedMF),
                                         'CT':list(CT)
                                         })
            except:
                print(landing)
        
        outcomes.to_csv("C:\\Users\\Daniel.Feeney\\Dropbox (Boa)\\Hike Work Research\\Work Pilot 2021/pressuresComb3.csv", mode= 'a', header=False)

    except:
        print(file)       


## plotting
if plottingEnabled == 1:
    landingToPlot = 2
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    
    color = 'tab:red'
    ax1.set_xlabel('time')
    ax1.set_ylabel('TotalForce(N)', color=color)
    ax1.plot(dat.Force[landings[landingToPlot]:takeoffs[landingToPlot]], color=color, label = 'Total Force')
    ax1.tick_params(axis='y', labelcolor=color)
    
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    
    ax2.set_ylabel('Dorsal Max Pressures')  # we already handled the x-label with ax1
    ax2.plot(dat.RDLatFF_MaxP[landings[landingToPlot]:takeoffs[landingToPlot]], label = 'Max Lat FF Pressure')
    ax2.plot(dat.RDMedFF_MaxP[landings[landingToPlot]:takeoffs[landingToPlot]], label = 'Max Med FF Pressure')
    ax2.plot(dat.RDLatMF_MaxP[landings[landingToPlot]:takeoffs[landingToPlot]], label = 'Max Lat MF Pressure')
    ax2.plot(dat.RDMedMF_MaxP[landings[landingToPlot]:takeoffs[landingToPlot]], label = 'Max Med MF Pressure')
    # ask matplotlib for the plotted objects and their labels
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1+h2, l1+l2, loc=2)
    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.show()
    #
    
    ## mean pressure
    fig2 = plt.figure()
    ax1 = fig2.add_subplot(111)
    
    color = 'tab:red'
    ax1.set_xlabel('time')
    ax1.set_ylabel('TotalForce(N)', color=color)
    ax1.plot(dat.Force[landings[landingToPlot]:takeoffs[landingToPlot]], color=color, label = 'Total Force')
    ax1.tick_params(axis='y', labelcolor=color)
    
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    
    ax2.set_ylabel('Dorsal Mean Pressures (kPa)')  # we already handled the x-label with ax1
    ax2.plot(dat.RDLatFF_MeanP[landings[landingToPlot]:takeoffs[landingToPlot]], label = 'Mean Lat FF Pressure')
    ax2.plot(dat.RDMedFF_MeanP[landings[landingToPlot]:takeoffs[landingToPlot]], label = 'Mean Med FF Pressure')
    ax2.plot(dat.RDLatMF_MeanP[landings[landingToPlot]:takeoffs[landingToPlot]], label = 'Mean Lat MF Pressure')
    ax2.plot(dat.RDMedMF_MeanP[landings[landingToPlot]:takeoffs[landingToPlot]], label = 'Mean Med MF Pressure')
    # ask matplotlib for the plotted objects and their labels
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1+h2, l1+l2, loc=2)
    fig2.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.show()
    
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
