# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 16:19:30 2020
When you have both sides of data from loadsol
@author: Daniel.Feeney
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import os
from tkinter.filedialog import askopenfilenames


def findToeTurns(toeForce, heelForce):
    toeTurns = []
    for step in range(len(toeForce)-1):
        if toeForce[step] <= heelForce[step] and toeForce[step + 1] > heelForce[step + 1]:
            toeTurns.append(step)
    return toeTurns

#Find takeoff from FP when force goes from above thresh to 0
def findHeelTurns(toeForce, heelForce):
    heelTurns = []
    for step in range(len(toeForce)-1):
        if heelForce[step] <= toeForce[step] and heelForce[step + 1] > toeForce [step + 1]:
            heelTurns.append(step)
    return heelTurns


# Read in files
# only read .asc files for this work
fPath = 'C:/Users/kate.harrison/Boa Technology Inc/PFL - Documents/General/Segments/Snow Performance/SB_2DialTakeDown_Mar2022/Forces/'

# entries = os.listdir(fPath)
# Select files for a single subject
entries = askopenfilenames(initialdir = fPath)

### Time Series

toeTurns_FrontHeel = []
toeTurns_FrontToes = []
toeTurns_RearHeel = []
toeTurns_RearToes = []
heelTurns_FrontHeel = []
heelTurns_FrontToes = []
heelTurns_RearHeel = []
heelTurns_RearToes = []

### Toe Turns
maxToeF_BothToes = []
maxTotalF_BothToes = []
maxToeRFDup_BothToes = []
maxToeRFDdn_BothToes = []
maxTotalRFDup_BothToes = []
maxTotalRFDdn_BothToes = []
timeToToePeak_BothToes = []
timeToTotalPeak_BothToes = []
avgToeRFD_BothToes = []
avgTotalRFD_BothToes = []
meanHeelContact_BothToes = []
noPeaks_FrontToe = []
peakProminence_FrontToe = []
noPeaks_RearToe = []
peakProminence_RearToe = []


# naming #
Subject = []
Config = []
Trial = []

for fName in entries:
    try:
        #fName = entries[0] 
        dat = pd.read_csv(fName,sep='\t', skiprows = 4, header = None, index_col = False, usecols = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
        dat.columns = ['Time', 'FrontHeel', 'FrontMedial','FrontLateral','FrontTotal', 'Time2', 'RearLateral','RearMedial','RearHeel','RearTotal']
        info = fName.split(sep = "/")[-1]
        subName = info.split(sep = "_")[0]
        stance = info.split(sep = "_")[1]
        configName = info.split(sep = "_")[2]
        trialNo = info.split(sep = "_")[3]
        
        if stance == 'Goofy':
            dat.rename(columns = {'FrontHeel':'RearHeel', 'FrontMedial':'RearMedial', 'FrontLateral':'RearLateral', 'FrontTotal':'RearTotal', 'RearLateral':'FrontLateral', 'RearMedial':'FrontMedial', 'RearHeel':'FrontHeel','RearTotal':'FrontTotal'}, inplace = True)
        
        dat = dat.apply(pd.to_numeric)
        dat['FrontToes'] = dat.FrontMedial + dat.FrontLateral
        dat['RearToes'] = dat.RearMedial + dat.RearLateral
        
        dat['bothToes'] = dat.FrontToes + dat.RearToes
        dat['bothHeels'] = dat.FrontHeel + dat.RearHeel
        fig, ax = plt.subplots()
        ax.plot(dat.FrontTotal, label = 'Front Total Force')
        ax.plot(dat.RearTotal, label = 'Rear Total Force')
        fig.legend()
        print('Select start and end of analysis trial')
        pts = np.asarray(plt.ginput(2, timeout=-1))
        plt.close()
        # downselect the region of the dataframe you'd like
        dat = dat.iloc[int(np.floor(pts[0,0])) : int(np.floor(pts[1,0])),:]
        
        dat['TotalFront'] = dat.FrontHeel + dat.FrontToes
        dat['TotalRear'] = dat.RearHeel + dat.RearToes
        
        dat['Total'] = dat.TotalFront + dat.TotalRear
        
        ## Filter data for more accurate RFD)
        fs = 100 
        fc = 2
        w = fc / (fs / 2)
        b, a = signal.butter(4, w, 'low')
        dat['FrontHeel_Filt'] = signal.filtfilt(b, a, dat.FrontHeel)
        dat.FrontHeel_Filt[dat.FrontHeel_Filt<0] = 0
        
        dat['FrontToes_Filt'] = signal.filtfilt(b, a, dat.FrontToes)
        dat.FrontToes_Filt[dat.FrontToes_Filt<0] = 0
        
        dat['RearHeel_Filt'] = signal.filtfilt(b, a, dat.RearHeel)
        dat.RearHeel_Filt[dat.RearHeel_Filt<0] = 0
        
        dat['RearToes_Filt'] = signal.filtfilt(b, a, dat.RearToes)
        dat.RearToes_Filt[dat.RearToes_Filt<0] = 0
        
        dat['TotalFront_Filt'] = signal.filtfilt(b, a, dat.TotalFront)
        dat.TotalFront_Filt[dat.TotalFront_Filt<0] = 0
        
        dat['TotalRear_Filt'] = signal.filtfilt(b, a, dat.TotalRear)
        dat.TotalRear_Filt[dat.TotalRear_Filt<0] = 0
        
        dat['bothToes_Filt'] = signal.filtfilt(b, a, dat.bothToes)
        dat.bothToes_Filt[dat.bothToes_Filt<0] = 0
        
        dat['bothHeels_Filt'] = signal.filtfilt(b, a, dat.bothHeels)
        dat.bothHeels_Filt[dat.bothHeels_Filt<0] = 0
        
        dat['Total_Filt'] = signal.filtfilt(b, a, dat.Total)
        dat.Total_Filt[dat.Total_Filt<0] = 0
        
        dat['ToeProp'] = dat.bothToes_Filt/dat.Total_Filt
        
        # plt.plot(dat.TotalFront, label = 'Total Front Foot Force')
        # plt.plot(dat.FrontToes_Filt, label = 'Front Toe Force')
        # plt.plot(dat.FrontHeel_Filt, label = 'Front Heel Force')
        # plt.legend()
      
        
        # plt.plot(dat.bothToes_Filt)
        # plt.plot(dat.bothHeels_Filt)
        
        tmpToes = np.array(dat.bothToes_Filt)
        tmpHeel = np.array(dat.bothHeels_Filt)
        
        realToeStart = findToeTurns(tmpToes, tmpHeel)
        realHeelStart = findHeelTurns(tmpToes, tmpHeel)
        
        realHeelStart[:] = [x for x in realHeelStart if x > realToeStart[0]] # we want first heel turn after first toe turn
        realToeStart[:] = [x for x in realToeStart if x < realHeelStart[-1]] # we want to end on a heel turn
        
                  
        ###### Extract variables from each turn initiation ######
        
        for i in range(len(realToeStart)-1):
            
            #i = 10
            turnTime = realHeelStart[i] - realToeStart[i] 
            
            if turnTime >= 100:
                
                try:
                    #i = 1
                    ### Time Series
            
                    toeTurns_FrontHeel.append(dat.FrontHeel[realToeStart[i]:realToeStart[i] + 100])
                    toeTurns_FrontToes.append(dat.FrontToes[realToeStart[i]:realToeStart[i] + 100])
                    toeTurns_RearHeel.append(dat.RearHeel[realToeStart[i]:realToeStart[i] + 100])
                    toeTurns_RearToes.append(dat.RearToes[realToeStart[i]:realToeStart[i] + 100])
                    heelTurns_FrontHeel.append(dat.FrontHeel[realHeelStart[i]:realHeelStart[i] + 100])
                    heelTurns_FrontToes.append(dat.FrontToes[realHeelStart[i]:realHeelStart[i] + 100])
                    heelTurns_RearHeel.append(dat.RearHeel[realHeelStart[i]:realHeelStart[i] + 100])
                    heelTurns_RearToes.append(dat.RearToes[realHeelStart[i]:realHeelStart[i] + 100])
            
            
                    ### Toe Turns
                
                    maxToeF_BothToes.append(np.max(dat.bothToes_Filt[realToeStart[i]:realToeStart[i]+100])/np.max(dat.Total_Filt[realToeStart[i]:realToeStart[i]+100]))
                    maxTotalF_BothToes.append(np.max(dat.Total_Filt[realToeStart[i]:realToeStart[i]+100]))
                    maxToeRFDup_BothToes.append(np.nanmax(dat.bothToes_Filt[realToeStart[i]:realToeStart[i]+100].diff())) 
                    maxToeRFDdn_BothToes.append(np.nanmin(dat.bothToes_Filt[realToeStart[i]:realToeStart[i]+100].diff())) 
                    maxTotalRFDup_BothToes.append(np.nanmax(dat.Total_Filt[realToeStart[i]:realToeStart[i]+100].diff())) 
                    maxTotalRFDdn_BothToes.append(np.nanmin(dat.Total_Filt[realToeStart[i]:realToeStart[i]+100].diff())) 
                    
                    try:
                        timeToToePeak_BothToes.append(list(dat.ToeProp[realToeStart[i]:realToeStart[i]+100]).index(max(dat.ToeProp[realToeStart[i]:realToeStart[i]+100])))
                    except:
                        timeToToePeak_BothToes.append('nan')
                    
                    try:
                        timeToTotalPeak_BothToes.append(list(dat.Total_Filt[realToeStart[i]:realToeStart[i]+100]).index(max(dat.Total_Filt[realToeStart[i]:realToeStart[i]+100])))
                    except:
                        timeToTotalPeak_BothToes.append('nan')
                    
                    try:
                        avgToeRFD_BothToes.append(maxToeF_BothToes[-1]/timeToToePeak_BothToes[-1])
                    except:
                         avgToeRFD_BothToes.append('nan')
                         
                    try:
                        avgTotalRFD_BothToes.append(maxTotalF_BothToes[-1]/timeToTotalPeak_BothToes[-1])
                    except:
                        avgTotalRFD_BothToes.append('nan') 
                        
                    try:
                        meanHeelContact_BothToes.append((np.mean(dat.bothHeels_Filt[realToeStart[i]:realToeStart[i]+100]))/np.mean(dat.Total_Filt[realToeStart[i]:realToeStart[i]+100]))
                    except:
                        meanHeelContact_BothToes.append('nan')
                        
                        
                    peaks = signal.find_peaks(dat.FrontHeel[realToeStart[i]:realToeStart[i]+100])[0]
                    noPeaks_FrontToe.append(len(peaks))
                    prom = signal.peak_prominences(dat.FrontHeel[realToeStart[i]:realToeStart[i]+100], peaks)[0]
                    peakProminence_FrontToe.append(np.mean(signal.peak_prominences(dat.FrontHeel[realToeStart[i]:realToeStart[i]+100], peaks)[0]))
                    
                    peaks = signal.find_peaks(dat.RearHeel[realToeStart[i]:realToeStart[i]+100])[0]
                    noPeaks_RearToe.append(len(peaks))
                    prom = signal.peak_prominences(dat.RearHeel[realToeStart[i]:realToeStart[i]+100], peaks)[0]
                    peakProminence_RearToe.append(np.mean(signal.peak_prominences(dat.RearHeel[realToeStart[i]:realToeStart[i]+100], peaks)[0]))
        
                    # naming #
                    Subject.append(subName)
                    Config.append(configName)
                    Trial.append(trialNo)
        
                except:
                    print(fName + " Turn Number " + str(i))

    except:
        print(fName)
     

outcomes = pd.DataFrame({'Subject':list(Subject), 'Config': list(Config), 'Trial':list(Trial),
                         
                         'MaxToeF_BothToes':list(maxToeF_BothToes), 'MaxTotalF_BothToes':list(maxTotalF_BothToes),
                         'maxToeRFDup_BothToes':list(maxToeRFDup_BothToes), 'maxToeRFDdn_BothToes':list(maxToeRFDdn_BothToes), 'maxTotalRFDup_BothToes':list(maxTotalRFDup_BothToes), 'maxTotalRFDdn_BothToes':list(maxTotalRFDdn_BothToes), 
                         'timeToToePeak_BothtToes':list(timeToToePeak_BothToes), 'timeToTotalPeak_BothToes':list(timeToTotalPeak_BothToes),
                         'avgToeRFD_BothToes':list(avgToeRFD_BothToes), 'avgTotalRFD_BothToes':list(avgTotalRFD_BothToes),
                         'meanHeelContact_BothToes':list(meanHeelContact_BothToes)
                         
                         })

outfileName = fPath + 'CompiledResults6.csv'

if os.path.exists(outfileName) == False:
    
    outcomes.to_csv(outfileName, mode='a', header=True, index = False)

else:
    outcomes.to_csv(outfileName, mode='a', header=False, index = False)



### Plotting time series

configs = np.unique(Config)
toeTurn_FrontHeel = np.stack(toeTurns_FrontHeel)
toeTurn_FrontToes = np.stack(toeTurns_FrontToes)
toeTurn_RearHeel = np.stack(toeTurns_RearHeel)
toeTurn_RearToes = np.stack(toeTurns_RearToes)
heelTurn_FrontHeel = np.stack(heelTurns_FrontHeel)
heelTurn_FrontToes = np.stack(heelTurns_FrontToes)
heelTurn_RearHeel = np.stack(heelTurns_RearHeel)
heelTurn_RearToes = np.stack(heelTurns_RearToes)

Config = np.array(Config)
toeTurn_FrontHeel = np.column_stack((Config, toeTurn_FrontHeel))
toeTurn_FrontToes = np.column_stack((Config, toeTurn_FrontToes))
toeTurn_RearHeel = np.column_stack((Config, toeTurn_RearHeel))
toeTurn_RearToes = np.column_stack((Config, toeTurn_RearToes))
heelTurn_FrontHeel = np.column_stack((Config, heelTurn_FrontHeel))
heelTurn_FrontToes = np.column_stack((Config, heelTurn_FrontToes))
heelTurn_RearHeel = np.column_stack((Config, heelTurn_RearHeel))
heelTurn_RearToes = np.column_stack((Config, heelTurn_RearToes))

plt.figure(1)
plt.title('ToeTurn_FrontFoot')

plt.figure(2)
plt.title('ToeTurn_RearFoot')

plt.figure(3)
plt.title('HeelTurn_FrontFoot')

plt.figure(4)
plt.title('HeelTurn_RearFoot')

for config in configs:
    
    config = configs[0]
    tmpToeTurn_FrontHeel = toeTurn_FrontHeel[toeTurn_FrontHeel[:,0] == config]
    toeTurn_FrontHeelmean = np.mean(toeTurn_FrontHeel, axis = 0)
    toeTurn_RearHeelmean = np.mean(toeTurn_RearHeel, axis = 0)
    toeTurn_FrontToesmean = np.mean(toeTurn_FrontToes, axis = 0)
    toeTurn_RearToesmean = np.mean(toeTurn_RearToes, axis = 0)

heelTurn_FrontHeelmean = np.mean(heelTurn_FrontHeel, axis = 0)
heelTurn_RearHeelmean = np.mean(heelTurn_RearHeel, axis = 0)
heelTurn_FrontToesmean = np.mean(heelTurn_FrontToes, axis = 0)
heelTurn_RearToesmean = np.mean(heelTurn_RearToes, axis = 0)

plt.figure(1)
plt.title('ToeTurn_FrontFoot')
plt.plot(toeTurn_FrontHeelmean, label = 'Heel_' + configName) 
plt.plot(toeTurn_FrontToesmean, label = 'Toes_' + configName)
plt.legend()