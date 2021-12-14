# -*- coding: utf-8 -*-
"""
Created on Fri Mar 19 08:55:02 2021
For when only the right loadsol records in snowboarding
@author: Daniel.Feeney
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# Read in files
# only read .asc files for this work
fPath = 'C:/Users/Daniel.Feeney/Dropbox (Boa)/Snow Protocol/SnowboardProtocol/'
fPath = 'C:/Users/Daniel.Feeney/Dropbox (Boa)/Snow Protocol/4onthefloor/loadsoldata/'
entries = os.listdir(fPath)


for file in entries[0:2]:
    try:
        fName = file
        dat = pd.read_csv(fPath+fName,sep='         ', skiprows = 3, header = 0, index_col = False)
        dat.columns = ['Time', 'RLateral','RMedial','RHeel','RTotal']
        subName = fName.split(sep = "_")[0]
        configName = fName.split(sep = "_")[1]
        
        dat['RToes'] = dat.RMedial + dat.RLateral
        
        fig, ax = plt.subplots()
        ax.plot(dat.RTotal, label = 'Right Total Force')
        fig.legend()
        print('Select start and end of analysis trial')
        pts = np.asarray(plt.ginput(2, timeout=-1))
        plt.close()
        # downselect the region of the dataframe you'd like
        dat = dat.iloc[int(np.floor(pts[0,0])) : int(np.floor(pts[1,0])),:]
        
        #### Initial analysis plans: segment heel and toe turns, calcualte smoothness
        ## as the CV or SD of force during the turn, calcualte turn time, calculate symmetry
        
        # Find indices of toe and heel turns but leave original data untouched 
        heelThresh = 150 #below this value will be set to 0 temporarily to find indices to start/end turns
        toeThresh = 75
        
        tmpToes = np.array(dat.RToes)
        tmpHeel = np.array(dat.RHeel)
        
        tmpHeel[tmpHeel < heelThresh] = 0
        tmpToes[tmpToes < toeThresh] = 0
        
        # Select toe turn starts
        plt.plot(tmpToes, label='Toes')
        plt.plot(tmpHeel, label='Heel')
        plt.legend()
        print('Select start of Toe Turns')
        toestart = np.asarray(plt.ginput(5, timeout=-1))
        realToeStart = [int(toestart[0,0]), int(toestart[1,0]), (int(toestart[2,0])), int(toestart[3,0]),int(toestart[4,0])]
        plt.close()
        
        # Select heel turn starts
        plt.plot(tmpToes, label='Toes')
        plt.plot(tmpHeel, label='Heel')
        plt.legend()
        print('Select start of heel turns')
        heelstart = np.asarray(plt.ginput(5, timeout=-1))
        realHeelStart = [int(heelstart[0,0]), int(heelstart[1,0]), (int(heelstart[2,0])), int(heelstart[3,0]),int(heelstart[4,0])]
        plt.close()   
        
        ###### Extract variables from each turn initiation ######
        
        # Right Toes #
        maxFR = [ np.max(dat.RToes[toeTurnStart:toeTurnStart+100]) for toeTurnStart in realToeStart ]
        maxRFDupR = [ np.max(dat.RToes[toeTurnStart:toeTurnStart+100].diff()) for toeTurnStart in realToeStart ]
        maxRFDdnR = [ np.min(dat.RToes[toeTurnStart:toeTurnStart+100].diff()) for toeTurnStart in realToeStart ]
        timeToPeakR = [ list(dat.RToes[toeTurnStart:toeTurnStart+100]).index(max(dat.RToes[toeTurnStart:toeTurnStart+100])) for toeTurnStart in realToeStart ]
        stdPeakR = [ np.std(dat.RToes[times-10:times+10]) for times in timeToPeakR ]
        noRight = list(np.repeat('R', len(stdPeakR)))
        noToes = list(np.repeat('Toes', len(noRight)))
                
        # Right Heel #
        maxFRHeel = [ np.max(dat.RHeel[heelTurnStart:heelTurnStart+100]) for heelTurnStart in realHeelStart ]
        maxRFDupRHeel = [ np.max(dat.RHeel[heelTurnStart:heelTurnStart+100].diff()) for heelTurnStart in realHeelStart ]
        maxRFDdnRHeel = [ np.min(dat.RHeel[heelTurnStart:heelTurnStart+100].diff()) for heelTurnStart in realHeelStart ]
        timeToPeakRHeel = [ list(dat.RHeel[heelTurnStart:heelTurnStart+100]).index(max(dat.RHeel[heelTurnStart:heelTurnStart+100])) for heelTurnStart in realHeelStart ]
        stdPeakRHeel = [ np.std(dat.RHeel[times-10:times+10]) for times in timeToPeakRHeel ]
        noRightHeel = list(np.repeat('R', len(stdPeakR)))
        noHeel = list(np.repeat('Heel', len(noRightHeel)))
        
        # naming #
        totalLength =  len(noRight)+ len(noRightHeel)
        longSubject = list( np.repeat(subName, totalLength) )
        longConfig = list( np.repeat(configName, totalLength) )
        
        outcomes = pd.DataFrame({'Subject':list(longSubject), 'Config': list(longConfig), 'Side':list( noRight + noRightHeel),
                                 'TurnType': list(noToes + noHeel),'MaxForceToes':list( maxFR  + maxFRHeel),
                                 'MaxRFDUp': list(maxRFDupR + maxRFDupRHeel),
                                 'MaxRFDdn': list(maxRFDdnR  + maxRFDdnRHeel), 
                                 'timeToPeak':list(timeToPeakR  + timeToPeakRHeel),
                                 'stdPeak': list(stdPeakR + stdPeakRHeel)})
        
        outcomes.to_csv('C:/Users/Daniel.Feeney/Dropbox (Boa)/Snow Protocol/4onthefloor/Results.csv', mode='a', header=False)
#plt.plot(dat.LToes[realToeStart[1]:realToeStart[1]+100])

#avgF2 = [movAvgForce(forceZ, landing, landing+100, 10) for toeTurnStart in realToeStart]

    except:
        print(file)
     
    
