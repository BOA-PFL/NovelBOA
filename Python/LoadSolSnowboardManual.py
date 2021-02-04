# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 16:19:30 2020

@author: Daniel.Feeney
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# Read in files
# only read .asc files for this work
fPath = 'C:/Users/Daniel.Feeney/Dropbox (Boa)/Snow Protocol/SnowboardProtocol/'
entries = os.listdir(fPath)


for file in entries:
    try:
        fName = file
        dat = pd.read_csv(fPath+fName,sep='         ', skiprows = 3, header = 0, index_col = False)
        dat.columns = ['Time', 'LHeel', 'LMedial','LLateral','LTotal', 'Time2', 'RLateral','RMedial','RHeel','RTotal']
        subName = fName.split(sep = "_")[0]
        configName = fName.split(sep = "_")[1]
        
        dat['LToes'] = dat.LMedial + dat.LLateral
        dat['RToes'] = dat.RMedial + dat.RLateral
        
        fig, ax = plt.subplots()
        ax.plot(dat.LTotal, label = 'Left Total Force')
        ax.plot(dat.RTotal, label = 'Right Total Force')
        fig.legend()
        print('Select start and end of analysis trial')
        pts = np.asarray(plt.ginput(2, timeout=15))
        plt.close()
        # downselect the region of the dataframe you'd like
        dat = dat.iloc[int(np.floor(pts[0,0])) : int(np.floor(pts[1,0])),:]
        
        #### Initial analysis plans: segment heel and toe turns, calcualte smoothness
        ## as the CV or SD of force during the turn, calcualte turn time, calculate symmetry
        
        # Find indices of toe and heel turns but leave original data untouched 
        heelThresh = 150 #below this value will be set to 0 temporarily to find indices to start/end turns
        toeThresh = 75
        
        tmpToes = np.array(dat.LToes)
        tmpHeel = np.array(dat.LHeel)
        
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
        
        ###### After finding starts of turns, find avg, SD, CV, etc. for each turn ####
        dat = dat.reset_index()
        turnToPlot = 3
        fwdLook = 200
        
        # fig, ax = plt.subplots(2)
        # fig.suptitle('Toe Turn')
        # ax[0].plot(dat.LToes[realToeStart[turnToPlot]:realToeStart[turnToPlot]+fwdLook], label='L Toes')
        # ax[0].plot(dat.RToes[realToeStart[turnToPlot]:realToeStart[turnToPlot]+fwdLook], label = 'R Toes')
        # ax[0].set_title('Toes')
        # ax[0].set_ylim([0,800])
        # fig.legend()
        # ax[1].plot(dat.LHeel[realToeStart[turnToPlot]:realToeStart[turnToPlot]+fwdLook], 'tab:green', label = 'L Heel')
        # ax[1].plot(dat.RHeel[realToeStart[turnToPlot]:realToeStart[turnToPlot]+fwdLook], 'tab:red', label = 'R Heel')
        # ax[1].set_title('Heel')
        # ax[1].set_ylim([0,800])
        # fig.legend()
        # fig.tight_layout()
        
        # fig, ax = plt.subplots(2)
        # fig.suptitle('Heel Turn')
        # ax[0].plot(dat.LToes[realHeelStart[turnToPlot]:realHeelStart[turnToPlot]+fwdLook], label='L Toes')
        # ax[0].plot(dat.RToes[realHeelStart[turnToPlot]:realHeelStart[turnToPlot]+fwdLook], label = 'R Toes')
        # ax[0].set_ylim([0,800])
        # fig.legend()
        # ax[1].plot(dat.LHeel[realHeelStart[turnToPlot]:realHeelStart[turnToPlot]+fwdLook], 'tab:green', label = 'L Heel')
        # ax[1].plot(dat.RHeel[realHeelStart[turnToPlot]:realHeelStart[turnToPlot]+fwdLook], 'tab:red', label = 'R Heel')
        # ax[1].set_ylim([0,800])
        # fig.legend()
        # fig.tight_layout()
        
        
        ###### Extract variables from each turn initiation ######
        # left toes #
        maxFL = [ np.max(dat.LToes[toeTurnStart:toeTurnStart+100]) for toeTurnStart in realToeStart ]
        maxRFDupL = [ np.max(dat.LToes[toeTurnStart:toeTurnStart+100].diff()) for toeTurnStart in realToeStart ]
        maxRFDdnL = [ np.min(dat.LToes[toeTurnStart:toeTurnStart+100].diff()) for toeTurnStart in realToeStart ]
        timeToPeakL = [ list(dat.LToes[toeTurnStart:toeTurnStart+100]).index(max(dat.LToes[toeTurnStart:toeTurnStart+100])) for toeTurnStart in realToeStart ]
        stdPeakL = [ np.std(dat.LToes[times-10:times+10]) for times in timeToPeakL ]
        noLeft = list(np.repeat('L', len(stdPeakL)))
        
        # Right Toes #
        maxFR = [ np.max(dat.RToes[toeTurnStart:toeTurnStart+100]) for toeTurnStart in realToeStart ]
        maxRFDupR = [ np.max(dat.RToes[toeTurnStart:toeTurnStart+100].diff()) for toeTurnStart in realToeStart ]
        maxRFDdnR = [ np.min(dat.RToes[toeTurnStart:toeTurnStart+100].diff()) for toeTurnStart in realToeStart ]
        timeToPeakR = [ list(dat.RToes[toeTurnStart:toeTurnStart+100]).index(max(dat.RToes[toeTurnStart:toeTurnStart+100])) for toeTurnStart in realToeStart ]
        stdPeakR = [ np.std(dat.RToes[times-10:times+10]) for times in timeToPeakL ]
        noRight = list(np.repeat('R', len(stdPeakR)))
        noToes = list(np.repeat('Toes', len(noLeft + noRight)))
        
        # left Heel #
        maxFLHeel = [ np.max(dat.LHeel[heelTurnStart:heelTurnStart+100]) for heelTurnStart in realHeelStart ]
        maxRFDupLHeel = [ np.max(dat.LHeel[heelTurnStart:heelTurnStart+100].diff()) for heelTurnStart in realHeelStart ]
        maxRFDdnLHeel = [ np.min(dat.LHeel[heelTurnStart:heelTurnStart+100].diff()) for heelTurnStart in realHeelStart ]
        timeToPeakLHeel = [ list(dat.LHeel[heelTurnStart:heelTurnStart+100]).index(max(dat.LHeel[heelTurnStart:heelTurnStart+100])) for heelTurnStart in realHeelStart ]
        stdPeakLHeel = [ np.std(dat.LHeel[times-10:times+10]) for times in timeToPeakLHeel ]
        noLeftHeel = list(np.repeat('L', len(stdPeakL)))
        
        # Right Heel #
        maxFRHeel = [ np.max(dat.RHeel[heelTurnStart:heelTurnStart+100]) for heelTurnStart in realHeelStart ]
        maxRFDupRHeel = [ np.max(dat.RHeel[heelTurnStart:heelTurnStart+100].diff()) for heelTurnStart in realHeelStart ]
        maxRFDdnRHeel = [ np.min(dat.RHeel[heelTurnStart:heelTurnStart+100].diff()) for heelTurnStart in realHeelStart ]
        timeToPeakRHeel = [ list(dat.RHeel[heelTurnStart:heelTurnStart+100]).index(max(dat.RHeel[heelTurnStart:heelTurnStart+100])) for heelTurnStart in realHeelStart ]
        stdPeakRHeel = [ np.std(dat.RHeel[times-10:times+10]) for times in timeToPeakRHeel ]
        noRightHeel = list(np.repeat('R', len(stdPeakR)))
        noHeel = list(np.repeat('Heel', len(noLeftHeel + noRightHeel)))
        
        # naming #
        totalLength = len(noLeft) + len(noRight) + len(noLeftHeel) + len(noRightHeel)
        longSubject = list( np.repeat(subName, totalLength) )
        longConfig = list( np.repeat(configName, totalLength) )
        
        outcomes = pd.DataFrame({'Subject':list(longSubject), 'Config': list(longConfig), 'Side':list(noLeft + noRight),
                                 'TurnType': list(noToes + noHeel),'MaxForceToes':list(maxFL + maxFR),
                                 'MaxRFDUp': list(maxRFDupL + maxRFDupR),'MaxRFDdn': list(maxRFDdnL + maxRFDdnR), 
                                 'timeToPeak':list(timeToPeakL + timeToPeakR),'stdPeak': list(stdPeakL + stdPeakR)})
        
        outcomes.to_csv('C:/Users/Daniel.Feeney/Dropbox (Boa)/Snow Protocol/SnowboardProtocol/Results/snowboardResults.csv', mode='a', header=False)
#plt.plot(dat.LToes[realToeStart[1]:realToeStart[1]+100])

#avgF2 = [movAvgForce(forceZ, landing, landing+100, 10) for toeTurnStart in realToeStart]

    except:
        print(file)