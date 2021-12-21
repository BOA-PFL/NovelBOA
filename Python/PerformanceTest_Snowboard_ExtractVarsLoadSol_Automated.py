# -*- coding: utf-8 -*-
"""
Created on Dec 20, 2021
When you have both sides of data from loadsol
@author: Daniel.Feeney

This is still a beta test of the automation based on some improvements
I found in the ski testing to find the peaks and use a backward and forward 
window to delimit the turns 
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.signal import find_peaks


# Read in files
# only read .asc files for this work
fPath = 'C:/Users/Daniel.Feeney/Dropbox (Boa)/Snow Protocol/4onthefloor/loadsoldata/'
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
        pts = np.asarray(plt.ginput(2, timeout=-1))
        plt.close()
        # downselect the region of the dataframe you'd like
        dat = dat.iloc[int(np.floor(pts[0,0])) : int(np.floor(pts[1,0])),:]
        dat = dat.reset_index()
        
        #### Initial analysis plans: segment heel and toe turns, calcualte smoothness
        ## as the CV or SD of force during the turn, calcualte turn time, calculate symmetry
        
        #### Working on peak and turn initiation detection algorithm here ####
        ### This will work on the toe turns initially and abstract to the left ## 
        plt.plot(dat.LToes)
        print('Select threshold force where peak forces will be above')
        fThresh = np.asarray(plt.ginput(1, timeout=-1))[0][1]
        plt.close()
        # select minimal force for turn initiation #
        plt.plot(dat.LToes)
        print('Select threshold force where min forces will be BELOW')
        minForce = np.asarray(plt.ginput(1, timeout=-1))[0][1]
        plt.close()
        
        peaks, _ = find_peaks(dat.LToes, height = fThresh, prominence= 100, distance = 200)      # prominance is min height to descend from summit to get to any higher terrain
        # x = dat.LToes
        # plt.plot(x)
        # for xc in peaks:
        #     plt.axvline(x=xc, color = 'red')

        # Toe Turns defined by the left side. calculating metrics based on 
        maxFL = []
        maxRFDupL = []
        maxRFDdnL = []
        cvPeakL = []
        noLeft = []
        maxFR = []
        maxRFDupR = []
        maxRFDdnR = []
        timeToPeakR = []
        cvPeakR = []
        noRight = []
        noToes = []
        timeToPeakL = []
        timeToPeakR = []
        totalRightToe = [] 
        totalLeftToe = []
        turnType = []
        for turn in peaks:
            tmpLForce = dat.LToes[turn-300:turn]
            tmpForceRev = np.flip(tmpLForce)
            try:
                # Defining start and end of turn #
                turnStart = turn - next(x for x, val in enumerate(tmpForceRev) if val < minForce) # first index > minForce  
                tmpForce2 = dat.LToes[turn:turn+300]
                turnEnd = turn + next(x for x, val in enumerate(tmpForce2) if val < minForce)
                # Extracting left variables #
                maxFL.append( np.max(dat.LToes[turnStart:turnEnd]) )
                maxRFDupL.append( np.max(dat.LToes[turnStart:turn].diff()) )
                maxRFDdnL.append( np.min(dat.LToes[turn:turnEnd].diff()) )
                cvPeakL.append(  (np.std(dat.LToes[turn-100 : turn+100]) / np.mean(dat.LToes[turn-100:turn+100])) * 100  )
                timeToPeakL.append( turn - turnStart )
                # Find right peak #
                rightPeak = np.argmax(dat.RToes[turnStart:turnEnd]) + turnStart
                # Extracting right variables #
                maxFR.append( np.max(dat.RToes[turnStart:turnEnd]) )
                maxRFDupR.append( np.max(dat.RToes[turnStart:rightPeak].diff()) )
                maxRFDdnR.append( np.min(dat.RToes[rightPeak:turnEnd].diff()) )
                cvPeakR.append(  (np.std(dat.RToes[rightPeak-100 : rightPeak+100]) / np.mean(dat.RToes[rightPeak-100:rightPeak+100])) * 100  )
                timeToPeakR.append( rightPeak - turnStart )
                totalRightToe.append('right')
                totalLeftToe.append('left')
                turnType.append('Toe')
            except:
                print(turn)

        leftOutcomes = pd.DataFrame({'TurnSide':list(totalLeftToe),'MaxForce':list(maxFL), 'RFDUp':list(maxRFDupL), 'RFDdn':list(maxRFDdnL),
                                     'cvPeak':list(cvPeakL), 'TimeToPk':list(timeToPeakL)})
        rightOutcomes = pd.DataFrame({'TurnSide':list(totalRightToe),'MaxForce':list(maxFR), 'RFDUp':list(maxRFDupR), 'RFDdn':list(maxRFDdnR),
                                     'cvPeak':list(cvPeakR), 'TimeToPk':list(timeToPeakR)})
        
        # meta naming #
        totalLength = len(totalRightToe) + len(totalRightToe) 
        longSubject = list( np.repeat(subName, totalLength) )
        longConfig = list( np.repeat(configName, totalLength) )
        
        totalOutput = pd.concat( [leftOutcomes, rightOutcomes] )
        totalOutput['Subject'] = longSubject
        totalOutput['Config'] = longConfig
        
        totalOutput.to_csv('C:/Users/Daniel.Feeney/Dropbox (Boa)/Snow Protocol/SnowboardProtocol/Results/snowboardResults.csv', mode='a', header=False)

    except:
        print(file)
     
   