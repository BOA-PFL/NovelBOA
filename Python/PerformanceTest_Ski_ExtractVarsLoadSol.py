# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 16:19:30 2020
When you have both sides of data from loadsol
@author: Daniel.Feeney
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.signal import find_peaks

def makeTurnPlot(RF,RL,RM,RH,LF,LL,LM,LH):
    """
    Function takes in time series of force on left and right sides
    and plots them for each turn using equal height y-axes
    """
    fig, (ax1, ax2) = plt.subplots(1,2)
    ax1.plot(RF, color = 'r', label = 'Right Force')
    ax1.plot(RL, color = 'k', label = 'Right Lateral')
    ax1.plot(RM, color = 'b', label = 'Right Medial')
    ax1.plot(RH, color = 'g', label = 'Right Heel')
    ax1.set_ylim([0,800])
    ax1.legend()
    ax2.plot(LF, color = 'r', label = 'Left Force')
    ax2.plot(LL, color = 'k', label = 'Left Lateral')
    ax2.plot(LM, color = 'b', label = 'Left Medial')
    ax2.plot(LH, color = 'g', label = 'Left Heel')
    ax2.set_ylim([0,800])
    ax2.legend()
    
# Read in files
# only read .asc files for this work
fPath = 'C:\\Users\\daniel.feeney\\iCloudDrive\\iCloud~de~novel~loadsols\\SkiTesting\\'
entries = os.listdir(fPath)


for file in entries:
    try:
        fName = file
        dat = pd.read_csv(fPath+fName, sep = '	',skiprows = 3, header = 0, index_col = False)
        dat.columns = ['Time', 'LHeel', 'LMedial','LLateral','LTotal', 'Time2', 
                       'RLateral','RMedial','RHeel','RTotal', 'time2','accX','axxY','accZ','pass']
        #subName = fName.split(sep = "_")[0]
        #configName = fName.split(sep = "_")[1]
        
        dat['LToes'] = dat.LMedial + dat.LLateral
        dat['RToes'] = dat.RMedial + dat.RLateral
        
        ### Subset the trial to a portion that does not include standing ###
        fig, ax = plt.subplots()
        ax.plot(dat.LTotal, label = 'Left Total Force')
        ax.plot(dat.RTotal, label = 'Right Total Force')
        fig.legend()
        print('Select start and end of analysis trial')
        pts = np.asarray(plt.ginput(2, timeout=-1))
        plt.close()
        # downselect the region of the dataframe you selected from above 
        dat = dat.iloc[int(np.floor(pts[0,0])) : int(np.floor(pts[1,0])),:]
        dat = dat.reset_index()

        #### Working on peak and turn initiation detection algorithm here ####
        ### This will work on the right side initially and abstract to the left ## 
        plt.plot(dat.RTotal)
        print('Select threshold force where peak forces will be above')
        fThresh = np.asarray(plt.ginput(1, timeout=-1))[0][1]
        plt.close()
        # select minimal force for turn initiation #
        plt.plot(dat.RTotal)
        print('Select threshold force where min forces will be BELOW')
        minForce = np.asarray(plt.ginput(1, timeout=-1))[0][1]
        plt.close()
        
        peaks, _ = find_peaks(dat.RTotal, height = fThresh, prominence= 100, distance = 100)      # prominance is min height to descend from summit to get to any higher terrain
        # #Visualization to check peaks
        # x = dat['RTotal']
        # plt.plot(x)
        # for xc in peaks:
        #     plt.axvline(x=xc, color = 'red')
        
        ## Left turns with greatest force on right side (right is outside ski) ## 
        lPks = []
        lPkLatEarly = []
        lPkMedEarly = []
        avgHeelStart = []
        pkHeelLate = []
        pkMedInsideLate = []
        cvForce = []
        lTurn = []
        for turn in peaks[1:len(peaks)]:
            # subset the force signal, flip np array and find first index below threshold (turnStart)
            tmpForce = np.array(dat.RTotal[turn-80:turn])
            tmpForceRev = np.flip(tmpForce)
            try:
                turnStart = turn - next(x for x, val in enumerate(tmpForceRev) if val < minForce) # first index > minForce  
                tmpForce2 = dat.RTotal[turn:turn+100]
                turnEnd = turn + next(x for x, val in enumerate(tmpForce2) if val < minForce)
                ### CHECK THESE TODO ##
                ts = 0
                te = turnEnd - turnStart
                tp = turn - turnStart
                
                tmpRF = np.array( dat.RTotal[turnStart:turnEnd] )
                tmpRL = np.array( dat.RLateral[turnStart:turnEnd] )
                tmpRM = np.array( dat.RMedial[turnStart:turnEnd] )
                tmpRH = np.array( dat.RHeel[turnStart:turnEnd] )
                tmpLF = np.array( dat.LTotal[turnStart:turnEnd] )
                tmpLL = np.array( dat.LLateral[turnStart:turnEnd] )
                tmpLM = np.array( dat.LMedial[turnStart:turnEnd] )
                tmpLH = np.array( dat.LHeel[turnStart:turnEnd]  )
                #makeTurnPlot(tmpRF, tmpRL, tmpRM, tmpRH, tmpLF, tmpLL, tmpLM, tmpLH)
                
                ## Extract relevent parameters from a turn here ##
                lPks.append( np.max(tmpRF) ) #peak force
                lPkLatEarly.append( np.max(tmpRL[ts:tp]) )
                lPkMedEarly.append( np.max(tmpRM[ts:tp]) )
                avgHeelStart.append( np.mean(tmpRH[ts:tp]) )
                pkHeelLate.append( np.max(tmpRH[tp:te]) )
                pkMedInsideLate.append( np.max(tmpLM[tp:te]) )
                cvForce.append( (np.std(tmpRF[tp-20:tp+20]) / np.mean(tmpRF[tp-20:tp+20])) * 100 )
                lTurn.append('left')
            except:
                print('never reached min force turn ' + str(turn))
            
        outcomes = pd.DataFrame({'TurnType': list(lTurn),'PeakForce':list(lPks),
                                 'PkLatForceEarly':list(lPkLatEarly), 'PkMedForceEarly':list(lPkMedEarly),
                                 'PkHeelLate':list(pkHeelLate),'pkMedInsideLate':list(pkMedInsideLate),
                                 'CVForce':list(cvForce)})
         
        
        
        ## Right Turns with greater force on left side ##
        leftPeaks, _ = find_peaks(dat.LTotal, height = fThresh, prominence= 100, distance = 100)      # prominance is min height to descend from summit to get to any higher terrain
        # #Visualization to check peaks
        # x = dat['LTotal']
        # plt.plot(x)
        # for xc in leftPeaks:
        #     plt.axvline(x=xc, color = 'red')
        
        rPks = []
        rPkLatEarly = []
        rPkMedEarly = []
        ravgHeelStart = []
        rpkHeelLate = []
        rpkMedInsideLate = []
        rcvForce = []
        rTurn = []
        for turn in leftPeaks[1:len(leftPeaks)]:
            # subset the force signal, flip np array and find first index below threshold (turnStart)
            tmpForce = np.array(dat.LTotal[turn-80:turn])
            tmpForceRev = np.flip(tmpForce)
            try:
                turnStart = turn - next(x for x, val in enumerate(tmpForceRev) if val < minForce) # first index > minForce  
                 
                tmpForce2 = dat.LTotal[turn:turn+100]
                turnEnd = turn + next(x for x, val in enumerate(tmpForce2) if val < minForce)
                ts = 0
                te = turnEnd - turnStart
                tp = turn - turnStart
                            
                tmpRF = np.array( dat.RTotal[turnStart:turnEnd] )
                tmpRL = np.array( dat.RLateral[turnStart:turnEnd] )
                tmpRM = np.array( dat.RMedial[turnStart:turnEnd] )
                tmpRH = np.array( dat.RHeel[turnStart:turnEnd] )
                tmpLF = np.array( dat.LTotal[turnStart:turnEnd] )
                tmpLL = np.array( dat.LLateral[turnStart:turnEnd] )
                tmpLM = np.array( dat.LMedial[turnStart:turnEnd] )
                tmpLH = np.array( dat.LHeel[turnStart:turnEnd] )
                makeTurnPlot(tmpRF, tmpRL, tmpRM, tmpRH, tmpLF, tmpLL, tmpLM, tmpLH)
                
                rPks.append( np.max(tmpLF) ) #peak force
                rPkLatEarly.append( np.max(tmpLL[ts:tp]) )
                rPkMedEarly.append( np.max(tmpLM[ts:tp]) )
                ravgHeelStart.append( np.mean(tmpLH[ts:tp]) )
                rpkHeelLate.append( np.max(tmpLH[tp:te]) )
                rpkMedInsideLate.append( np.max(tmpRM[tp:te]) )
                rcvForce.append( (np.std(tmpLF[tp-20:tp+20]) / np.mean(tmpLF[tp-20:tp+20])) * 100 )
                rTurn.append('Right')
            except:
                print('never reached min force turn ' + str(turn))
        ### end of turn detection algorithm ###
        
        outcomesRight = pd.DataFrame({'TurnType': list(rTurn),'PeakForce':list(rPks),
                         'PkLatForceEarly':list(rPkLatEarly), 'PkMedForceEarly':list(rPkMedEarly),
                         'PkHeelLate':list(rpkHeelLate),'pkMedInsideLate':list(rpkMedInsideLate),
                         'CVForce':list(rcvForce)})

        totalOutput = pd.concat( [outcomes, outcomesRight] )
        totalOutput.to_csv('C:\\Users\\daniel.feeney\\Boa Technology Inc\\PFL - General\\Snow_Alpine_Pilot\\LoadsolAlpineTesting\\testresults.csv', mode='a', header=False)

    except:
        print(file)
     
  
