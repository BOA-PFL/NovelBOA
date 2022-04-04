# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 16:19:30 2020
When you have both sides of data from loadsol
@author: Daniel.Feeney

Overview: using the peak force values as a starting point & looking backward 
for the force to drop below a threshold (defined with ginput) to define as the 
turn start. Looking forward from peak to the next index where force drops below 
the threshold to delimit turn completion. 

Notes from 2/21 with larger dataset:
    DF - left turns not well triggered with algorithm
    BC - some turns were too short by algo, added a continue statement to help
    RCF - only right side collected during V2 (entries 20 and 21)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import signal
from tkinter.filedialog import askopenfilenames

def findRightTurns(RForce, LForce):
    RTurns = []
    for step in range(len(RForce)-1):
        if LForce[step] <= RForce[step] and LForce[step + 1] > RForce[step + 1]:
            RTurns.append(step)
    return RTurns

def findLeftTurns(RForce, LForce):
    LTurns = []
    for step in range(len(RForce)-1):
        if RForce[step] <= LForce[step] and RForce[step + 1] > LForce[step + 1]:
            LTurns.append(step)
    return LTurns

def makeTurnPlot(RF,RL,RM,RH,LF,LL,LM,LH):
    """
    Function takes in time series of force on left and right sides
    and plots them for each turn using equal height y-axes.
    Used primarily when creating the functions but is mostly turned off for 
    processing. 
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
#fPath = 'C:\\Users\\daniel.feeney\\iCloudDrive\\iCloud~de~novel~loadsols\\SkiTesting\\'
## Example data lives in there and work well ## 
fPath = 'C:/Users/kate.harrison/Boa Technology Inc/PFL - Documents/General/Testing Segments/Snow Performance/Alpine_CuffOnSnow_Apr2022/XSENSORdata/'
entries = askopenfilenames(initialdir = fPath)

OutsideFootForce = []
OutsideFootMedialForce = []
avgOutsideHeelStart = []
propHeelLate = []
absPropHeelLate = []
cvForce = []
lTurn = []
sName = []
cName = []

for fName in entries:
    try:
        #fName = entries[2]
        dat = pd.read_csv(fName, sep = '	',skiprows = 3, header = 0, index_col = False)
        dat.columns = ['Time', 'LHeel', 'LMedial','LLateral','LTotal', 'Time2', 
                       'RLateral','RMedial','RHeel','RTotal', 'time2','accX','axxY','accZ','pass']
        
        #dat.columns = ['Time','RHeel','RLateral','RMedial','RTotal','Time2','se','ei','ni','te']
        #use above if one side only
        info = fName.split(sep = "/")[-1]
        
        subName = info.split(sep = "_")[0]
        configName = info.split(sep = "_")[1]
        
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
        
        fs = 100 
        fc = 2
        w = fc / (fs / 2)
        b, a = signal.butter(4, w, 'low')
        dat['LTotal_Filt'] = signal.filtfilt(b, a, dat.LTotal)
        dat['RTotal_Filt'] = signal.filtfilt(b, a, dat.RTotal)
        dat['RMedial_Filt'] = signal.filtfilt(b, a, dat.RMedial)
        dat['RLateral_Filt'] = signal.filtfilt(b, a, dat.RLateral)
        dat['LMedial_Filt'] = signal.filtfilt(b, a, dat.LMedial)
        dat['LLateral_Filt'] = signal.filtfilt(b, a, dat.LLateral)
        dat['LHeel_Filt'] = signal.filtfilt(b, a, dat.LHeel)
        dat['RHeel_Filt'] = signal.filtfilt(b, a, dat.RHeel)
        
        
        
        # plt.figure()
        # plt.plot(dat.RTotal_Filt, label = "Right Foot")
        # plt.plot(dat.LTotal_Filt, label = "Left Foot")
        # plt.legend()
        
        # plt.figure()
        # plt.plot(dat.LMedial_Filt, label = "Left Medial")
        # plt.plot(dat.LLateral_Filt, label = "Left Lateral")
        # plt.plot(dat.LHeel_Filt, label = "Left Heel")
        # plt.legend()
        
        plt.figure()
        #plt.plot(dat.RTotal_Filt, label = "Right Total")
        plt.plot(dat.RMedial_Filt, label = "Right Medial")
        plt.plot(dat.RLateral_Filt, label = "Right Lateral")
        plt.plot(dat.RHeel_Filt, label = "Right Heel")
        plt.legend()
        plt.title(info)
       
        RTurns = findRightTurns(dat.RTotal_Filt, dat.LTotal_Filt)
        LTurns = findLeftTurns(dat.RTotal_Filt, dat.LTotal_Filt)
        
        RTurns[:] = [x for x in RTurns if x > LTurns[0]] # we want first right turn after first left turn
        LTurns[:] = [x for x in LTurns if x < RTurns[-1]] # we want to end with a Right Turn 
        
        
        
        for i in range(len(LTurns)):
            
            #i = 0
            turnTime = RTurns[i]-LTurns[i]
            
            
            if turnTime > 100:
            
                try:
                    
                    pkIdx = np.argmax(dat.RTotal_Filt[LTurns[i]:RTurns[i]])
                    
                    ## Extract relevent parameters from a turn here ##
                    OutsideFootForce.append( dat.RTotal_Filt[LTurns[i]+pkIdx]/(dat.LTotal_Filt[LTurns[i] + pkIdx] + dat.RTotal_Filt[LTurns[i] + pkIdx]) )#proportion of force on outside foot
                    #LPeakLatProp.append( np.max(dat.LLateral_Filt[RTurns[i]:RTurns[i] + pkIdx]) ) # FOOT ROLL. Lower is better? proportion of force on outside lateral foot 
                    OutsideFootMedialForce.append( dat.RMedial_Filt[LTurns[i] + pkIdx]/dat.RTotal_Filt[LTurns[i] + pkIdx] ) # FOOT ROLL. Higher is better? proportion of force on outside medial 
                    avgOutsideHeelStart.append( np.mean(dat.RHeel_Filt[LTurns[i]:LTurns[i] + pkIdx])/np.mean(dat.RTotal_Filt[LTurns[i]:LTurns[i] + pkIdx] )) # FORWARD STANCE. Lower is better? proportion of heel force during early turn
                    
                    tmpHeel = dat["RHeel_Filt"].tolist()
                    tmpToes = dat.RMedial_Filt + dat.RLateral_Filt
                    tmpToes = tmpToes.tolist()
                    
                    pkOutsideHeelLate = np.max(dat.RHeel_Filt[LTurns[i] + pkIdx: RTurns[i]])
                    pkOutsideHeelLateIdx = tmpHeel.index(pkOutsideHeelLate) # BALANCE? Peak heel force late in turn. Weight should shift from toe towards heel late in turn - closer to 50% is better? 
                    propHeelLate.append(( pkOutsideHeelLate - tmpToes[pkOutsideHeelLateIdx])/dat.RTotal_Filt[pkOutsideHeelLateIdx])
                    absPropHeelLate.append(abs(propHeelLate[-1]))
                    lTurn.append('Left') 
                    sName.append(subName)
                    cName.append(configName)
                except:
                    print(fName + str(i))
            
    except:
        print(fName + "WHOLE TRIAL SUCKS")

outcomes = pd.DataFrame({'Subject':list(sName),'Config':list(cName),'TurnType': list(lTurn),
                                 'OutsideFootForce':list(OutsideFootForce), 'OutsideFootMedialForce':list(OutsideFootMedialForce),'avgOutsideHeelStart':list(avgOutsideHeelStart),
                                 'propHeelLate':list(propHeelLate), 'absPropHeelLate':list(absPropHeelLate)
                                 })
         
        
outfileName = fPath + 'CompiledResults2.csv'

if os.path.exists(outfileName) == False:
    
    outcomes.to_csv(outfileName, mode='a', header=True, index = False)

else:
    outcomes.to_csv(outfileName, mode='a', header=False, index = False)        
    #     ## Right Turns with greater force on left side ##
    #     leftPeaks, _ = signal.find_peaks(dat.LTotal, height = fThresh, prominence= 100, distance = 100)      # prominance is min height to descend from summit to get to any higher terrain
    #     # #Visualization to check peaks
    #     # x = dat['LTotal']
    #     # plt.plot(x)
    #     # for xc in leftPeaks:
    #     #     plt.axvline(x=xc, color = 'red')
        
    #     rPks = []
    #     rPkLatEarly = []
    #     rPkMedEarly = []
    #     ravgHeelStart = []
    #     rpkHeelLate = []
    #     rpkMedInsideLate = []
    #     rcvForce = []
    #     rTurn = []
    #     sName2 = []
    #     cName2 = []
    #     for turn in leftPeaks[1:len(leftPeaks)]:
    #         # subset the force signal, flip np array and find first index below threshold (turnStart)
    #         tmpForce = np.array(dat.LTotal[turn-80:turn])
    #         tmpForceRev = np.flip(tmpForce)
    #         try:
    #             turnStart = turn - next(x for x, val in enumerate(tmpForceRev) if val < minForce) # first index > minForce  
                 
    #             tmpForce2 = dat.LTotal[turn:turn+100]
    #             turnEnd = turn + next(x for x, val in enumerate(tmpForce2) if val < minForce)
    #             turnLen = turnEnd - turnStart
    #             if turnLen < 50:
    #                 continue
                
    #             ts = 0
    #             te = turnEnd - turnStart
    #             tp = turn - turnStart
                            
    #             tmpRF = np.array( dat.RTotal[turnStart:turnEnd] )
    #             tmpRL = np.array( dat.RLateral[turnStart:turnEnd] )
    #             tmpRM = np.array( dat.RMedial[turnStart:turnEnd] )
    #             tmpRH = np.array( dat.RHeel[turnStart:turnEnd] )
    #             tmpLF = np.array( dat.LTotal[turnStart:turnEnd] )
    #             tmpLL = np.array( dat.LLateral[turnStart:turnEnd] )
    #             tmpLM = np.array( dat.LMedial[turnStart:turnEnd] )
    #             tmpLH = np.array( dat.LHeel[turnStart:turnEnd] )
    #             #makeTurnPlot(tmpRF, tmpRL, tmpRM, tmpRH, tmpLF, tmpLL, tmpLM, tmpLH)
                
    #             rPks.append( np.max(tmpLF) ) #peak force
    #             rPkLatEarly.append( np.max(tmpLL[ts:tp]) )
    #             rPkMedEarly.append( np.max(tmpLM[ts:tp]) )
    #             ravgHeelStart.append( np.mean(tmpLH[ts:tp]) )
    #             rpkHeelLate.append( np.max(tmpLH[tp:te]) )
    #             rpkMedInsideLate.append( np.max(tmpRM[tp:te]) )
    #             rcvForce.append( (np.std(tmpLF[tp-20:tp+20]) / np.mean(tmpLF[tp-20:tp+20])) * 100 )
    #             rTurn.append('Right')
    #             sName2.append(subName)
    #             cName2.append(configName)
    #         except:
    #             print('never reached min force turn ' + str(turn))
    #     ### end of turn detection algorithm ###
        
    #     outcomesRight = pd.DataFrame({'Subject':list(sName2),'Config':list(cName2),'TurnType': list(rTurn),'PeakForce':list(rPks),
    #                      'PkLatForceEarly':list(rPkLatEarly), 'PkMedForceEarly':list(rPkMedEarly),
    #                      'PkHeelLate':list(rpkHeelLate),'pkMedInsideLate':list(rpkMedInsideLate),
    #                      'CVForce':list(rcvForce)})

    #     totalOutput = pd.concat( [outcomes, outcomesRight] )
    #     totalOutput.to_csv('C:\\Users\\daniel.feeney\\Boa Technology Inc\\PFL - General\\Snow Performance\\Alpine_V1vsV2_Internal_Feb2022\\testresults.csv', mode='a', header=False)

    # except:
    #     print(file)