# -*- coding: utf-8 -*-3
"""
Created on Wed Dec 16 16:19:30 2020
When you have both sides of data from loadsol

Note: 'left turns' means skier is turning to the left but more force should
be on the downhill (or right) ski. 
'right turns' means skier is turning to the right but more force should be on 
the downhill (or left) ski.  - DF


"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import signal
from tkinter.filedialog import askopenfilenames
from tkinter import messagebox


fPath = 'C:\\Users\\daniel.feeney\\Boa Technology Inc\\PFL Team - General\\Testing Segments\\Snow Performance\\\SkiValidation_Dec2022\Loadsol\\'
fileExt = r".txt"
entries = [fName for fName in os.listdir(fPath) if fName.endswith(fileExt)]
#entries = askopenfilenames(initialdir = fPath)



def findRightTurns(RForce, LForce):
    """
    Find start of right turns. Defined as when force on LEFT ski (i.e. downhill ski) 
    exceeds force on right ski.

    Parameters
    ----------
    RForce : Pandas Column
        Time series of force data under right foot.
    LForce : Pandas Column
        Time series of force data under left foot. 

    Returns
    -------
    RTurns : List
        Index of frames where right turns started. 

    """
    RTurns = []
    for step in range(len(RForce)-1):
        if LForce[step] <= RForce[step] and LForce[step + 1] > RForce[step + 1] and np.mean(LForce[step:step+200]) > np.mean(RForce[step:step+200]):
            RTurns.append(step)
    return RTurns

def findLeftTurns(RForce, LForce):
    """
    Find start of left turns. Defined as when force on RIGHT ski 
    (i.e. downhill ski) exceeds for on the left ski.

    Parameters
    ----------
    RForce : Pandas Column
        Time series of force data under the right foot. 
    LForce : Pandas Column
        Time series of force data under the left foot. 

    Returns
    -------
    LTurns : list
        Index of frames where left turns started. 

    """
    LTurns = []
    for step in range(len(RForce)-1):
        if RForce[step] <= LForce[step] and RForce[step + 1] > LForce[step + 1] and np.mean(RForce[step:step+200]) > np.mean(LForce[step:step+200]) :
            LTurns.append(step)
    return LTurns



def cleanTurns(turnIndices):
    """
    Finds if there are two turn initiations from the same foot close together
    and removes the first one
    """
    lowThreshold = np.median(np.diff(turnIndices)) * 0.5
    trueIndices = [idx for idx, x in enumerate(list(np.diff(turnIndices))) if x < lowThreshold]
    
    cleanList = [x for idx, x in enumerate(turnIndices) if idx not in trueIndices]
    
    return(cleanList)


def makeTurnPlot(inputDF, turnIndices, turnSide):
    """
    Parameters
    ----------
    inputDF : Pandas df
        Loadsol Export with low-pass filtered force data
    turnIndices : List
        Indices of turns.
    turnSide : String
        Left or Right.

    Returns
    -------
    Plot

    """
    plt.figure()
    plt.plot(inputDF.RTotal_Filt, label = "Right Foot")
    plt.plot(inputDF.LTotal_Filt, label = "Left Foot")
    plt.vlines(x = turnIndices, ymin = 0, ymax = np.max(600),
              color = 'k', label = turnSide, linewidth=3.0, ls='--')
    plt.legend()
    plt.title(turnSide)
    plt.suptitle(subName)


# Initiate discrete outcome variables
OutsideFootProp = []
OutsideFootForce = []
OutsideFootMedialProp = []
OutsideFootMedialForce = []
avgOutsideHeelStartForce = []
avgOutsideHeelStartProp = []
propHeelLate = []
absPropHeelLate = []
RFD = []
cvForce = []
turnSide = []
sName = []
cName = []
timeToPeak = []
badFileList = []

for fName in entries:
    try:
        # Loop through files and use time series force data to identify turns
        #fName = entries[2]
        dat = pd.read_csv(fPath+fName, sep = '	',skiprows = 3, header = 0, index_col = False)
        dat.columns = ['Time', 'LHeel', 'LMedial','LLateral','LTotal', 'Time2', 
                       'RLateral','RMedial','RHeel','RTotal', 'time2','accX','axxY','accZ','pass']
        
        #dat.columns = ['Time','RHeel','RLateral','RMedial','RTotal','Time2','se','ei','ni','te']
        #use above if one side only
        info = fName.split(sep = "/")[-1]
        
        subName = info.split(sep = "_")[0]
        configName = info.split(sep = "_")[1]
        
        dat['LToes'] = dat.LMedial + dat.LLateral
        dat['RToes'] = dat.RMedial + dat.RLateral
        
        
        
        # Load in the trial segmentation variable if it is in the directory
        if os.path.exists(fPath+fName+'TrialSeg.npy') == True:
            trial_segment_old = np.load(fPath+fName+'TrialSeg.npy',allow_pickle=True)
            trialStart = trial_segment_old[1][0,0]
            trialEnd = trial_segment_old[1][1,0]
            dat = dat.iloc[int(np.floor(trialStart)) : int(np.floor(trialEnd)),:]
            dat = dat.reset_index()
            
        else:
                
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
        fc = 6
        w = fc / (fs / 2)
        b, a = signal.butter(2, w, 'low')
        dat['LTotal_Filt'] = signal.filtfilt(b, a, dat.LTotal)
        dat['RTotal_Filt'] = signal.filtfilt(b, a, dat.RTotal)
        dat['RMedial_Filt'] = signal.filtfilt(b, a, dat.RMedial)
        dat['RLateral_Filt'] = signal.filtfilt(b, a, dat.RLateral)
        dat['LMedial_Filt'] = signal.filtfilt(b, a, dat.LMedial)
        dat['LLateral_Filt'] = signal.filtfilt(b, a, dat.LLateral)
        dat['LHeel_Filt'] = signal.filtfilt(b, a, dat.LHeel)
        dat['RHeel_Filt'] = signal.filtfilt(b, a, dat.RHeel)
        
       
        RTurns = findRightTurns(dat.RTotal_Filt, dat.LTotal_Filt)
        LTurns = findLeftTurns(dat.RTotal_Filt, dat.LTotal_Filt)

        RTurns = cleanTurns(RTurns)
        LTurns = cleanTurns(LTurns)
        
        makeTurnPlot(dat, LTurns, 'Left Turns')
        makeTurnPlot(dat, RTurns, 'Right Turns')
        answer = messagebox.askyesno("Question","Is data clean?")
        
        if answer == False:
            plt.close('all')
            print('Adding file to bad file list')
            badFileList.append(fName)
        
        if answer == True:
            plt.close('all')
            print('Estimating point estimates')
        
        
            for i, value in enumerate(LTurns):
                # Loop through all cleaned turns to calculate discrete outcome measures.
                # using right side only (left turns). No longer assuming left to right
                # transitions are constant and only using data from one side. 
                # Extend this to Left Turns as well
                # variables of interest: outside force (downhill foot peak force) higher is better
                # Outside medial peak force (higher is better)
                # Outside heel force lower is better with an average higher force on forefoot Tricky**
    
                try:
                    
                    pkIdx = np.argmax(dat.RTotal_Filt[LTurns[i]:LTurns[i+1]])
                    pkIdx = value + pkIdx
                    timeToPeak.append(pkIdx - value)
                    RFD.append((dat.RTotal_Filt[pkIdx] - dat.RTotal_Filt[value]) / (pkIdx - value))
                    
                    ## Extract relevent parameters from a turn here ##
                    # EARLY TURN DH FORCE: Proportion of force on downhill foot. Higher is better
                    OutsideFootProp.append( dat.RTotal_Filt[pkIdx]/(dat.LTotal_Filt[pkIdx] + dat.RTotal_Filt[pkIdx]) )
                    OutsideFootForce.append( dat.RTotal_Filt[pkIdx] )
                    # FOOT ROLL. Proportion or total of force on outside medial 
                    OutsideFootMedialProp.append( dat.RMedial_Filt[pkIdx]/dat.RTotal_Filt[pkIdx] ) 
                    OutsideFootMedialForce.append( dat.RMedial_Filt[pkIdx] )
                    # FORWARD STANCE. Proportion of heel force during early turn
                    avgOutsideHeelStartProp.append( np.mean(dat.RHeel_Filt[value:pkIdx])/np.mean(dat.RTotal_Filt[value:pkIdx] )) 
                    avgOutsideHeelStartForce.append( np.mean(dat.RHeel_Filt[value:pkIdx]) )
                                        
                    turnSide.append('Left') 
                    sName.append(subName)
                    cName.append(configName)
                    # tmpHeel = dat["RHeel_Filt"].tolist()
                    # tmpToes = dat.RMedial_Filt + dat.RLateral_Filt
                    # tmpToes = tmpToes.tolist()
                    
                    # #More work needed to figure out heel force before we can analyze
                    # pkOutsideHeelLate = np.max(dat.RHeel_Filt[pkIdx: LTurns[i+1]])
                    # pkOutsideHeelLateIdx = tmpHeel.index(pkOutsideHeelLate) 
                    # propHeelLate.append(( pkOutsideHeelLate - tmpToes[pkOutsideHeelLateIdx])/dat.RTotal_Filt[pkOutsideHeelLateIdx])
                    # # BALANCE - proportion of force on heel vs. toes late in turn (50% is target)
                    # absPropHeelLate.append(abs(propHeelLate[-1])) 

                    
                except Exception as e: print(e)
                    
            for i, value in enumerate(RTurns):
                # Loop through all cleaned Right Turns to calculate discrete outcome measures.
                # using right side only (left turns). No longer assuming left to right
                # transitions are constant and only using data from one side. 
                # Extend this to Left Turns as well
                # variables of interest: outside force (downhill foot peak force) higher is better
                # Outside medial peak force (higher is better)
                # Outside heel force lower is better with an average higher force on forefoot Tricky**
    
                try:
                    
                    pkIdx = np.argmax(dat.LTotal_Filt[LTurns[i]:LTurns[i+1]])
                    pkIdx = value + pkIdx
                    timeToPeak.append(pkIdx - value)
                    RFD.append((dat.LTotal_Filt[pkIdx] - dat.LTotal_Filt[value]) / (pkIdx - value))

                    ## Extract relevent parameters from a turn here ##
                    # EARLY TURN DH FORCE: Proportion of force on downhill foot. Higher is better
                    OutsideFootProp.append( dat.LTotal_Filt[pkIdx]/(dat.LTotal_Filt[pkIdx] + dat.RTotal_Filt[pkIdx]) )
                    OutsideFootForce.append( dat.LTotal_Filt[pkIdx] )
                    # FOOT ROLL. Proportion or total of force on outside medial 
                    OutsideFootMedialProp.append( dat.LMedial_Filt[pkIdx]/dat.LTotal_Filt[pkIdx] ) 
                    OutsideFootMedialForce.append( dat.LMedial_Filt[pkIdx] )
                    # FORWARD STANCE. Proportion of heel force during early turn
                    avgOutsideHeelStartProp.append( np.mean(dat.LHeel_Filt[value:pkIdx])/np.mean(dat.LTotal_Filt[value:pkIdx] )) 
                    # heel force should be lower at initial start
                    avgOutsideHeelStartForce.append( np.mean(dat.LHeel_Filt[value:pkIdx]) )
                    # heel force should be higher by end
                    #avgOutsideHeelEndForce.append( np.mean(dat.LHeel_Filt[pkIdx:pkIdx+50]) )
                    turnSide.append('Right') 
                    sName.append(subName)
                    cName.append(configName)
                    
                    tmpHeel = dat["LHeel_Filt"].tolist()
                    tmpToes = dat.LMedial_Filt + dat.LLateral_Filt
                    tmpToes = tmpToes.tolist()
                    
                    #More work needed to figure out heel force before we can analyze
                    # pkOutsideHeelLate = np.max(dat.RHeel_Filt[pkIdx: LTurns[i+1]])
                    # pkOutsideHeelLateIdx = tmpHeel.index(pkOutsideHeelLate) 
                    # propHeelLate.append(( pkOutsideHeelLate - tmpToes[pkOutsideHeelLateIdx])/dat.RTotal_Filt[pkOutsideHeelLateIdx])
                    # # BALANCE - proportion of force on heel vs. toes late in turn (50% is ideal)
                    # absPropHeelLate.append(abs(propHeelLate[-1])) 

                    
                except Exception as e: print(e)
                        
    except:
        print(fName)


# Create data frame with all outcome measures and export to csv. Will create a new csv if one does not exist for this dataset. 
# Otherwise will append.

outcomes = pd.DataFrame({'Subject':list(sName),'Config':list(cName),'TurnType': list(turnSide),
                                 'OutsideFootForce':list(OutsideFootForce), 'OutsideFootMedialForce':list(OutsideFootMedialForce),'OutsideFootProp':list(OutsideFootProp),
                                 'avgOutsideHeelStartForce':list(avgOutsideHeelStartForce),'OutsideFootMedialProp':list(OutsideFootMedialProp)
                                 })
         
  
outfileName = fPath + 'CompiledResultsTest.csv'

if os.path.exists(outfileName) == False:
    
    outcomes.to_csv(outfileName, mode='a', header=True, index = False)

else:
    outcomes.to_csv(outfileName, mode='a', header=False, index = False) 


