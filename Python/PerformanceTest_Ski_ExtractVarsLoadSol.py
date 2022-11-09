# -*- coding: utf-8 -*-
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


def findRightTurns(RForce, LForce):
    """
    Find start of right turns. Defined as when force on LEFT ski (i.e. downhill ski) 
    exceeds force on right ski.

    Parameters
    ----------
    RForce : numpy array
        Time series of force data under right foot.
    LForce : numpy array
        Time series of force data under left foot. 

    Returns
    -------
    RTurns : numpy array
        Index of frames where right turns started. 

    """
    RTurns = []

    for step in range(len(RForce)-1):
        if LForce[step] <= RForce[step] and LForce[step + 1] > RForce[step + 1]:
            RTurns.append(step)
    return RTurns

def findLeftTurns(RForce, LForce):
    """
    Find start of left turns. Defined as when force on RIGHT ski 
    (i.e. downhill ski) exceeds for on the left ski.

    Parameters
    ----------
    RForce : numpy array
        Time series of force data under the right foot. 
    LForce : numpy array
        Time series of force data under the left foot. 

    Returns
    -------
    LTurns : numpy array
        Index of frames where left turns started. 

    """
    LTurns = []
    for step in range(len(RForce)-1):
        if RForce[step] <= LForce[step] and RForce[step + 1] > LForce[step + 1]:
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

#makeTurnPlot(dat, RTurns, 'Right')
## Example data lives in there and work well ## 

fPath = 'C:/Users/daniel.feeney/Boa Technology Inc/PFL - Documents/General/Testing Segments/Snow Performance/Alpine_CuffOnSnow_Apr2022/XSENSORdata/'
entries = askopenfilenames(initialdir = fPath)

# Initiate discrete outcome variables
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
        # Loop through files and use time series force data to identify turns
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
        
       
        RTurns = findRightTurns(dat.RTotal_Filt, dat.LTotal_Filt)
        LTurns = findLeftTurns(dat.RTotal_Filt, dat.LTotal_Filt)
        
        RTurns[:] = [x for x in RTurns if x > LTurns[0]] # we want first right turn after first left turn
        LTurns[:] = [x for x in LTurns if x < RTurns[-1]] # we want to end with a Right Turn 
        
        RTurns = cleanTurns(RTurns)
        LTurns = cleanTurns(LTurns)
        
        makeTurnPlot(dat, RTurns, 'Right Turns')
        makeTurnPlot(dat, LTurns, 'Left Turns')
    

        for i, value in enumerate(LTurns):
            # Loop through all cleaned turns to calculate discrete outcome measures.
            # using right side only (left turns). No longer assuming left to right
            # transitions are constant and only using data from one side. 
            # Could port this to do R Turns as well. 

            try:
                
                pkIdx = np.argmax(dat.RTotal_Filt[LTurns[i]:LTurns[i+1]])
                pkIdx = value + pkIdx
                
                ## Extract relevent parameters from a turn here ##
                # EARLY TURN DH FORCE: Proportion of force on downhill foot. Higher is better
                OutsideFootForce.append( dat.RTotal_Filt[pkIdx]/(dat.LTotal_Filt[pkIdx] + dat.RTotal_Filt[pkIdx]) )
                # FOOT ROLL. Proportion of force on outside medial 
                OutsideFootMedialForce.append( dat.RMedial_Filt[pkIdx]/dat.RTotal_Filt[pkIdx] ) 
                # FORWARD STANCE. Proportion of heel force during early turn
                avgOutsideHeelStart.append( np.mean(dat.RHeel_Filt[value:pkIdx])/np.mean(dat.RTotal_Filt[value:pkIdx] )) 
                
                tmpHeel = dat["RHeel_Filt"].tolist()
                tmpToes = dat.RMedial_Filt + dat.RLateral_Filt
                tmpToes = tmpToes.tolist()
                
                pkOutsideHeelLate = np.max(dat.RHeel_Filt[pkIdx: LTurns[i+1]])
                pkOutsideHeelLateIdx = tmpHeel.index(pkOutsideHeelLate) 
                propHeelLate.append(( pkOutsideHeelLate - tmpToes[pkOutsideHeelLateIdx])/dat.RTotal_Filt[pkOutsideHeelLateIdx])
                # BALANCE - proportion of force on heel vs. toes late in turn (50% is ideal)
                absPropHeelLate.append(abs(propHeelLate[-1])) 
                lTurn.append('Left') 
                sName.append(subName)
                cName.append(configName)
                
            except:
                print(fName + str(i))
        
    except:
        print(fName)


# Create data frame with all outcome measures and export to csv. Will create a new csv if one does not exist for this dataset. 
# Otherwise will append.

outcomes = pd.DataFrame({'Subject':list(sName),'Config':list(cName),'TurnType': list(lTurn),
                                 'OutsideFootForce':list(OutsideFootForce), 'OutsideFootMedialForce':list(OutsideFootMedialForce),'avgOutsideHeelStart':list(avgOutsideHeelStart),
                                 'propHeelLate':list(propHeelLate), 'absPropHeelLate':list(absPropHeelLate)
                                 })
         
        
# outfileName = fPath + 'CompiledResults2.csv'

# if os.path.exists(outfileName) == False:
    
#     outcomes.to_csv(outfileName, mode='a', header=True, index = False)

# else:
#     outcomes.to_csv(outfileName, mode='a', header=False, index = False) 


