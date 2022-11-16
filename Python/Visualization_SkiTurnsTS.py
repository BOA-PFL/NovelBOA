# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 15:07:14 2022

Code to generate time series plots from two models for ski turns

@author: Dan.Feeney
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import os
from tkinter import messagebox
from dataclasses import dataclass

### main inputs ###

# select files

fPath = 'C:\\Users\\daniel.feeney\\Boa Technology Inc\\PFL Team - General\\Testing Segments\\Snow Performance\\Alpine_CuffOnSnow_Apr2022\\XSENSORdata\\'
fileExt = r".txt"
entries = [fName for fName in os.listdir(fPath) if fName.endswith(fileExt)]
entries = os.listdir(fPath)

# Choose two files to compare
fName1 = entries[2]
fName2 = entries[3]

stepLen = 400

### set plot font size ###
SMALL_SIZE = 14
MEDIUM_SIZE = 16
BIGGER_SIZE = 18

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

pd.options.mode.chained_assignment = None  # default='warn' set to warn for a lot of warnings
save_on = 0
TS_plots = 0 # turn on for time series plotting but not encouraged for running thru all files

# instantiate a dataclass to be used below
@dataclass
class TurnStarts:
    name: str
    config: str
    leftstarts: list
    rightstarts: list
    df: pd.DataFrame

@dataclass
class EnsembleData:
    meanData: np.ndarray
    sdData: np.ndarray
    
def calcEnsembleData(Force, TurnStarts, stepLength):
    """ 
    Calculates ensemble averaged force and SD data from a time series
    """
    rightForceMat = forceMatrix(Force, TurnStarts, len(TurnStarts), stepLength)
    meanRightForce = np.mean(rightForceMat, axis = 0)
    sdRightForce = np.std(rightForceMat, axis = 0)
    
    result = EnsembleData(meanRightForce, sdRightForce)
    
    return(result)
    

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
        if LForce[step] <= RForce[step] and LForce[step + 1] > RForce[step + 1] and np.mean(LForce[step:step+200] > np.mean(RForce[step:step+200])):
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

# preallocate matrix for force and fill in with force data
def forceMatrix(inputForce, landings, noSteps, stepLength):
    """
    input a force signal, return matrix with n rows (for each landing) by m col
    #for each point in stepLen
    """
    preForce = np.zeros((noSteps,stepLength))
    
    for iterVar, landing in enumerate(landings):
        try:
            preForce[iterVar,] = inputForce[landing:landing+stepLength]
        except:
            print(landing)
            
    return preForce
                
def calcTurnStarts(fName):
    """
    
    Parameters
    ----------
    fName : str
        the filename to a relevant txt file with loadsol data.

    Returns an instance of dataclass with config, subname, turn start indices
    that will be used in subsequent plotting
    -------
    
    """
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

    RTurns = cleanTurns(RTurns)
    LTurns = cleanTurns(LTurns)
    
    makeTurnPlot(dat, LTurns, 'Left Turns')
    answer = messagebox.askyesno("Question","Is data clean?")
    
    if answer == False:
        plt.close('all')
        print('Adding file to bad file list')
        #badFileList.append(fName)
    
    if answer == True:
        plt.close('all')
        print('Estimating point estimates')

   
    result = TurnStarts(subName, configName, LTurns, RTurns, df = dat)
               
    return(result)

turnsRound1 = calcTurnStarts(fName1)
turnsRound2 = calcTurnStarts(fName2)
x = np.linspace(0,stepLen,stepLen)

## Calculating averaged data ## 
rightForceLTurn = calcEnsembleData(turnsRound1.df.RTotal_Filt, turnsRound1.leftstarts,stepLen)
leftForceRTurn = calcEnsembleData(turnsRound1.df.LTotal_Filt, turnsRound1.rightstarts, stepLen)
rightForceLTurn2 = calcEnsembleData(turnsRound2.df.RTotal_Filt, turnsRound2.leftstarts,stepLen)
leftForceRTurn2 = calcEnsembleData(turnsRound2.df.LTotal_Filt, turnsRound2.rightstarts, stepLen)

### plotting ### 
fig, (ax1, ax2) = plt.subplots(2)
ax1.plot(x, rightForceLTurn.meanData, 'k', color='#DC582A')
ax1.fill_between(x,rightForceLTurn.meanData-rightForceLTurn.sdData, rightForceLTurn.meanData+rightForceLTurn.sdData,
    alpha=0.5, edgecolor='#DC582A', facecolor='#DC582A', label = turnsRound1.config)
ax1.plot(x, rightForceLTurn2.meanData, 'k', color='#00966C')
ax1.fill_between(x,rightForceLTurn2.meanData-rightForceLTurn2.sdData, rightForceLTurn2.meanData+rightForceLTurn2.sdData,
    alpha=0.5, edgecolor='#00966C', facecolor='#00966C', label = turnsRound2.config)
ax1.legend()
ax1.set_title('Turns to the Left')
ax1.set_ylabel('Right Force (N)')

ax2.plot(x, leftForceRTurn.meanData, 'k', color='#DC582A')
ax2.fill_between(x,leftForceRTurn.meanData-leftForceRTurn.sdData, leftForceRTurn.meanData+leftForceRTurn.sdData,
    alpha=0.5, edgecolor='#DC582A', facecolor='#DC582A', label = turnsRound1.config)
ax2.plot(x, leftForceRTurn2.meanData, 'k', color='#00966C')
ax2.fill_between(x,leftForceRTurn2.meanData-leftForceRTurn2.sdData, leftForceRTurn2.meanData+leftForceRTurn2.sdData,
    alpha=0.5, edgecolor='#00966C', facecolor='#00966C', label = turnsRound2.config)
ax2.set_title('Turns to the Right')
ax2.legend()
ax2.set_ylabel('Left Force (N)')

