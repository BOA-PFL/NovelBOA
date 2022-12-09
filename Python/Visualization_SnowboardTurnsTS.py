# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 14:03:26 2022

Visualization tool to compare the force profiles from turns for two models

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

#fPath = 'C:\\Users\\daniel.feeney\\Boa Technology Inc\\PFL Team - General\\Testing Segments\\Snow Performance\\SB_2DialTakeDown_Mar2022\\Forces\\'
fPath = 'C:\\Users\\daniel.feeney\\Boa Technology Inc\\PFL Team - General\\Testing Segments\\Snowboard Pilot2022\\BurtonTestDec22\\'
fileExt = r".txt"
entries = [fName for fName in os.listdir(fPath) if fName.endswith(fileExt)]
entries = os.listdir(fPath)

# Choose two files to compare
fName1 = entries[0]
fName2 = entries[2]

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
    toestarts: list
    heelstarts: list
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

def findToeTurns(toeForce, heelForce):
    """
    This function finds the frames where toe turns begin in snowboarding, 
    which is assumed to be at the
    point where force under the toes exceeds force under the heel.The mean of 
    the following 200 indices for toes must also be greater than heel
    
    Parameters
    ----------
    toeForce : Pandas column
        Time series of force data under the toes measured from insoles
    heelForce : Pandas column
        Time sereis of force data under the heel measured from insoles

    Returns
    -------
    toeTurns : list
        Index of frames in where toe turns started 

    """
    toeTurns = []
    for step in range(len(toeForce)-1):
        if toeForce[step] <= heelForce[step] and toeForce[step + 1] > heelForce[step + 1] and np.mean(toeForce[step:step+200]) > np.mean(heelForce[step:step+200]):
            toeTurns.append(step)
    return toeTurns

#Find takeoff from FP when force goes from above thresh to 0
def findHeelTurns(toeForce, heelForce):
    """
    This function finds the frames where heel turns begin in snowboarding, which is assumed to be at the
    point where force under the heels exceeds force under the toes.

    Parameters
    ----------
    toeForce : numpy array
        Time series of force data under the toes measured from insoles
    heelForce : numpy array
        Time sereis of force data under the heel measured from insoles

    Returns
    -------
    heelTurns : list
        Index of frames in where heel turns started

    """
    heelTurns = []
    for step in range(len(toeForce)-1):
        if heelForce[step] <= toeForce[step] and heelForce[step + 1] > toeForce [step + 1]and np.mean(heelForce[step:step+200]) > np.mean(toeForce[step:step+200]):
            heelTurns.append(step)
    return heelTurns


def cleanTurns(turnIndices):
    
    """
    Finds if there are two turn initiations from the same foot close together
    and removes the first one. Assumes 0.5* median difference in turn start is
    enough to differentiate
    """
    lowThreshold = np.median(np.diff(turnIndices)) * 0.5
    trueIndices = [idx for idx, x in enumerate(list(np.diff(turnIndices))) if x < lowThreshold]
    
    cleanList = [x for idx, x in enumerate(turnIndices) if idx not in trueIndices]
    
    return(cleanList)

def makeVizPlot(inputDF, inputToeTurns, inputHeelTurns):
    """
    makes plot to check data
    """
    fig, (ax, ax1) = plt.subplots(1,2)
    ax.plot(inputDF.bothToes_Filt, label = 'Total Toe Force')
    ax.plot(inputDF.bothHeels_Filt, label = 'Total Heel Force')
    ax.vlines(x = inputToeTurns, ymin = 0, ymax = 1000,
      color = 'k', label = 'Toe Turn Start', linewidth=3.0, ls='--')
    ax.legend() 
    ax.set_title('Toe Turns')
    ax1.plot(inputDF.bothToes_Filt, label = 'Front Toe Force')
    ax1.plot(inputDF.bothHeels_Filt, label = 'Front Heel Force')
    ax1.vlines(x = inputHeelTurns, ymin = 0, ymax = 1000,
      color = 'k', label = 'Heel Turn Start', linewidth=3.0, ls='--')
    ax1.legend()       
    ax1.set_title('Heel Turns')
    


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
    dat = pd.read_csv(fPath + fName,sep='\t', skiprows = 4, header = None, index_col = False, usecols = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
    dat.columns = ['Time', 'FrontHeel', 'FrontMedial','FrontLateral','FrontTotal', 'Time2', 'RearLateral','RearMedial','RearHeel','RearTotal']
    info = fName.split(sep = "/")[-1]
    subName = info.split(sep = "_")[0]
    stance = info.split(sep = "_")[1]
    configName = info.split(sep = "_")[2]
    #trialNo = info.split(sep = "_")[3]
    
    if stance == 'Goofy':
        dat.rename(columns = {'FrontHeel':'RearHeel', 'FrontMedial':'RearMedial', 'FrontLateral':'RearLateral', 'FrontTotal':'RearTotal', 'RearLateral':'FrontLateral', 'RearMedial':'FrontMedial', 'RearHeel':'FrontHeel','RearTotal':'FrontTotal'}, inplace = True)
    
    dat = dat.apply(pd.to_numeric)
    dat['FrontToes'] = dat.FrontMedial + dat.FrontLateral
    dat['RearToes'] = dat.RearMedial + dat.RearLateral
    
    dat['bothToes'] = dat.FrontToes + dat.RearToes
    dat['bothHeels'] = dat.FrontHeel + dat.RearHeel
    fig, ax = plt.subplots()
    
    # Select the period you want to pull data from (i.e. when they were riding)
    ax.plot(dat.FrontTotal, label = 'Front Total Force')
    ax.plot(dat.RearTotal, label = 'Rear Total Force')
    fig.legend()
    print('Select start and end of analysis trial')
    pts = np.asarray(plt.ginput(2, timeout=-1))
    plt.close()
    dat = dat.iloc[int(np.floor(pts[0,0])) : int(np.floor(pts[1,0])),:]
    dat = dat.reset_index()
    
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
    
    realToeStart = findToeTurns(dat.bothToes_Filt, dat.bothHeels_Filt)
    realHeelStart = findHeelTurns(dat.bothToes_Filt, dat.bothHeels_Filt)
    # clean below by not taking false start turns #        
    realHeelStart = cleanTurns(realHeelStart)
    realToeStart = cleanTurns(realToeStart)
    
    makeVizPlot(dat, realToeStart, realHeelStart)
    answer = messagebox.askyesno("Question","Is data clean?")
    
    if answer == False:
        plt.close('all')
        print('Adding file to bad file list')
        #badFileList.append(fName)
    
    if answer == True:
        plt.close('all')
        print('Proceeding to plotting')
    
    #test = TurnStarts('sub1','DD',[1,2],['a','b'])
    result = TurnStarts(subName, configName, realToeStart, realHeelStart, df = dat)
               
    return(result)

turnsRound1 = calcTurnStarts(fName1)
turnsRound2 = calcTurnStarts(fName2)
x = np.linspace(0,stepLen,stepLen)

## Calculating averaged data ## 

toeForceAvg = calcEnsembleData(turnsRound1.df.bothToes_Filt,  turnsRound1.toestarts, stepLen)
heelForceAvg = calcEnsembleData(turnsRound1.df.bothHeels_Filt, turnsRound1.heelstarts, stepLen)
toeForceAvg2 = calcEnsembleData(turnsRound2.df.bothToes_Filt,  turnsRound2.toestarts, stepLen)
heelForceAvg2 = calcEnsembleData(turnsRound2.df.bothHeels_Filt, turnsRound2.heelstarts, stepLen)

heelForceTT = calcEnsembleData(turnsRound1.df.bothHeels_Filt, turnsRound1.toestarts, stepLen)
heelForceTT2 = calcEnsembleData(turnsRound2.df.bothHeels_Filt, turnsRound2.toestarts, stepLen)


### plotting ### 
fig, (ax1, ax2) = plt.subplots(2)
ax1.plot(x, toeForceAvg.meanData, 'k', color='#DC582A')
ax1.fill_between(x,toeForceAvg.meanData-0.5*toeForceAvg.sdData, toeForceAvg.meanData+0.5*toeForceAvg.sdData,
    alpha=0.5, edgecolor='#DC582A', facecolor='#DC582A', label = turnsRound1.config)
ax1.plot(x, toeForceAvg2.meanData, 'k', color='#00966C')
ax1.fill_between(x,toeForceAvg2.meanData-0.5*toeForceAvg2.sdData, toeForceAvg2.meanData+0.5*toeForceAvg2.sdData,
    alpha=0.5, edgecolor='#00966C', facecolor='#00966C', label = turnsRound2.config)
ax1.legend()
ax1.set_title('Toe Turns')
ax1.set_ylabel('Toe Force (N)')
ax1.set_xlabel('ms')
ax2.plot(x, heelForceAvg.meanData, 'k', color='#DC582A')
ax2.fill_between(x,heelForceAvg.meanData-0.5*heelForceAvg.sdData, heelForceAvg.meanData+0.5*heelForceAvg.sdData,
    alpha=0.5, edgecolor='#DC582A', facecolor='#DC582A', label = turnsRound1.config)
ax2.plot(x, heelForceAvg2.meanData, 'k', color='#00966C')
ax2.fill_between(x,heelForceAvg2.meanData-0.5*heelForceAvg2.sdData, heelForceAvg2.meanData+0.5*heelForceAvg2.sdData,
    alpha=0.5, edgecolor='#00966C', facecolor='#00966C', label = turnsRound2.config)
ax2.set_title('Heel Turns')
ax2.legend()
ax2.set_ylabel('Heel Force (N)')
ax2.set_xlabel('ms')
plt.tight_layout()


fig2, ax3 = plt.subplots(1)
ax3.plot(x, heelForceTT.meanData, 'k', color='#DC582A')
ax3.fill_between(x,heelForceTT.meanData-0.5*heelForceTT.sdData, heelForceTT.meanData+0.5*heelForceTT.sdData,
    alpha=0.5, edgecolor='#DC582A', facecolor='#DC582A', label = turnsRound1.config)
ax3.plot(x, heelForceTT2.meanData, 'k', color='#00966C')
ax3.fill_between(x,heelForceTT2.meanData-0.5*heelForceTT2.sdData, heelForceTT2.meanData+0.5*heelForceTT2.sdData,
    alpha=0.5, edgecolor='#00966C', facecolor='#00966C', label = turnsRound2.config)
ax3.legend()
ax3.set_title('Heel Force During Toe Turns')
ax3.set_ylabel('Heel Force (N)')
ax3.set_xlabel('ms')
plt.tight_layout()
