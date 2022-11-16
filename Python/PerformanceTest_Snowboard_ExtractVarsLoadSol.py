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
from tkinter import messagebox
import scipy

pd.options.mode.chained_assignment = None  # default='warn' set to warn for a lot of warnings
save_on = 0

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
    
def intp_strides(var,landings,GS):
    """
    Function to interpolate the variable of interest across a stride
    (from foot contact to subsiquent foot contact) in order to plot the 
    variable of interest over top each other

    Parameters
    ----------
    var : list or numpy array
        Variable of interest. Can be taken from a dataframe or from a numpy array
    landings : list
        Foot contact indicies
    GS: list 
        Good strides that have been passed through the post-hoc filter process

    Returns
    -------
    intp_var : numpy array
        Interpolated variable to 101 points with the number of columns dictated
        by the number of strides.

    """
    # Preallocate
    intp_var = np.zeros((101,len(GS)-1))
    # Index through the strides
    for ii in range(len(GS)-1):
        dum = var[landings[GS[ii]]:landings[GS[ii]+1]]
        f = scipy.interpolate.interp1d(np.arange(0,len(dum)),dum)
        intp_var[:,ii] = f(np.linspace(0,len(dum)-1,101))
        
    return intp_var

# Read in files

fPath = 'C:\\Users\\daniel.feeney\\Boa Technology Inc\\PFL Team - General\\Testing Segments\\Snow Performance\\SB_2DialTakeDown_Mar2022\\Forces\\'
fileExt = r".txt"
entries = [fName for fName in os.listdir(fPath) if fName.endswith(fileExt)]
entries = os.listdir(fPath)
# Select files for a single subject
#entries = askopenfilenames(initialdir = fPath)

### Initiate Time Series

toeTurns_FrontHeel = []
toeTurns_FrontToes = []
toeTurns_RearHeel = []
toeTurns_RearToes = []
heelTurns_FrontHeel = []
heelTurns_FrontToes = []
heelTurns_RearHeel = []
heelTurns_RearToes = []

### Initiate discrete metrics
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

badFileList = []

# naming #
Subject = []
Config = []
Trial = []

for fName in entries:
    
    ### Loop through each file, use time series of force data to identify turns.
    try:
        #fName = entries
        dat = pd.read_csv(fPath + fName,sep='\t', skiprows = 4, header = None, index_col = False, usecols = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
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
            badFileList.append(fName)
        
        if answer == True:
            plt.close('all')
            print('Estimating point estimates')

        
            for i in range(len(realToeStart)-1):
                ### Loop through toe turns to extract variables. 
                
                try:
                    # #i = 1
                    # ### Time Series
            
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
     

# Create dataframe with all outcome variables and write to csv. Will create a new csv if one doesn't already 
# exist for this data set. Otherwise will append.
outcomes = pd.DataFrame({'Subject':list(Subject), 'Config': list(Config), 'Trial':list(Trial),
                         
                         'MaxToeF_BothToes':list(maxToeF_BothToes), 'MaxTotalF_BothToes':list(maxTotalF_BothToes),
                         'maxToeRFDup_BothToes':list(maxToeRFDup_BothToes), 'maxToeRFDdn_BothToes':list(maxToeRFDdn_BothToes), 'maxTotalRFDup_BothToes':list(maxTotalRFDup_BothToes), 'maxTotalRFDdn_BothToes':list(maxTotalRFDdn_BothToes), 
                         'timeToToePeak_BothtToes':list(timeToToePeak_BothToes), 'timeToTotalPeak_BothToes':list(timeToTotalPeak_BothToes),
                         'avgToeRFD_BothToes':list(avgToeRFD_BothToes), 'avgTotalRFD_BothToes':list(avgTotalRFD_BothToes),
                         'meanHeelContact_BothToes':list(meanHeelContact_BothToes)
                         
                         })

if save_on == 1:
    outfileName = fPath + 'CompiledResults6.csv'
    
    if os.path.exists(outfileName) == False:
        
        outcomes.to_csv(outfileName, mode='a', header=True, index = False)
    
    else:
        outcomes.to_csv(outfileName, mode='a', header=False, index = False)

