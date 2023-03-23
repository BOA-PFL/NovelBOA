# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 11:34:24 2023

Visualization of two representative turns from alpine skiing to show difference in BOA vs. Buckles


@author: Dan.Feeney
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import scipy.signal as sig
from scipy import signal

fPath = 'C:\\Users\\daniel.feeney\\Boa Technology Inc\\PFL Team - General\\Testing Segments\\Snow Performance\\\SkiValidation_Dec2022\Loadsol\\'
fNameBuckle = 'C:\\Users\\daniel.feeney\\Boa Technology Inc\\PFL Team - General\\Testing Segments\\Snow Performance\\\SkiValidation_Dec2022\Loadsol\\S04_Buckle_4.txt'
fNameBOA = 'C:\\Users\\daniel.feeney\\Boa Technology Inc\\PFL Team - General\\Testing Segments\\Snow Performance\\\SkiValidation_Dec2022\Loadsol\\S04_BOA_1.txt'


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
    #plt.suptitle(subName)

def EnsureTurnsAlternate(turndet1,turndet2,peaks1,peaks2):
    """
    This function takes 2 signals that alternate (such as left and right) with
    events that have been detected within those signals. It makes sure that the
    detected events (or peaks) oscillates between the two signals. If there are
    multiple events detected between events from the other signal, the event
    with the highest peak will be kept and the others eliminated.

    Parameters
    ----------
    turndet1 : numpy array
        Signal associated with the detected peaks from peaks1
    turndet2 : numpy array
        Signal associated with the detected peaks from peaks2
    peaks1 : numpy array
        Peaks detected from the function "find_peaks" (scipy library) from turndet1
    peaks2 : numpy array
        Peaks detected from the function "find_peaks" (scipy library) from turndet2

    Returns
    -------
    peaks1 : numpy array
        Cleaned peaks1 that oscillates with peaks2
    peaks2 : numpy array
        Cleaned peaks1 that oscillates with peaks1

    """
    # If the are multiple detected peaks after the last turn detection
    if peaks1[-2] > peaks2[-1]:
        idx_multiple = np.where(peaks1>peaks2[-1])[0]
        idx = np.argmax(turndet1[peaks1[idx_multiple]])
        # Remove the other detected peaks
        idx_remove = np.delete(idx_multiple,idx)
        peaks1 = np.delete(peaks1,idx_remove)
    
    if peaks2[-2] > peaks1[-1]:
        idx_multiple = np.where(peaks2>peaks1[-1])[0]
        idx = np.argmax(turndet2[peaks2[idx_multiple]])
        # Remove the other detected peaks
        idx_remove = np.delete(idx_multiple,idx)
        peaks2 = np.delete(peaks2,idx_remove)
    
    # Ensure that the entire signal oscillates between peaks1 and peaks2
    jj = 0
    while jj <= len(peaks2)-2:
        if peaks1[jj] < peaks2[jj] and peaks1[jj+1] > peaks2 [jj] and peaks1[jj+1] < peaks2 [jj+1]:
            # A normal turn
            jj = jj+1
        else:
            if peaks2[jj+1] < peaks1[jj+1]:
                # Multiple detected peaks from the following signal before a peak from the leading signal
                idx_multiple = np.where((peaks2>peaks1[jj])*(peaks2<peaks1[jj+1])==True)[0]
                # Figure out which of multiple peaks is higher
                idx = np.argmax(turndet2[peaks2[idx_multiple]])
                # Remove the other detected peaks
                idx_remove = np.delete(idx_multiple,idx)
                peaks2 = np.delete(peaks2,idx_remove)
            else:
                # Multiple detected peaks from the leading signal before a peak from the following signal
                if jj == 0:
                    idx_multiple = np.where(peaks1<peaks2[jj])[0]  
                else:
                    idx_multiple = np.where((peaks1>peaks2[jj-1])*(peaks1<peaks2[jj])==True)[0]    
                idx = np.argmax(turndet1[peaks1[idx_multiple]])
                # Remove the other detected peaks
                idx_remove = np.delete(idx_multiple,idx)
                peaks1 = np.delete(peaks1,idx_remove)
    return(peaks1,peaks2)


boaDat = pd.read_csv(fNameBOA, sep = '	',skiprows = 3, header = 0, index_col = False)
boaDat.columns = ['Time', 'LHeel', 'LMedial','LLateral','LTotal', 'Time2', 
                       'RLateral','RMedial','RHeel','RTotal', 'time2','accX','axxY','accZ','pass']

buckDat = pd.read_csv(fNameBuckle, sep = '	',skiprows = 3, header = 0, index_col = False)
buckDat.columns = ['Time', 'LHeel', 'LMedial','LLateral','LTotal', 'Time2', 
                       'RLateral','RMedial','RHeel','RTotal', 'time2','accX','axxY','accZ','pass']
            
        
boaDat['LToes'] = boaDat.LMedial + boaDat.LLateral
boaDat['RToes'] = boaDat.RMedial + boaDat.RLateral
buckDat['RToes'] = buckDat.RMedial + buckDat.RLateral
buckDat['LToes'] = buckDat.LMedial + buckDat.LLateral

## load files with segmentation for boa and buckles
shortboa = fNameBOA.split('.')[-2].split('\\')[-1]
trial_segment_old = np.load(fPath+shortboa+'TrialSeg.npy',allow_pickle=True)
trialStart = trial_segment_old[1][0,0]
trialEnd = trial_segment_old[1][1,0]
boaDat = boaDat.iloc[int(np.floor(trialStart)) : int(np.floor(trialEnd)),:]
boaDat = boaDat.reset_index()

shortbuck = fNameBuckle.split('.')[-2].split('\\')[-1]
trial_segment_old2 = np.load(fPath+shortbuck+'TrialSeg.npy',allow_pickle=True)
trialStart2 = trial_segment_old2[1][0,0]
trialEnd2 = trial_segment_old2[1][1,0]
buckDat = buckDat.iloc[int(np.floor(trialStart2)) : int(np.floor(trialEnd2)),:]
buckDat = buckDat.reset_index()

## filter data
fs = 100 
fc = 6
w = fc / (fs / 2)
b, a = signal.butter(2, w, 'low')
boaDat['LTotal_Filt'] = signal.filtfilt(b, a, boaDat.LTotal)
boaDat['RTotal_Filt'] = signal.filtfilt(b, a, boaDat.RTotal)
boaDat['RMedial_Filt'] = signal.filtfilt(b, a, boaDat.RMedial)
boaDat['RLateral_Filt'] = signal.filtfilt(b, a, boaDat.RLateral)
boaDat['LMedial_Filt'] = signal.filtfilt(b, a, boaDat.LMedial)
boaDat['LLateral_Filt'] = signal.filtfilt(b, a, boaDat.LLateral)
boaDat['LHeel_Filt'] = signal.filtfilt(b, a, boaDat.LHeel)
boaDat['RHeel_Filt'] = signal.filtfilt(b, a, boaDat.RHeel)
boaDat['RToe_Filt'] = boaDat.RMedial_Filt + boaDat.RLateral_Filt
boaDat['LToe_Filt'] = boaDat.LMedial_Filt + boaDat.LLateral_Filt

#buckles
buckDat['LTotal_Filt'] = signal.filtfilt(b, a, buckDat.LTotal)
buckDat['LTotal_Filt'] = signal.filtfilt(b, a, buckDat.LTotal)
buckDat['RTotal_Filt'] = signal.filtfilt(b, a, buckDat.RTotal)
buckDat['RMedial_Filt'] = signal.filtfilt(b, a, buckDat.RMedial)
buckDat['RLateral_Filt'] = signal.filtfilt(b, a, buckDat.RLateral)
buckDat['LMedial_Filt'] = signal.filtfilt(b, a, buckDat.LMedial)
buckDat['LLateral_Filt'] = signal.filtfilt(b, a, buckDat.LLateral)
buckDat['LHeel_Filt'] = signal.filtfilt(b, a, buckDat.LHeel)
buckDat['RHeel_Filt'] = signal.filtfilt(b, a, buckDat.RHeel)
buckDat['RToe_Filt'] = buckDat.RMedial_Filt + buckDat.RLateral_Filt
buckDat['LToe_Filt'] = buckDat.LMedial_Filt + buckDat.LLateral_Filt

## turn detection

# Turn detection
fs = 100 
fc = 0.5
w = fc / (fs / 2)
b, a = signal.butter(2, w, 'low')

# BOA turns
Lturn_detect = signal.filtfilt(b, a, boaDat.LTotal)
Rturn_detect = signal.filtfilt(b, a, boaDat.RTotal)
Lpeaksboa,_ = sig.find_peaks(Lturn_detect, prominence=25)
Rpeaksboa,_ = sig.find_peaks(Rturn_detect, prominence=25)    
# Clean up the turn detection: ensure they oscillate
if Lpeaksboa[0] < Rpeaksboa[0]:
    Lpeaksboa, Rpeaksboa = EnsureTurnsAlternate(Lturn_detect,Rturn_detect,Lpeaksboa,Rpeaksboa)

elif Lpeaksboa[0] > Rpeaksboa[0]:
    Rpeaksboa, Lpeaksboa = EnsureTurnsAlternate(Rturn_detect,Lturn_detect,Rpeaksboa,Lpeaksboa)

# buckle turns
Lturn_detect = signal.filtfilt(b, a, buckDat.LTotal)
Rturn_detect = signal.filtfilt(b, a, buckDat.RTotal)
Lpeaksbuck,_ = sig.find_peaks(Lturn_detect, prominence=25)
Rpeaksbuck,_ = sig.find_peaks(Rturn_detect, prominence=25)   
# Clean up the turn detection: ensure they oscillate
if Lpeaksbuck[0] < Rpeaksbuck[0]:
    Lpeaksbuck, Rpeaksbuck = EnsureTurnsAlternate(Lturn_detect,Rturn_detect,Lpeaksbuck,Rpeaksbuck)

elif Lpeaksboa[0] > Rpeaksboa[0]:
    Rpeaksbuck, Lpeaksbuck = EnsureTurnsAlternate(Rturn_detect,Lturn_detect,Rpeaksbuck,Lpeaksbuck)


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


## plotting ##
lboastart = 50
lboastop = 350
lbuckstart = 0
lbuckstop = 300
plt.figure
#plt.subplot(1,2,1)
plt.plot(np.arange(0,(lboastop-lboastart)),boaDat.LToe_Filt[lboastart:lboastop],color = '#00966C',label = 'Left BOA Forefoot')
#plt.plot(np.arange(0,(lboastop-lboastart)),boaDat.LHeel_Filt[lboastart:lboastop],color = '#00966C', linestyle = 'dashed',label = 'Left BOA Heel')
plt.plot(np.arange(0,(lbuckstop-lbuckstart)),buckDat.LToe_Filt[lbuckstart:lbuckstop],color = '#53565A', label = 'Left Buckle Forefoot')
#plt.plot(np.arange(0,(lbuckstop-lbuckstart)),buckDat.LHeel_Filt[lbuckstart:lbuckstop],color = '#53565A', linestyle = 'dashed',label = 'Left Buckle Heel')
plt.ylabel('Force (N)')
plt.xlabel('Time (cs)')
plt.title('Right Turn')
plt.legend()


## option 2
lboastart = 150
lboastop = 550
lbuckstart = 100
lbuckstop = 500
plt.figure
#plt.subplot(1,2,1)
plt.plot(np.arange(0,(lboastop-lboastart)),boaDat.LToe_Filt[lboastart:lboastop],color = '#00966C',label = 'Left BOA Forefoot')
plt.plot(np.arange(0,(lboastop-lboastart)),boaDat.RToe_Filt[lboastart:lboastop],color = '#DC582A', linestyle = 'dashed',label = 'Right BOA Forefoot')
plt.plot(np.arange(0,(lbuckstop-lbuckstart)),buckDat.LToe_Filt[lbuckstart:lbuckstop],color = '#53565A', label = 'Left Buckle Forefoot')
plt.plot(np.arange(0,(lbuckstop-lbuckstart)),buckDat.RToe_Filt[lbuckstart:lbuckstop],color = '#000000', linestyle = 'dashed',label = 'Right Buckle Forefoot')
plt.ylabel('Force (N)')
plt.xlabel('Time (cs)')
plt.title('Two Consecutive Turn')
plt.legend()


boastart = 300
boastop = 600
buckstart = 250
buckstop = 550
plt.subplot(1,2,2)
plt.plot(np.arange(0,(boastop-boastart)),boaDat.RToe_Filt[boastart:boastop],'#DC582A', label = 'Right BOA Forefoot')
plt.plot(np.arange(0,(boastop-boastart)),boaDat.RHeel_Filt[boastart:boastop],'#DC582A', linestyle = 'dashed',label = 'Right BOA Forefoot')
plt.plot(np.arange(0,(buckstop-buckstart)),buckDat.RToe_Filt[buckstart:buckstop],color = '#53565A', label = 'Right Buckle Forefoot')
plt.plot(np.arange(0,(buckstop-buckstart)),buckDat.RHeel_Filt[buckstart:buckstop],color = '#53565A', linestyle = 'dashed', label = 'Right Buckle Forefoot')
plt.ylabel('Force (N)')
plt.xlabel('Time (cs)')
plt.legend()
plt.title('Left Turn')
plt.suptitle('Representative Subject Left and Right Turns')



