# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 13:53:12 2023

Creating video of force traces during skiing

@author: Dan.Feeney
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import signal
from matplotlib.animation import FuncAnimation
from itertools import count
import time

fPath = 'C:\\Users\\daniel.feeney\\Boa Technology Inc\\PFL Team - General\\Testing Segments\\Snow Performance\\\SkiValidation_Dec2022\Loadsol\\'
fileExt = r".txt"
entries = [fName for fName in os.listdir(fPath) if fName.endswith(fileExt)]

fName = entries[105]

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


dat = pd.read_csv(fPath+fName, sep = '	',skiprows = 3, header = 0, index_col = False)
dat.columns = ['Time', 'LHeel', 'LMedial','LLateral','LTotal', 'Time2', 
               'RLateral','RMedial','RHeel','RTotal', 'time2','accX','axxY','accZ','pass']

       
subName = fName.split(sep = "_")[0]
configName = fName.split(sep = "_")[1]
trialNoTmp = fName.split(sep = "_")[2].split(sep=".")[0]

dat['LToes'] = dat.LMedial + dat.LLateral
dat['RToes'] = dat.RMedial + dat.RLateral


shortFName = fName.split('.')[0]
# Load in the trial segmentation variable if it is in the directory
if os.path.exists(fPath+shortFName+'TrialSeg.npy') == True:
    trial_segment_old = np.load(fPath+shortFName+'TrialSeg.npy',allow_pickle=True)
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
    # Save the trial segmentation
    trial_segment = np.array([shortFName,pts])
    np.save(fPath+shortFName+'TrialSeg.npy',trial_segment)

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

dat['RToe_Filt'] = dat.RMedial_Filt + dat.RLateral_Filt
dat['LToe_Filt'] = dat.LMedial_Filt + dat.LLateral_Filt

# plt.figure(1)
# plt.subplot(1,2,1)
# plt.plot(dat.LMedial_Filt)
# plt.plot(dat.LLateral_Filt)
# plt.plot(dat.LHeel_Filt)
# plt.legend(['Medial','Lateral','Heel'])

# plt.subplot(1,2,2)
# plt.plot(dat.RMedial_Filt)
# plt.plot(dat.RLateral_Filt)
# plt.plot(dat.RHeel_Filt)
# plt.legend(['Medial','Lateral','Heel'])

RForce = []
for ind, val in enumerate(dat.RTotal_Filt):
    #print(val)
    if not ind % 10:
        RForce.append(val)
LForce = []
for ind, val in enumerate(dat.LTotal_Filt):
    #print(val)
    if not ind % 10:
        LForce.append(val)
    
# plt.plot(RForce)
# plt.plot(LForce)


fig, ax = plt.subplots(dpi=120)
plt.style.use('fivethirtyeight')
x_values = []
y_values = []
z_values = []
counter = 0
index = count()

def animate(i):
    
    #print(counter)
    
    x = next(index) # counter or x variable -> index
    #counter = next(index)
    x_values.append(x)
    y = RForce[x]
    z = LForce[x]

    # append values to keep graph dynamic
    # this can be replaced with reading values from a csv files also
    # or reading values from a pandas dataframe
    y_values.append(y)
    z_values.append(z)
    
    ax.clear()
    ax.plot(x_values, y_values,linestyle='--', color = 'red')
    ax.plot(x_values, z_values,linestyle='--', color = 'blue')
    
    ax.set_xlabel("Index")
    ax.set_ylabel("Force (N)")
    ax.legend(['Right Force', 'Left Force'])
    plt.title('Ski Carving Turns')
    plt.tight_layout()
    
    #time.sleep(.05) # keep refresh rate of 0.25 seconds
    
anim = FuncAnimation(plt.gcf(), animate, 1000)
plt.tight_layout()
plt.show()
    
# # saving to m4 using ffmpeg writer
# writervideo = FuncAnimation.to_html5_video()
# anim.save('increasingStraightLine.mp4', writer=writervideo)
# plt.close()
    
    
    