# -*- coding: utf-8 -*-3
"""
Created on Wed Dec 16 16:19:30 2020
When you have both sides of data from loadsol

Creating code to generate peak forces and force asymmetry for athlete reports
in ski data. Note: the downhill ski should have more force than the uphill ski
so when the peak force on the left foot is higher, this is a turn with the 
body facing to the right. This is the assumption in this code & holds for good
skiers. 

"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import scipy.signal as sig
from scipy import signal
from tkinter.filedialog import askopenfilename
from tkinter import messagebox


fPath = 'C:\\Users\\eric.honert\\Boa Technology Inc\\PFL Team - General\\Testing Segments\\Snow Performance\\SkiValidation_Dec2022\Loadsol\\'
fileExt = r".txt"
fName = askopenfilename(initialdir = fPath)


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

#makeTurnPlot(dat, RTurns, 'Right')
## Example data lives in there and work well ## 

# Initiate discrete outcome variables
peakRForce = []
peakLForce = []


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

trimName = fName.split('/')[-1].split('.')[0]

# modify the default parameters of np.load

# Load in the trial segmentation variable if it is in the directory
if os.path.exists(fPath+trimName+'TrialSeg.npy') == True:
    trial_segment_old = np.load(fPath+trimName+'TrialSeg.npy',allow_pickle=True)
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
b, a = signal.butter(4, w, 'low')
dat['LTotal_Filt'] = signal.filtfilt(b, a, dat.LTotal)
dat['RTotal_Filt'] = signal.filtfilt(b, a, dat.RTotal)
dat['RMedial_Filt'] = signal.filtfilt(b, a, dat.RMedial)
dat['RLateral_Filt'] = signal.filtfilt(b, a, dat.RLateral)
dat['LMedial_Filt'] = signal.filtfilt(b, a, dat.LMedial)
dat['LLateral_Filt'] = signal.filtfilt(b, a, dat.LLateral)
dat['LHeel_Filt'] = signal.filtfilt(b, a, dat.LHeel)
dat['RHeel_Filt'] = signal.filtfilt(b, a, dat.RHeel)

# Turn detection
fs = 100 
fc = 0.5
w = fc / (fs / 2)
b, a = signal.butter(2, w, 'low')

Lturn_detect = signal.filtfilt(b, a, dat.LTotal)
Rturn_detect = signal.filtfilt(b, a, dat.RTotal)

Lpeaks,_ = sig.find_peaks(Lturn_detect, prominence=25)
Rpeaks,_ = sig.find_peaks(Rturn_detect, prominence=25)

        
# Clean up the turn detection: ensure they oscillate
if Lpeaks[0] < Rpeaks[0]:
    Lpeaks, Rpeaks = EnsureTurnsAlternate(Lturn_detect,Rturn_detect,Lpeaks,Rpeaks)

elif Lpeaks[0] > Rpeaks[0]:
    Rpeaks, Lpeaks = EnsureTurnsAlternate(Rturn_detect,Lturn_detect,Rpeaks,Lpeaks)
   



for i, value in enumerate(Rpeaks):
# Loop through all cleaned turns to calculate discrete outcome measures.
# using right side only (left turns) first
    ## Extract relevent parameters from a turn here ##
    peakRForce.append(np.max(dat.RTotal_Filt[value-100:value+100]))
    
for i, value in enumerate(Lpeaks):
# Loop through all cleaned turns to calculate discrete outcome measures.
# using right side only (left turns) first
    peakLForce.append(np.max(dat.LTotal_Filt[value-100:value+100]))

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
bar_colors = ['tab:blue','tab:red']


rightAvgForce = np.mean(peakRForce)
leftAvgForce = np.mean(peakLForce)

pctDiff = abs(100 * ((rightAvgForce - leftAvgForce) / ( (rightAvgForce+leftAvgForce) / 2)))
if rightAvgForce > leftAvgForce:
    print("Skier has ", round(pctDiff), " % greater force on right foot (left turns with right foot downhill are stronger")

if leftAvgForce > rightAvgForce:
    print("Skier has ", round(pctDiff), "% greater force on left foot (right turns with left foot downhill are stronger")


fig, (ax) = plt.subplots(1,1)
variab = ('Left', 'Right')
performance = [leftAvgForce, rightAvgForce]
error = [np.std(peakLForce), np.std(peakRForce)]
ax.bar(variab, performance,yerr=error, color = bar_colors)
ax.set_xlabel('Side')
ax.set_title('Peak Force During Turn')
ax.set_ylabel('Force (N)')
variab = ('Left', 'Right')
plt.tight_layout()
plt.show()

if os.path.exists(fPath+trimName+'TrialSeg.npy') == False:
    # Save the trial segmentation
    trial_segment = np.array([trimName,pts])
    np.save(fPath+trimName+'TrialSeg.npy',trial_segment)
