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
from scipy import signal
from tkinter.filedialog import askopenfilename
from tkinter import messagebox


fPath = 'C:\\Users\\daniel.feeney\\Boa Technology Inc\\PFL Team - General\\Testing Segments\\Snow Performance\\SkiValidation_Dec2022\Loadsol\\'
fileExt = r".txt"
fName = askopenfilename(initialdir = fPath)


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

   
RTurns = findRightTurns(dat.RTotal_Filt, dat.LTotal_Filt)
LTurns = findLeftTurns(dat.RTotal_Filt, dat.LTotal_Filt)

RTurns = cleanTurns(RTurns)
LTurns = cleanTurns(LTurns)

makeTurnPlot(dat, LTurns, 'Left Turns')
makeTurnPlot(dat, RTurns, 'Reft Turns')
answer = messagebox.askyesno("Question","Is data clean?")

if answer == False:
    plt.close('all')
    print('Adding file to bad file list')

if answer == True:
    plt.close('all')
    print('Estimating point estimates')


for i, value in enumerate(LTurns):
# Loop through all cleaned turns to calculate discrete outcome measures.
# using right side only (left turns) first

    if i < len(LTurns) - 2:
        # index of peak turn. saved for later in case we want to use it
        pkIdx = np.argmax(dat.RTotal_Filt[LTurns[i]:LTurns[i+1]])
        pkIdx = value + pkIdx
        ## Extract relevent parameters from a turn here ##
        peakRForce.append(np.max(dat.RTotal_Filt[LTurns[i]:LTurns[i+1]]))
    

    
for i, value in enumerate(RTurns):
# Loop through all cleaned turns to calculate discrete outcome measures.
# using right side only (left turns) first

    if i < len(RTurns) - 2:
    # index of peak turn. saved for later in case we want to use it
        pkIdx = np.argmax(dat.LTotal_Filt[RTurns[i]:RTurns[i+1]])
        pkIdx = value + pkIdx
        ## Extract relevent parameters from a turn here ##
        peakLForce.append(np.max(dat.LTotal_Filt[RTurns[i]:RTurns[i+1]]))



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

# Save the trial segmentation
trial_segment = np.array([trimName,pts])
np.save(fPath+trimName+'TrialSeg.npy',trial_segment)
