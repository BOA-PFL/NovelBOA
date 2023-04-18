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
import scipy
import scipy.signal as sig
from scipy import signal
import addcopyfighandler
from tkinter.filedialog import askopenfilenames
from tkinter import messagebox

fPath = 'C:\\Users\\eric.honert\\Boa Technology Inc\\PFL Team - General\\Testing Segments\\Snow Performance\\\SkiValidation_Dec2022\Loadsol\\'
fileExt = r".txt"
entries = [fName for fName in os.listdir(fPath) if fName.endswith(fileExt)]
check_data = 0
freq = 100
save_on = 0
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

def intp_turns(var,turnidx):
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

    Returns
    -------
    intp_var : numpy array
        Interpolated variable to 101 points with the number of columns dictated
        by the number of strides.

    """
    # Preallocate
    intp_var = np.zeros((101,len(turnidx)-1))
    # Index through the strides
    for ii in range(len(turnidx)-1):
        dum = var[turnidx[ii]:turnidx[ii+1]]
        f = scipy.interpolate.interp1d(np.arange(0,len(dum)),dum)
        intp_var[:,ii] = f(np.linspace(0,len(dum)-1,101))
        
    return intp_var


# Initiate discrete outcome variables
OutTotMaxForce = []
OutMedMaxForce = []
OutLatMaxForce = []
OutHeelMaxForce = []
OutToeMaxForce = []
OutTotAvgForce = []

OutMedFracImpulse = []
OutLatFracImpulse = []
OutToeFracImpulse = []
OutHeelFracImpulse = []
OutToeFracImpulseEarly = []
OutHeelFracImpulseLate = []

RFD = []
RFDtime = []

InsTotMinForce = []



OutsideFootProp = []
OutsideFootForce = []
OutsideFootMedialProp = []
OutsideFootMedialForce = []
avgOutsideHeelStartForce = []
avgOutsideHeelStartProp = []
propHeelLate = []
absPropHeelLate = []

cvForce = []

sName = []
cName = []
TurnDir = []
Side = []
timeToPeak = []
badFileList = []
trialNo = []

for ii, entry in enumerate(entries[17:-1]):
    try:
        fName = entry
        print(fName)
        # Loop through files and use time series force data to identify turns
        #fName = entries[2]
        dat = pd.read_csv(fPath+fName, sep = '	',skiprows = 3, header = 0, index_col = False)
        dat.columns = ['Time', 'LHeel', 'LMedial','LLateral','LTotal', 'Time2', 
                       'RLateral','RMedial','RHeel','RTotal', 'time2','accX','axxY','accZ','pass']
        
        #dat.columns = ['Time','RHeel','RLateral','RMedial','RTotal','Time2','se','ei','ni','te']
        #use above if one side only
                
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
        
        # plt.figure(ii)
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
        
        # plt.close()
        
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
         
        
        if check_data == 1:
            plt.figure
            plt.subplot(2,1,1)
            plt.plot(Lturn_detect)
            plt.plot(Lpeaks,Lturn_detect[Lpeaks],'ks')
            plt.legend(['Left Turn Sig','Det Turn'])
            plt.subplot(2,1,2)
            plt.plot(Rturn_detect,'r')
            plt.plot(Rpeaks,Rturn_detect[Rpeaks],'ks')
            plt.legend(['Right Turn Sig','Det Turn'])
            
            answer = messagebox.askyesno("Question","Is data clean?")
            plt.close('all')
            
            if answer == False:
                plt.close('all')
                print('Adding file to bad file list')
                badFileList.append(fName)
        else:
            answer = True
            
        if answer == True:
            print('Estimating point estimates')
        
            for jj, value in enumerate(Lpeaks):
                # Loop through all cleaned turns to calculate discrete outcome measures.
                # using right side only (left turns). No longer assuming left to right
                # transitions are constant and only using data from one side. 
                # Extend this to Left Turns as well
                # variables of interest: outside force (downhill foot peak force) higher is better
                # Outside medial peak force (higher is better)
                # Outside heel force lower is better with an average higher force on forefoot Tricky**
    
                try:
                    # Look at a 200 frame window for the true max force
                    # Look at the Outside (or downhill) ski
                    if value < 100:
                        window = value
                    else:
                        window = 100
                    
                    OutTotMaxForce.append(np.nanmax(dat.LTotal_Filt[value-window:value+window]))
                    OutMedMaxForce.append(np.nanmax(dat.LMedial_Filt[value-window:value+window]))
                    OutLatMaxForce.append(np.nanmax(dat.LLateral_Filt[value-window:value+window]))
                    OutToeMaxForce.append(np.nanmax(dat.LToe_Filt[value-window:value+window]))
                    OutHeelMaxForce.append(np.nanmax(dat.LHeel_Filt[value-window:value+window]))
                    
                    # Look at the Inside (or uphill) ski
                    InsTotMinForce.append(np.nanmax(dat.RTotal_Filt[value-window:value+window]))
                    
                    # Examine a 0.5 second window around the peak of the turn
                    # for fraction of force  exerted on the foot of the outside (downhill) ski
                    TotImpulse = sum(dat.LTotal_Filt[value-25:value+25])
                    OutMedFracImpulse.append(sum(dat.LMedial_Filt[value-25:value+25])/TotImpulse)
                    OutLatFracImpulse.append(sum(dat.LLateral_Filt[value-25:value+25])/TotImpulse)
                    OutToeFracImpulse.append(sum(dat.LToe_Filt[value-25:value+25])/TotImpulse)
                    OutHeelFracImpulse.append(sum(dat.LHeel_Filt[value-25:value+25])/TotImpulse)
                    
                    if jj < len(Lpeaks)-1:
                        # Proportion of the force on the Posterior Portion during
                        # the end of the turn (last 25%)
                        idx_max = np.argmax(dat.LTotal_Filt[value-window:value+window])+value-window
                        idx_nextmin = np.argmin(dat.LTotal_Filt[idx_max:Lpeaks[jj+1]])+idx_max
                        
                        idx75 = round(0.5*(idx_nextmin-idx_max))+idx_max
                        OutHeelFracImpulseLate.append(sum(dat.LHeel_Filt[idx75:idx_nextmin])/sum(dat.LTotal_Filt[idx75:idx_nextmin]))
                    else:
                        OutHeelFracImpulseLate.append(np.nan)
                        
                    if jj>0:
                        # Rate of force development
                        # Find the index of the true local maxima
                        idx_max = np.argmax(dat.LTotal_Filt[value-window:value+window])+value-window
                        idx_min = np.argmin(dat.LTotal_Filt[Lpeaks[jj-1]:idx_max])+Lpeaks[jj-1]
                        RFD.append(np.mean(np.gradient(dat.LTotal_Filt[idx_min:idx_max],1/freq)))
                        RFDtime.append((idx_max-idx_min)/freq)
                        
                        # Proportion of the force on the Anterior Portion during
                        # The start of the turn (first 25%)
                        idx25 = round(0.5*(idx_max-idx_min))+idx_min
                        OutToeFracImpulseEarly.append(sum(dat.LToe_Filt[idx_min:idx25])/sum(dat.LTotal_Filt[idx_min:idx25]))
                        
                    else:
                        RFD.append(np.nan)
                        RFDtime.append(np.nan)
                        OutToeFracImpulseEarly.append(np.nan)
                        
                    
                    if jj>0 and jj < len(Lpeaks)-1:
                        OutTotAvgForce.append(np.mean(dat.LTotal_Filt[idx_min:idx_nextmin]))
                    else:
                        OutTotAvgForce.append(np.nan)
                        
                    TurnDir.append('R')
                    sName.append(subName)
                    cName.append(configName)
                    trialNo.append(trialNoTmp)
                    
                    
                except Exception as e: print(e)
                    
            for jj, value in enumerate(Rpeaks):
                # Loop through all cleaned Right Turns to calculate discrete outcome measures.
                # using right side only (left turns). No longer assuming left to right
                # transitions are constant and only using data from one side. 
                # Extend this to Left Turns as well
                # variables of interest: outside force (downhill foot peak force) higher is better
                # Outside medial peak force (higher is better)
                # Outside heel force lower is better with an average higher force on forefoot Tricky**
    
                try:
                    if value < 100:
                        window = value
                    else:
                        window = 100
                    # Look at a 200 frame window for the true max force
                    OutTotMaxForce.append(np.nanmax(dat.RTotal_Filt[value-window:value+window]))
                    OutMedMaxForce.append(np.nanmax(dat.RMedial_Filt[value-window:value+window]))
                    OutLatMaxForce.append(np.nanmax(dat.RLateral_Filt[value-window:value+window]))
                    OutToeMaxForce.append(np.nanmax(dat.RToe_Filt[value-window:value+window]))
                    OutHeelMaxForce.append(np.nanmax(dat.RHeel_Filt[value-window:value+window]))
                    
                    # Look at the Inside (or uphill) ski
                    InsTotMinForce.append(np.nanmax(dat.LTotal_Filt[value-window:value+window]))
                    
                    # Examine a 0.5 second window around the peak of the turn
                    # for fraction of force  exerted on the foot of the outside (downhill) ski
                    TotImpulse = sum(dat.RTotal_Filt[value-25:value+25])
                    OutMedFracImpulse.append(sum(dat.RMedial_Filt[value-25:value+25])/TotImpulse)
                    OutLatFracImpulse.append(sum(dat.RLateral_Filt[value-25:value+25])/TotImpulse)
                    OutToeFracImpulse.append(sum(dat.RToe_Filt[value-25:value+25])/TotImpulse)
                    OutHeelFracImpulse.append(sum(dat.RHeel_Filt[value-25:value+25])/TotImpulse)
                    
                    if jj < len(Rpeaks)-1:
                        # Proportion of the force on the Posterior Portion during
                        # the end of the turn (last 25%)
                        idx_max = np.argmax(dat.RTotal_Filt[value-window:value+window])+value-window
                        idx_nextmin = np.argmin(dat.RTotal_Filt[idx_max:Rpeaks[jj+1]])+idx_max
                        
                        idx75 = round(0.5*(idx_nextmin-idx_max))+idx_max
                        OutHeelFracImpulseLate.append(sum(dat.RHeel_Filt[idx75:idx_nextmin])/sum(dat.RTotal_Filt[idx75:idx_nextmin]))
                    else:
                        OutHeelFracImpulseLate.append(np.nan)
                    
                    # Rate of force development
                    if jj>0:
                        # Find the index of the true local maxima
                        idx_max = np.argmax(dat.RTotal_Filt[value-window:value+window])+value-window
                        idx_min = np.argmin(dat.RTotal_Filt[Rpeaks[jj-1]:idx_max])+Rpeaks[jj-1]
                        RFD.append(np.mean(np.gradient(dat.RTotal_Filt[idx_min:idx_max],1/freq)))
                        RFDtime.append((idx_max-idx_min)/freq)
                        
                        # Proportion of the force on the Anterior Portion during
                        # The start of the turn (first 25%)
                        idx25 = round(0.5*(idx_max-idx_min))+idx_min
                        OutToeFracImpulseEarly.append(sum(dat.RToe_Filt[idx_min:idx25])/sum(dat.RTotal_Filt[idx_min:idx25]))
                        
                    else:
                        RFD.append(np.nan)
                        RFDtime.append(np.nan)
                        OutToeFracImpulseEarly.append(np.nan)
                    
                    if jj>0 and jj < len(Rpeaks)-1:
                        OutTotAvgForce.append(np.mean(dat.RTotal_Filt[idx_min:idx_nextmin]))
                    else:
                        OutTotAvgForce.append(np.nan)

                    TurnDir.append('L')
                    sName.append(subName)
                    cName.append(configName)
                    trialNo.append(trialNoTmp)
                    
                except Exception as e: print(e)
                        
    except Exception as e: print(e)


# Create data frame with all outcome measures and export to csv. Will create a new csv if one does not exist for this dataset. 
# Otherwise will append.

outcomes = pd.DataFrame({'Subject':list(sName),'Config':list(cName),'TurnDir': list(TurnDir), 'TrialNo':list(trialNo),
                                 'OutTotMaxForce':list(OutTotMaxForce), 'OutMedMaxForce':list(OutMedMaxForce),
                                 'OutLatMaxForce':list(OutLatMaxForce), 'OutHeelMaxForce':list(OutHeelMaxForce),
                                 'OutToeMaxForce':list(OutToeMaxForce), 'OutMedFracImpulse':list(OutMedFracImpulse),
                                 'OutLatFracImpulse':list(OutLatFracImpulse), 'OutToeFracImpulse':list(OutToeFracImpulse),
                                 'OutHeelFracImpulse':list(OutHeelFracImpulse), 'RFD':list(RFD), 'RFDtime':list(RFDtime),
                                 'OutHeelFracImpulseLate': list(OutHeelFracImpulseLate), 'OutToeFracImpulseEarly': list(OutToeFracImpulseEarly),
                                 'OutTotAvgForce': list(OutTotAvgForce)
                                 })

         
  
outfileName = fPath + 'CompiledResultsTest3.csv'

if save_on == 1:
    if os.path.exists(outfileName) == False:
        
        outcomes.to_csv(outfileName, mode='a', header=True, index = False)
    
    else:
        outcomes.to_csv(outfileName, mode='a', header=False, index = False) 
else:
    print('Complete')


