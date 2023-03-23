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
import scipy.signal as sig
from scipy import signal
import addcopyfighandler
from tkinter.filedialog import askopenfilenames
from tkinter import messagebox

fPath = 'C:/Users/Kate.Harrison/Boa Technology Inc/PFL Team - General/Testing Segments/Snow Performance/PFLMechanistic_StepOn_March2023/dataPDO/'
fileExt = r".txt"
entries = [fName for fName in os.listdir(fPath) if fName.endswith(fileExt)]
freq = 100

check_data = 1
save_on = 1
#entries = askopenfilenames(initialdir = fPath)



def findToeTurns(ToeForce, HeelForce):
    """
    Find start of toe turns. Defined as when force on Toe  
    exceeds force on heel.

    Parameters
    ----------
    ToeForce : Pandas Column
        Time series of force data under the toes.
    HeelForce : Pandas Column
        Time series of force data under the heel. 

    Returns
    -------
    ToeTurns : List
        Index of frames where right toe turns started. 

    """
    ToeTurns = []
    for step in range(len(ToeForce)-1):
        if ToeForce[step] <= HeelForce[step] and ToeForce[step + 1] > HeelForce[step + 1] and np.mean(ToeForce[step:step+200]) > np.mean(HeelForce[step:step+200]):
            ToeTurns.append(step)
    return ToeTurns

def findHeelTurns(ToeForce, HeelForce):
    """
    Find start of Heel turns. Defined as when force on Heel 
    exceeds force on the Toes.

    Parameters
    ----------
    ToeForce : Pandas Column
        Time series of force data under the toes. 
    HeelForce : Pandas Column
        Time series of force data under the heels

    Returns
    -------
    HeelTurns : list
        Index of frames where heel turns started. 

    """
    HeelTurns = []
    for step in range(len(HeelForce)-1):
        if HeelForce[step] <= ToeForce[step] and HeelForce[step + 1] > ToeForce[step + 1] and np.mean(HeelForce[step:step+200]) > np.mean(ToeForce[step:step+200]) :
            HeelTurns.append(step)
    return HeelTurns



def cleanTurns(turnIndices):
    """
    Finds if there are two turn initiations on the same side (toe or heel) close together
    and removes the first one
    """
    lowThreshold = np.median(np.diff(turnIndices)) * 0.6
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
        Toe or Heel.

    Returns
    -------
    Plot

    """
    plt.figure()
    plt.plot(inputDF.RTotal_Filt, label = "Toe")
    plt.plot(inputDF.LTotal_Filt, label = "Heel")
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


# Initiate discrete outcome variables
ToeTotMaxForce = []
FrontToeMaxForce = []
BackToeMaxForce = []
FrontHeelMaxForce = []
FrontTotalMaxForce = []
BackTotalMaxForce = []
HeelTotMaxForce = []
BackHeelMaxForce = []


ToeTotAvgForce = []
HeelTotAvgForce = []

FrontToeFracImpulse = []
FrontHeelFracImpulse = []
TotalToeFracImpulse = []
TotalHeelFracImpulse = []
BackToeFracImpulse = []
BackHeelFracImpulse = []
FrontFootFracImpulseEarly = []
BackFootFracImpulseLate = []

RFD = []
RFDtime = []

# InsTotMinForce = []

# ToeProp = []
# ToeForce = []
# HeelProp = []
# HeelForce = []
# FrontToeProp = []
# FrontToeForce = []
# FrontHeelProp = []
# FrontHeelForce = []
# BackToeProp = []
# BackToeForce = []
# BackHeelProp = []
# BackHeelForce = []
# avgFrontFootStartForce = []
# avgFrontFootStartProp = []
# propBackFootLate = []
# absPropBackFootLate = []

cvForce = []

sName = []
cName = []
TurnDir = []
Side = []
timeToPeak = []
badFileList = []
trialNo = []



for ii, fName in enumerate(entries):
    try:
        ToeTurn_ToeTS = []
        ToeTurn_HeelTS = []
        HeelTurn_HeelTS = []
        HeelTurn_ToeTS = []
        
        #ii = 2
        #fName = entries[ii]
        print(fName)
        # Loop through files and use time series force data to identify turns
        subName = fName.split(sep = "_")[0]
        foot = fName.split(sep = '_')[1]
        configName = fName.split(sep = "_")[2]
        trialNoTmp = fName.split(sep = "_")[3].split(sep=".")[0]
        
        dat = pd.read_csv(fPath+fName, sep = '	',skiprows = 3, header = 0, index_col = False)
        if dat.shape[1] == 15:
            dat = dat.drop(columns = ['Unnamed: 14'])
        
        if foot == 'Regular':
            dat.columns = ['Time', 'FrontHeel', 'FrontMedial','FrontLateral','FrontTotal', 'Time2', 
                       'RearLateral','RearMedial','RearHeel','RearTotal', 'time2','accX','accY','accZ']
        else:
            dat.columns = ['Time', 'RearHeel', 'RearMedial', 'RearLateral', 'RearTotal', 'Time2',
                           'FrontLateral', 'FrontMedial', 'FrontHeel', 'FrontTotal', 'time2', 'accX', 'accY', 'accZ']   
        
        #dat.columns = ['Time','RHeel','RLateral','RMedial','RTotal','Time2','se','ei','ni','te']
        #use above if one side only
        dat['FrontToes'] = dat.FrontMedial + dat.FrontLateral
        dat['RearToes'] = dat.RearMedial + dat.RearLateral
        dat['TotalToes'] = dat.FrontToes + dat.RearToes
        dat['TotalHeel'] = dat.FrontHeel + dat.RearHeel
        dat['TotalForce'] = dat.RearTotal + dat.FrontTotal
        
        
          
        
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
            ax.plot(dat.TotalToes, label = 'Total Toe Force')
            ax.plot(dat.TotalHeel, label = 'Total Heel Force')
            fig.legend()
            print('Select start and end of analysis trial')
            pts = np.asarray(plt.ginput(2, timeout=-1))
            plt.close()
            # downselect the region of the dataframe you selected from above 
            dat = dat.iloc[int(np.floor(pts[0,0])) : int(np.floor(pts[1,0])),:]
            dat = dat.reset_index()
            # Save the trial segmentation
            trial_segment = np.array([shortFName,pts], dtype = object)
            np.save(fPath+shortFName+'TrialSeg.npy',trial_segment)
        
        fs = 100 
        fc = 6
        w = fc / (fs / 2)
        b, a = signal.butter(2, w, 'low')
        dat['FrontTotal_Filt'] = signal.filtfilt(b, a, dat.FrontTotal)
        dat['RearTotal_Filt'] = signal.filtfilt(b, a, dat.RearTotal)
        dat['FrontToes_Filt'] = signal.filtfilt(b, a, dat.FrontToes)
        dat['FrontHeel_Filt'] = signal.filtfilt(b, a, dat.FrontHeel)
        dat['RearToes_Filt'] = signal.filtfilt(b, a, dat.RearToes)
        dat['RearHeel_Filt'] = signal.filtfilt(b, a, dat.RearHeel)
        dat['TotalToes_Filt'] = signal.filtfilt(b, a, dat.TotalToes)
        dat['TotalHeel_Filt'] = signal.filtfilt(b, a, dat.TotalHeel)
        dat['TotalForce_Filt'] = signal.filtfilt(b, a, dat.TotalForce)
        
       
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
        
        Toeturn_detect = signal.filtfilt(b, a, dat.TotalToes)
        Heelturn_detect = signal.filtfilt(b, a, dat.TotalHeel)
        
        Toepeaks,_ = sig.find_peaks(Toeturn_detect, prominence=100)
        Heelpeaks,_ = sig.find_peaks(Heelturn_detect, prominence=100)
        
        
        
                
        # Clean up the turn detection: ensure they oscillate
        if Toepeaks[0] < Heelpeaks[0]:
            Toepeaks, Heelpeaks = EnsureTurnsAlternate(Toeturn_detect,Heelturn_detect,Toepeaks,Heelpeaks)
        
        elif Toepeaks[0] > Heelpeaks[0]:
            Heelpeaks, Toepeaks = EnsureTurnsAlternate(Heelturn_detect,Toeturn_detect,Heelpeaks,Toepeaks)
         
        
        if check_data == 1:
            plt.figure()
            plt.subplot(2,1,1)
            plt.plot(Toeturn_detect)
            plt.plot(Toepeaks,Toeturn_detect[Toepeaks],'ks')
            plt.legend(['Toe Turn Sig','Det Turn'])
            plt.subplot(2,1,2)
            plt.plot(Heelturn_detect,'r')
            plt.plot(Heelpeaks,Heelturn_detect[Heelpeaks],'ks')
            plt.legend(['Heel Turn Sig','Det Turn'])
            
            answer = messagebox.askyesno("Question","Is data clean?")
            plt.close()
            
            if answer == False:
                plt.close()
                print('Adding file to bad file list')
                badFileList.append(fName)
        else:
            answer = True
            
        if answer == True:
            print('Estimating point estimates')
        
            for jj, value in enumerate(Toepeaks):
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
                    
                    #jj = 8
                    #value = Toepeaks[jj]
                    
                    if value < 100:
                        window = value
                    else:
                        window = 100
                    
                    ToeTotMaxForce.append(np.nanmax(dat.TotalToes_Filt[value-window:value+window]))
                    FrontToeMaxForce.append(np.nanmax(dat.FrontToes_Filt[value-window:value+window]))
                    FrontHeelMaxForce.append(np.nanmax(dat.FrontHeel_Filt[value-window:value+window]))
                    HeelTotMaxForce.append(np.nanmax(dat.TotalHeel_Filt[value-window:value+window]))
                    BackToeMaxForce.append(np.nanmax(dat.RearToes_Filt[value-window:value+window]))
                    BackHeelMaxForce.append(np.nanmax(dat.RearHeel_Filt[value-window:value+window]))
                    FrontTotalMaxForce.append(np.nanmax(dat.FrontTotal_Filt[value-window:value+window]))
                    BackTotalMaxForce.append(np.nanmax(dat.RearTotal_Filt[value-window:value+window]))
                    
                                       
                    # Examine a 0.5 second window around the peak of the turn
                    # for fraction of force  exerted on the foot of the outside (downhill) ski
                    TotImpulse = sum(dat.TotalForce_Filt[value-25:value+25])
                    FrontImpulse = sum(dat.FrontTotal_Filt[value-25:value+25])
                    BackImpulse = sum(dat.RearTotal_Filt[value-25:value+25])
                    FrontToeFracImpulse.append(sum(dat.FrontToes_Filt[value-25:value+25])/FrontImpulse)
                    FrontHeelFracImpulse.append(sum(dat.FrontHeel_Filt[value-25:value+25])/FrontImpulse)
                    TotalToeFracImpulse.append(sum(dat.TotalToes_Filt[value-25:value+25])/TotImpulse)
                    TotalHeelFracImpulse.append(sum(dat.TotalHeel_Filt[value-25:value+25])/TotImpulse)
                    BackToeFracImpulse.append(sum(dat.RearToes_Filt[value-25:value+25])/BackImpulse)
                    BackHeelFracImpulse.append(sum(dat.RearHeel_Filt[value-25:value+25])/BackImpulse)
                  
                    
                    if jj < len(Toepeaks)-1:
                        # Proportion of the force on the Posterior Portion during
                        # the end of the turn (last 25%)
                        idx_max = np.argmax(dat.TotalToes_Filt[value-window:value+window])+value-window
                        idx_nextmin = np.argmin(dat.TotalToes_Filt[idx_max:Toepeaks[jj+1]])+idx_max
                        
                        idx75 = round(0.5*(idx_nextmin-idx_max))+idx_max
                        BackFootFracImpulseLate.append(sum(dat.RearTotal_Filt[idx75:idx_nextmin])/sum(dat.TotalForce_Filt[idx75:idx_nextmin]))
                    else:
                        BackFootFracImpulseLate.append(np.nan)
                        
                    if jj>0:
                        # Rate of force development
                        # Find the index of the true local maxima
                        idx_max = np.argmax(dat.TotalToes_Filt[value-window:value+window])+value-window
                        maxF = dat.TotalToes_Filt[idx_max]
                        #idx_min = np.argmin(dat.TotalToes_Filt[Toepeaks[jj-1]:idx_max])+Toepeaks[jj-1]
                        
                        
                        mins = signal.argrelmin(np.array(dat.TotalToes_Filt[idx_max-400:idx_max]))
                        
                        idx_min = []
                        for i in range(len(mins[0])):
                            if dat.TotalToes_Filt[mins[0][i]+idx_max-400]<maxF*0.1:
                                idx_min.append(mins[0][i]+idx_max-400)
                                
                        idx_min = idx_min[-1]
                        
                        idx_pk = signal.find_peaks(dat.TotalToes_Filt[idx_min:idx_max+ 25], height = 0.8*maxF)[0][0] + idx_min 
                            
                                               
                        RFD.append(np.mean(np.gradient(dat.TotalToes_Filt[idx_min:idx_pk],1/freq)))
                        RFDtime.append((idx_pk-idx_min)/freq)
                        
                        # Proportion of the force on the Anterior Portion during
                        # The start of the turn (first 25%)
                        idx25 = round(0.5*(idx_max-idx_min))+ idx_min
                        FrontFootFracImpulseEarly.append(sum(dat.TotalToes_Filt[idx_min:idx25])/sum(dat.TotalToes_Filt[idx_min:idx25]))
                        ToeTurn_ToeTS.append(dat.TotalToes_Filt[idx_min:idx_min + 300])
                        ToeTurn_HeelTS.append(dat.TotalToes_Filt[idx_min:idx_min + 300])
                    else:
                        RFD.append(np.nan)
                        RFDtime.append(np.nan)
                        FrontFootFracImpulseEarly.append(np.nan)
                        
                    
                    if jj>0 and jj < len(Toepeaks)-1:
                        ToeTotAvgForce.append(np.nanmean(dat.TotalToes_Filt[idx_min:idx_nextmin]))
                        HeelTotAvgForce.append(np.nanmean(dat.TotalHeel_Filt[idx_min:idx_nextmin]))
                    else:
                        ToeTotAvgForce.append(np.nan)
                        HeelTotAvgForce.append(np.nan)
                        
                    TurnDir.append('Toe')
                    sName.append(subName)
                    cName.append(configName)
                    trialNo.append(trialNoTmp)
                    
                    
                except:
                    print('Bad turn: ', fName, ' Toe Turn ', jj)
                    
            for jj, value in enumerate(Heelpeaks):
                # Loop through all cleaned Right Turns to calculate discrete outcome measures.
                # using right side only (left turns). No longer assuming left to right
                # transitions are constant and only using data from one side. 
                # Extend this to Left Turns as well
                # variables of interest: outside force (downhill foot peak force) higher is better
                # Outside medial peak force (higher is better)
                # Outside heel force lower is better with an average higher force on forefoot Tricky**
    
                try:
                    
                    #jj = 1
                    #value = Heelpeaks[jj]
                    if value < 100:
                        window = value
                    else:
                        window = 100
                    # Look at a 200 frame window for the true max force
                    ToeTotMaxForce.append(np.nanmax(dat.TotalToes_Filt[value-window:value+window]))
                    FrontToeMaxForce.append(np.nanmax(dat.FrontToes_Filt[value-window:value+window]))
                    FrontHeelMaxForce.append(np.nanmax(dat.FrontHeel_Filt[value-window:value+window]))
                    HeelTotMaxForce.append(np.nanmax(dat.TotalHeel_Filt[value-window:value+window]))
                    BackToeMaxForce.append(np.nanmax(dat.RearToes_Filt[value-window:value+window]))
                    BackHeelMaxForce.append(np.nanmax(dat.RearHeel_Filt[value-window:value+window]))
                    FrontTotalMaxForce.append(np.nanmax(dat.FrontTotal_Filt[value-window:value+window]))
                    BackTotalMaxForce.append(np.nanmax(dat.RearTotal_Filt[value-window:value+window]))
                    
                    
                    
                    # Examine a 0.5 second window around the peak of the turn
                    # for fraction of force  exerted on the foot of the outside (downhill) ski
                    TotImpulse = sum(dat.TotalForce_Filt[value-25:value+25])
                    FrontImpulse = sum(dat.FrontTotal_Filt[value-25:value+25])
                    BackImpulse = sum(dat.RearTotal_Filt[value-25:value+25])
                    FrontToeFracImpulse.append(sum(dat.FrontToes_Filt[value-25:value+25])/FrontImpulse)
                    FrontHeelFracImpulse.append(sum(dat.FrontHeel_Filt[value-25:value+25])/FrontImpulse)
                    TotalToeFracImpulse.append(sum(dat.TotalToes_Filt[value-25:value+25])/TotImpulse)
                    TotalHeelFracImpulse.append(sum(dat.TotalHeel_Filt[value-25:value+25])/TotImpulse)
                    BackToeFracImpulse.append(sum(dat.RearToes_Filt[value-25:value+25])/BackImpulse)
                    BackHeelFracImpulse.append(sum(dat.RearHeel_Filt[value-25:value+25])/BackImpulse)
                    
                    if jj < len(Heelpeaks)-1:
                        # Proportion of the force on the Posterior Portion during
                        # the end of the turn (last 25%)
                        idx_max = np.argmax(dat.TotalHeel_Filt[value-window:value+window])+value-window
                        idx_nextmin = np.argmin(dat.TotalHeel_Filt[idx_max:Heelpeaks[jj+1]])+idx_max
                        
                        idx75 = round(0.5*(idx_nextmin-idx_max))+idx_max
                        BackFootFracImpulseLate.append(sum(dat.RearTotal_Filt[idx75:idx_nextmin])/sum(dat.TotalForce_Filt[idx75:idx_nextmin]))
                    else:
                        BackFootFracImpulseLate.append(np.nan)
                    
                    # Rate of force development
                    if jj>0:
                        # Find the index of the true local maxima
                        idx_max = np.argmax(dat.TotalHeel_Filt[value-window:value+window])+value-window
                        #idx_min = np.argmin(dat.TotalHeel_Filt[Heelpeaks[jj-1]:idx_max])+Heelpeaks[jj-1]
                        
                        maxF = dat.TotalHeel_Filt[idx_max]
                        #idx_min = np.argmin(dat.TotalToes_Filt[Toepeaks[jj-1]:idx_max])+Toepeaks[jj-1]
                        
                        mins = signal.argrelmin(np.array(dat.TotalHeel_Filt[idx_max-400:idx_max]))
                        
                        idx_min = []
                        for i in range(len(mins[0])):
                            if dat.TotalHeel_Filt[mins[0][i]+idx_max-400]<maxF*0.2:
                                idx_min.append(mins[0][i] + idx_max-400)
                                
                        idx_min = idx_min[-1]
                                
                                                     
                        
                        
                        idx_pk = signal.find_peaks(dat.TotalHeel_Filt[idx_min:idx_max + 25], height = 0.8*maxF)[0][0] + idx_min
                                
                        RFD.append(np.mean(np.gradient(dat.TotalHeel_Filt[idx_min:idx_pk],1/freq)))
                        RFDtime.append((idx_pk-idx_min)/freq)
                        
                        # Proportion of the force on the Anterior Portion during
                        # The start of the turn (first 25%)
                        idx25 = round(0.5*(idx_max-idx_min))+idx_min
                        FrontFootFracImpulseEarly.append(sum(dat.FrontTotal_Filt[idx_min:idx25])/sum(dat.TotalForce_Filt[idx_min:idx25]))
                        
                        HeelTurn_HeelTS.append(dat.TotalHeel_Filt[idx_min:idx_min + 300])
                        HeelTurn_ToeTS.append(dat.TotalToes_Filt[idx_min:idx_min + 300])
                    else:
                        RFD.append(np.nan)
                        RFDtime.append(np.nan)
                        FrontFootFracImpulseEarly.append(np.nan)
                    
                    if jj>0 and jj < len(Heelpeaks)-1:
                        ToeTotAvgForce.append(np.nanmean(dat.TotalToes_Filt[idx_min:idx_nextmin]))
                        HeelTotAvgForce.append(np.nanmean(dat.TotalHeel_Filt[idx_min:idx_nextmin]))
                    else:
                        ToeTotAvgForce.append(np.nan)
                        HeelTotAvgForce.append(np.nan)
                        
                    

                    TurnDir.append('Heel')
                    sName.append(subName)
                    cName.append(configName)
                    trialNo.append(trialNoTmp)
                    
                except:
                    print('Bad Turn: ',fName, ' Heel Turn ', jj)
    
        ToeTurn_ToeTS = np.stack(ToeTurn_ToeTS)
        ToeTurn_HeelTS = np.stack(ToeTurn_HeelTS)
        HeelTurn_HeelTS = np.stack(HeelTurn_HeelTS)
        HeelTurn_ToeTS = np.stack(HeelTurn_ToeTS)
        
        plt.figure(subName + configName + ' Toe Turns')
        plt.plot(np.mean(ToeTurn_ToeTS, axis = 0), color = 'b', label = 'Toe Force', linewidth = 2)
        plt.legend()
        for t in range(ToeTurn_ToeTS.shape[0]):         
            plt.figure(subName + configName + ' Toe Turns')
            plt.plot(ToeTurn_ToeTS[t, :], color = 'b', linewidth = 0.2)
            #plt.plot(ToeTurn_HeelTS[t, :], color = 'r', label = 'Heel Force')
            
                
        plt.figure(subName + configName + ' Heel Turns')
        plt.plot(np.mean(HeelTurn_HeelTS, axis = 0), color = 'b', label = 'Heel Force', linewidth = 2)
        plt.legend()
        for t in range(HeelTurn_ToeTS.shape[0]):
            plt.figure(subName + configName + ' Heel Turns')
            plt.plot(HeelTurn_HeelTS[t, :], color = 'b', linewidth = 0.2)
                        
    except Exception as e: print(e)


# Create data frame with all outcome measures and export to csv. Will create a new csv if one does not exist for this dataset. 
# Otherwise will append.

outcomes = pd.DataFrame({'Subject':list(sName),'Config':list(cName),'TurnDir': list(TurnDir), 'TrialNo':list(trialNo),
                                 'ToeTotMaxForce':list(ToeTotMaxForce),
                                 'FrontToeMaxForce':list(FrontToeMaxForce),
                                 'BackToeMaxForce':list(BackToeMaxForce),
                                 'FrontHeelMaxForce':list(FrontHeelMaxForce),
                                 'FrontTotalMaxForce':list(FrontTotalMaxForce),
                                 'BackTotalMaxForce':list(BackTotalMaxForce),
                                 'HeelTotMaxForce':list(HeelTotMaxForce),
                                 'BackHeelMaxForce':list(BackHeelMaxForce),
                                 'ToeTotAvgForce':list(ToeTotAvgForce),
                                 'HeelTotAvgForce':list(HeelTotAvgForce),
                                 'FrontToeFracImpulse':list(FrontToeFracImpulse),
                                 'FrontHeelFracImpulse':list(FrontHeelFracImpulse),
                                 'TotalToeFracImpulse':list(TotalToeFracImpulse),
                                 'TotalHeelFracImpulse':list(TotalHeelFracImpulse),
                                 'BackToeFracImpulse':list(BackToeFracImpulse),
                                 'BackHeelFracImpulse':list(BackHeelFracImpulse),
                                 'FrontFootFracImpulseEarly':list(FrontFootFracImpulseEarly),
                                 'BackFootFracImpulseLate':list(BackFootFracImpulseLate),
                                 'RFD':list(RFD),
                                 'RFDtime':list(RFDtime)
                                 })

         
  
outfileName = fPath + 'CompiledResultsTest3.csv'

if save_on == 1:
    if os.path.exists(outfileName) == False:
        
        outcomes.to_csv(outfileName, mode='a', header=True, index = False)
    
    else:
        outcomes.to_csv(outfileName, mode='a', header=False, index = False) 
else:
    print('Complete')


