# -*- coding: utf-8 -*-
"""
Created on Mon Aug 30 13:28:44 2021
Script to process MVA files from cycling pilot test

@author: Daniel.Feeney
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.misc import electrocardiogram
from scipy.signal import find_peaks
import seaborn as sns

fwdLook = 30

# list of functions 
# finding landings on the force plate once the filtered force exceeds the force threshold
def findLandings(force, fThreshold):
    ric = []
    for step in range(len(force)-1):
        if force[step] == 0 and force[step + 1] >= fThreshold:
            ric.append(step)
    return ric

#Find takeoff from FP when force goes from above thresh to 0
def findTakeoffs(force, fThreshold):
    rto = []
    for step in range(len(force)-1):
        if force[step] >= fThreshold and force[step + 1] == 0:
            rto.append(step + 1)
    return rto

def trimTakeoffs(landings, takeoffs):
    if takeoffs[0] < landings[0]:
        del(takeoffs[0])
    return(takeoffs)

def trimLandings(landings, trimmedTakeoffs):
    if landings[len(landings)-1] > trimmedTakeoffs[len(trimmedTakeoffs)-1]:
        del(landings[-1])
    return(landings)

# Read in files
# only read .asc files for this work
fPath = 'C:\\Users\\daniel.feeney\\Boa Technology Inc\\PFL - General\\Cycling2021\\EH_CyclingPilot_2021\\Pressures\\'
fileExt = r".mva"
entries = [fName for fName in os.listdir(fPath) if fName.endswith(fileExt)]
#fThresh = 200; #below this value will be set to 0.
sub = []
config = []
condition = []
trial = []
initialPct = []
peakPct = []
endPct = []

for file in entries:
        try:
            fName = file #Load one file at a time
            dat = pd.read_csv(fPath+fName, sep='\t', skiprows = 11, header = 0)
            dat.rename(columns={'time[secs]    ': 'Time', 'force[N]   ': 'LForce',
                                'max pressure ':'LMaxP', 'mean pressure ':'LMeanP',
                                '+-% of mean ':'LPctMean','force[N]   .1':'RForce',
                                'max pressure .1':'RMaxP', 'mean pressure .1':'RMeanP',
                                '+-% of mean .1':'RPctMean'}, inplace=True)
            
            landings = findLandings(dat.LForce, 1)
            takeoffs = findTakeoffs(dat.LForce,1)
            trimmedTakeoffs = trimTakeoffs(landings, takeoffs)
            trimmedLandings = trimLandings(landings, trimmedTakeoffs)
            
            for countVar, landing in enumerate(trimmedLandings):
                
                tmpForce = dat.LForce[landing : landing + fwdLook]
                tmpPk = max(tmpForce)
                timePk = list(tmpForce).index(tmpPk) #indx of max force applied during that pedal stroke
                
                initialPct.append( dat.LPctMean[landing+1] )
                peakPct.append( dat.LPctMean[landing + timePk] )
                try:
                    endPct.append( dat.LPctMean[trimmedTakeoffs[countVar]-1] )
                except:
                    endPct.append(0)
            
                sub.append( fName.split('_')[0] )
                config.append( fName.split('_')[1] )
                condition.append( fName.split('_')[2] )
                trial.append( fName.split('_')[3].split('.')[0] )
        except:
            print(file)
        
        
outcomes = pd.DataFrame({ 'Subject':list(sub),'config':list(config), 'condition': list(condition), 'trial': list(trial),
                 'initialPct': list(initialPct), 'peakPct': list(peakPct),'endPct': list(endPct) })
        
        
        