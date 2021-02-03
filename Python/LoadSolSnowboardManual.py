# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 16:19:30 2020

@author: Daniel.Feeney
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import ctypes  

# Read in files
# only read .asc files for this work
fPath = 'C:/Users/Daniel.Feeney/Dropbox (Boa)/Snow Protocol/SnowboardProtocol/'
entries = os.listdir(fPath)
fName = entries[0]
dat = pd.read_csv(fPath+fName,sep='         ', skiprows = 3, header = 0, index_col = False)
#dat.columns = ['Time', 'LeftHeel', 'LeftMedial','LeftLateral','Total']
dat.columns = ['Time', 'LHeel', 'LMedial','LLateral','LTotal', 'Time2', 'RLateral','RMedial','RHeel','RTotal']

dat['LToes'] = dat.LMedial + dat.LLateral
dat['RToes'] = dat.RMedial + dat.RLateral

fig, ax = plt.subplots()
ax.plot(dat.LTotal, label = 'Left Total Force')
ax.plot(dat.RTotal, label = 'Right Total Force')
fig.legend()
print('Select start and end of analysis trial')
pts = np.asarray(plt.ginput(2, timeout=15))
plt.close()
# downselect the region of the dataframe you'd like
dat = dat.iloc[int(np.floor(pts[0,0])) : int(np.floor(pts[1,0])),:]

#fig, ax = plt.subplots()
#ax.plot(dat.LToes, label='Left Toes')
#ax.plot(dat.LHeel, label='Left Heel')
#fig.legend()

#### Initial analysis plans: segment heel and toe turns, calcualte smoothness
## as the CV or SD of force during the turn, calcualte turn time, calculate symmetry

# Find indices of toe and heel turns but leave original data untouched 
heelThresh = 150 #below this value will be set to 0 temporarily to find indices to start/end turns
toeThresh = 75

tmpToes = np.array(dat.LToes)
tmpHeel = np.array(dat.LHeel)

tmpHeel[tmpHeel < heelThresh] = 0
tmpToes[tmpToes < toeThresh] = 0

# Select toe turn starts
plt.plot(tmpToes, label='Toes')
plt.plot(tmpHeel, label='Heel')
plt.legend()
print('Select start of Toe Turns')
toestart = np.asarray(plt.ginput(5, timeout=-1))
realToeStart = [int(toestart[0,0]), int(toestart[1,0]), (int(toestart[2,0])), int(toestart[3,0]),int(toestart[4,0])]
plt.close()

# Select heel turn starts
plt.plot(tmpToes, label='Toes')
plt.plot(tmpHeel, label='Heel')
plt.legend()
print('Select start of heel turns')
heelstart = np.asarray(plt.ginput(5, timeout=-1))
realHeelStart = [int(heelstart[0,0]), int(heelstart[1,0]), (int(heelstart[2,0])), int(heelstart[3,0]),int(heelstart[4,0])]
plt.close()   

###### After finding starts of turns, find avg, SD, CV, etc. for each turn ####
dat = dat.reset_index()
turnToPlot = 3
fwdLook = 200

fig, ax = plt.subplots(2)
fig.suptitle('Toe Turn')
ax[0].plot(dat.LToes[realToeStart[turnToPlot]:realToeStart[turnToPlot]+fwdLook], label='L Toes')
ax[0].plot(dat.RToes[realToeStart[turnToPlot]:realToeStart[turnToPlot]+fwdLook], label = 'R Toes')
ax[0].set_title('Toes')
ax[0].set_ylim([0,800])
fig.legend()
ax[1].plot(dat.LHeel[realToeStart[turnToPlot]:realToeStart[turnToPlot]+fwdLook], 'tab:green', label = 'L Heel')
ax[1].plot(dat.RHeel[realToeStart[turnToPlot]:realToeStart[turnToPlot]+fwdLook], 'tab:red', label = 'R Heel')
ax[1].set_title('Heel')
ax[1].set_ylim([0,800])
fig.legend()

fig, ax = plt.subplots(2)
fig.suptitle('Heel Turn')
ax[0].plot(dat.LToes[realHeelStart[turnToPlot]:realHeelStart[turnToPlot]+fwdLook], label='L Toes')
ax[0].plot(dat.RToes[realHeelStart[turnToPlot]:realHeelStart[turnToPlot]+fwdLook], label = 'R Toes')
ax[0].set_ylim([0,800])
fig.legend()
ax[1].plot(dat.LHeel[realHeelStart[turnToPlot]:realHeelStart[turnToPlot]+fwdLook], 'tab:green', label = 'L Heel')
ax[1].plot(dat.RHeel[realHeelStart[turnToPlot]:realHeelStart[turnToPlot]+fwdLook], 'tab:red', label = 'R Heel')
ax[1].set_ylim([0,800])
fig.legend()


###### Extract variables from each turn initiation ######
maxF = [np.max(dat.LToes[toeTurnStart:toeTurnStart+100]) for toeTurnStart in realToeStart]
maxRFDup = [np.max(dat.LToes[toeTurnStart:toeTurnStart+100].diff()) for toeTurnStart in realToeStart]
maxRFDdn = [np.min(dat.LToes[toeTurnStart:toeTurnStart+100].diff()) for toeTurnStart in realToeStart]
timeToPeak = [ list(dat.LToes[toeTurnStart:toeTurnStart+100]).index(max(dat.LToes[toeTurnStart:toeTurnStart+100])) for toeTurnStart in realToeStart ]

list(dat.LToes[realToeStart[1]:realToeStart[1]+100]).index(max(dat.LToes[realToeStart[1]:realToeStart[1]+100]))

plt.plot(dat.LToes[realToeStart[1]:realToeStart[1]+100])

avgF2 = [movAvgForce(forceZ, landing, landing+100, 10) for toeTurnStart in realToeStart]