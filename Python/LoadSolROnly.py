# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 13:14:07 2021

@author: Daniel.Feeney
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

def calcVLR(force, startVal, lengthFwd):
    # function to calculate VLR from 80 and 20% of the max value observed in the first n
    # indices (n defined by lengthFwd). 
    maxF = np.max(force[startVal:startVal+lengthFwd])
    eightyPctMax = 0.8 * maxF
    twentyPctMax = 0.2 * maxF
    # find indices of 80 and 20 and calc loading rate as diff in force / diff in time (N/s)
    eightyIndex = next(x for x, val in enumerate(force[startVal:startVal+lengthFwd]) 
                      if val > eightyPctMax) 
    twentyIndex = next(x for x, val in enumerate(force[startVal:startVal+lengthFwd]) 
                      if val > twentyPctMax) 
    VLR = ((eightyPctMax - twentyPctMax) / ((eightyIndex/1000) - (twentyIndex/1000)))
    return(VLR)

# Define constants and options
fThresh = 85; #below this value will be set to 0.
minStepLen = 10; #minimal step length
writeData = 0; #will write to spreadsheet if 1 entered
desiredStepLength = 25; #length to look forward after initial contact

# Read in file and add names
fPath = 'C:/Users/Daniel.Feeney/Dropbox (Boa)/EndurancePerformance/Altra_MontBlanc_Jan2021/RawData/'
entries = os.listdir(fPath)

for file in entries[5:7]:
    try:
                
        fName = file
        
        dat = pd.read_csv(fPath+fName,sep='\s+', skiprows = 3, header = 0)
        dat = dat.drop(dat.columns[[5,6,7,8]], axis=1)  
        dat.columns = ['Time', 'RightLateral','RightMedial','RightHeel','RightTotal']
        
        subName = fName.split(sep = "_")[0]
        Test = fName.split(sep = "_")[1]
        Config = fName.split(sep = "_")[2].split(sep=".")[0]
        
        #### Do the same thing for the right side #### 
        ##### Filter force below threshold to 0 #####
        RForce = dat.RightTotal
        RForce[RForce<fThresh] = 0
        
        # delimit steps on right side
        ric = []
        count = 1;
        for step in range(len(RForce)-1):
            if RForce[step] == 0 and RForce[step + 1] >= fThresh:
                ric.append(step)
                count = count + 1
        #left right off
        rto = []
        count = 1;
        for step in range(len(RForce)-1):
            if RForce[step] >= fThresh and RForce[step + 1] == 0:
                rto.append(step + 1)
                count = count + 1
        
        # ## plotting left side
        # fig, ax = plt.subplots(4)
        # for i in range(noSteps):
        #     ax[0].plot(RightMat[i,:])
        #     ax[1].plot(RHeelMat[i,:])
        #     ax[2].plot(RLatMat[i,:])
        #     ax[3].plot(RMedMat[i,:])
            
        # ax[0].set_title('Total Right Force')
        # ax[1].set_title('Right Heel Force')
        # ax[2].set_title('Right Lateral Force')
        # ax[3].set_title('Right Medial Force')
        
        #%%
        
        # Right side
        MaxR = []
        totImpulseR = []
        pkHeelR = []
        heelImpulseR = []
        pkLatR = []
        latImpulseR = []
        pkMedR = []
        medImpulseR = []
        stanceTimeR = []
        rateTotR = []  
        nameR = []
        configR = []
         
        for step in ric:
            try:
                rateTotR.append(calcVLR(dat.RightTotal, step, desiredStepLength))
            except:
                rateTotR.append(0)    
            #Right
            MaxR.append(np.max(dat.RightTotal[step:step+desiredStepLength]))
            totImpulseR.append(np.sum(dat.RightTotal[step:step+desiredStepLength]))
            pkHeelR.append(np.max(dat.RightHeel[step:step+desiredStepLength]))
            heelImpulseR.append(np.sum(dat.RightHeel[step:step+desiredStepLength]))
            pkLatR.append(np.max(dat.RightLateral[step:step+desiredStepLength]))
            latImpulseR.append(np.sum(dat.RightLateral[step:step+desiredStepLength]))
            pkMedR.append(np.max(dat.RightMedial[step:step+desiredStepLength]))
            medImpulseR.append(np.sum(dat.RightMedial[step:step+desiredStepLength]))
            nameR.append(subName)
            configR.append(Config)
        
        rightDat = pd.DataFrame({'Sub':list(nameR), 'Config': list(configR),
                      'VLR':list(rateTotR),'MaxF': list(MaxR),'pkHeel': list(pkHeelR), 'HeelImpulse': list(heelImpulseR),
                      'PkLat': list(pkLatR), 'LatImp': list(latImpulseR), 'PkMed': list(pkMedR),
                      'MedImp': list(medImpulseR)})
        
        rightDat.to_csv('C:/Users/Daniel.Feeney/Dropbox (Boa)/EndurancePerformance/Altra_MontBlanc_Jan2021/RawData//SummarizedResults3.csv', mode='a', header=False)

    except:
        print(file)



