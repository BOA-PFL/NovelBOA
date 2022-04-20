# -*- coding: utf-8 -*-
"""
Created on Wed Jun 16 15:53:52 2021

@author: Kate.Harrison
"""

######COP analysis

# packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os


fPath = 'C:/Users/kate.harrison/Boa Technology Inc/PFL - Documents/General/PowerPerformance/NBJapan_Sept2021/Novel/'
fileExt = r".mva"
entries = [fName for fName in os.listdir(fPath) if fName.endswith(fileExt)]


Subject = []
Config = []
targetPeakFz = []


for fName in entries:

    try:
        
        
        subName = fName.split(sep = "_")[0]
        ConfigTmp = fName.split(sep="_")[2]

        dat = pd.read_csv(fPath+fName,sep='\t', skiprows = 8, header = 'infer')

        dat.columns = ['Time', 'RearHeel_Force', 'RearHeel_MaxP', 'RearHeel_MeanP', 'RearHeel_Pct', 
                       'RearMedial_Force', 'RearMedial_MaxP', 'RearMedial_MeanP', 'RearMedial_Pct',
                       'RearLateral_Force', 'RearLateral_MaxP', 'Rearlateral_MeanP', 'RearLateral_Pct',
                       'RearToe_Force', 'RearToe_MaxP', 'RearToe_MeanP', 'RearToe_Pct',
                       'FrontHeel_Force', 'FrontHeel_MaxP', 'FrontHeel_MeanP', 'FrontHeel_Pct', 
                       'FrontMedial_Force', 'FrontMedial_MaxP', 'FrontMedial_MeanP', 'FrontMedial_Pct',
                       'FrontLateral_Force', 'FrontLateral_MaxP', 'Frontlateral_MeanP', 'FrontLateral_Pct',
                       'FrontToe_Force', 'FrontToe_MaxP', 'FrontToe_MeanP', 'FrontToe_Pct'
                       ]
        
        # 2. Find midpoint of COPx (anterior-posterioe) by weighting resultant forces. 

        GRF_Rear = (dat.RearHeel_Force + dat.RearMedial_Force + dat.RearLateral_Force + dat.RearToe_Force)
        GRF_Target = (dat.FrontHeel_Force + dat.FrontMedial_Force + dat.FrontLateral_Force + dat.FrontToe_Force)
        GRF_Tot = GRF_Rear + GRF_Target
        bc = np.argmax(GRF_Tot)

        start = bc - 5

        finish = bc
        # Find GRF peak and 2ms before and after to ID period of interest

        

        #5. Calculate variables
        Subject.append(subName)
        Config.append(ConfigTmp)
        targetPeakFz.append(max(GRF_Target[start:finish+1]))
        
        
    except:
        print(fName)
    
    
outcomes = pd.DataFrame({'Subject':list(Subject),'Config':list(Config), 'targetPeakFz':list(targetPeakFz), 
                                         })

outFileName = fPath + 'CompiledGolfForceData.csv'

outcomes.to_csv(outFileName, index = False)
