# -*- coding: utf-8 -*-
"""
Created on Mon May  3 12:14:38 2021

@author: Kate.Harrison
"""

#load packages, including spm1d

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import spm1d

# set data directory
fPath = 'C:/Users/kate.harrison/Dropbox (Boa)/EndurancePerformance/TNF_Scrambler_Apr_21/Novel_Data/Running/'
fileExt = r".mva"
entries = [fName for fName in os.listdir(fPath) if fName.endswith(fileExt)]


#functions

fThresh = 0

def findLandings(force):
    lic = []
    for step in range(len(force)-1):
        if force[step] == 0 and force[step + 1] >= fThresh:
            lic.append(step)
    return lic

#Find takeoff from FP when force goes from above thresh to 0
def findTakeoffs(force):
    lto = []
    for step in range(len(force)-1):
        if force[step] >= fThresh and force[step + 1] == 0:
            lto.append(step + 1)
    return lto

# initialize data frames

spmDF = []

#import data

for file in entries[0]:
    try:
#parse data by HS and TO

# normalize to 100 time points

# average for each subject
#restructure data

# run spm test