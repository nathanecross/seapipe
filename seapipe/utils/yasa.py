#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 21:23:38 2024

@author: ncro8394
"""

import mne
import yasa
# Load an EDF file using MNE
path = '/Users/ncro8394/Documents/projects/seapipe/DATA/sub-IN001/ses-V1/eeg'
raw = mne.io.read_raw_edf(f"{path}/sub-IN001_ses-V1_eeg.edf", preload=True)
# Initialize the sleep staging instance
sls = yasa.SleepStaging(raw, eeg_name='C3:M1', eog_name='EOGl:M2',
                        emg_name="EMG")
# Get the predicted sleep stages
hypno = sls.predict()
# Get the predicted probabilities
proba = sls.predict_proba()
# Get the confidence
confidence = proba.max(axis=1)
# Plot the predicted probabilities
sls.plot_predict_proba()

key = {'Artefact':8,
       'W': 6,
       'N1': 3,
       'N2': 2,
       'N3': 1,
       'R': 4}

hypno_num = [6]
for x in hypno:
    hypno_num.append(key[x])

from wonambi.attr import Annotations
import matplotlib.pyplot as plt
from numpy import zeros

annot = Annotations('/Users/ncro8394/Documents/projects/seapipe/OUT/staging/sub-IN001/ses-V1/sub-IN001_ses-V1_eeg.xml')
hypno_gs = [x['stage'] for x in annot.get_epochs()]


key = {'Artefact':8,
       'Wake': 6,
       'NREM1': 3,
       'NREM2': 2,
       'NREM3': 1,
       'REM': 4}

hypno_gs_num = []
for x in hypno_gs:
    hypno_gs_num.append(key[x])
    
    
agreement = zeros((len(hypno_num)))   
for i,mx in enumerate(hypno_num):
    print(mx, hypno_gs_num[i])
    if mx == hypno_gs_num[i]:
        agreement[i] = 1
        
print((sum(agreement)/len(agreement))*100)    
    