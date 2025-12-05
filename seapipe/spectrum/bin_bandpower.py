#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 17 09:51:10 2025

@author: ncro8394
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

file = '/Users/ncro8394/Documents/projects/healthy/test_dataset/derivatives/bandpower_timecourse/bandpower_timeseries_0.5-4Hz_sliding/sub-HC001/ses-V1/sub-HC001_ses-V1_Fz_NREM2_sliding_0.5-4Hz.csv'

bp_df = pd.read_csv(file)



bp_ts = bp_df['band_power']

threshold = np.std(bp_ts)*10

bp_ts[bp_ts>threshold] = np.nan

plt.figure()
plt.plot(bp_ts)

n_bins = 8
bin_size = len(bp_ts) // n_bins

# Trim remainder so length divides evenly
trim = len(bp_ts) % n_bins
if trim:
    bp_ts = bp_ts[:-trim]

binned = np.nanmean(bp_ts.values.reshape(n_bins, bin_size), axis=1)



x = list(range(1,n_bins+1))
plt.bar(x, binned)




file2 = '/Users/ncro8394/Documents/projects/healthy/test_dataset/derivatives/clusterfluc/sub-HC001/ses-V0/sub-HC001_ses-V0_Fz_NREM2_SWA_fluctuations_timeseries.csv'

lf_df = pd.read_csv(file2).iloc[0]

bin_size = len(lf_df) // n_bins
# Trim remainder so length divides evenly
trim = len(lf_df) % n_bins
if trim:
    lf_df = lf_df[:-trim]

binned2 = np.nanmean(lf_df.values.reshape(n_bins, bin_size), axis=1)

