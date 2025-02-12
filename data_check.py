#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 12:21:48 2025

@author: ncro8394
"""
import numpy as np
import pandas as pd
import samplerate
from scipy.signal import decimate, find_peaks, periodogram, welch
from scipy.fft import fft, fftfreq
from sleepecg import detect_heartbeats, get_toy_ecg
from wonambi import Dataset
from wonambi.attr import Annotations

def gini(x, w=None):
    # Calculates the Gini coefficient:
    #   a measure of statistical dispersion calculated by comparing the actual 
    #   Lorenz curve to the diagonal line of equality.
    x = np.asarray(x)
    if w is not None:
        w = np.asarray(w)
        sorted_indices = np.argsort(x)
        sorted_x = x[sorted_indices]
        sorted_w = w[sorted_indices]
        # Force float dtype to avoid overflows
        cumw = np.cumsum(sorted_w, dtype=float)
        cumxw = np.cumsum(sorted_x * sorted_w, dtype=float)
        return (np.sum(cumxw[1:] * cumw[:-1] - cumxw[:-1] * cumw[1:]) / 
                (cumxw[-1] * cumw[-1]))
    else:
        sorted_x = np.sort(x)
        n = len(x)
        cumx = np.cumsum(sorted_x, dtype=float)
        # The above formula, with all weights equal to 1 simplifies to:
        return (n + 1 - 2 * np.sum(cumx) / cumx[-1]) / n

def moving_average_time(a, s_freq, win=5):
    n = win*s_freq*60
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

def moving_average(a, win=100):
    n = win
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

file = '/Users/ncro8394/Documents/projects/seapipe/DATA/sub-IN003/ses-V1/eeg/sub-IN003_ses-V1_eeg.edf'

d = Dataset(file)




chans = d.header['chan_name']


## QC chosen EEG channels
# Frontal
f_chans = [x for x in chans if 'F3' in x or 'F4' in x or 'Fz' in x]
fch_dat = d.read_data(chan=f_chans)
quals = pd.DataFrame(data=None, dtype = float, 
                      columns=['integral',
                              'stdev',
                              'time_not_in_low_range',
                              'time_not_flatline',
                              '99th_quantile',
                              'gini_coeff'])
for i, name in enumerate(f_chans):
    print(f'Checking quality of {name}')
    s_freq = d.header['orig']['n_samples_per_record'][chans.index(name)]
    e = moving_average_time(abs(fch_dat.data[0][i]),s_freq,5)
    quals.loc[i,'integral'] = sum(e) #intragral (AUC)
    quals.loc[i,'stdev'] = np.std(e) #standard deviation 
    quals.loc[i,'time_not_in_low_range'] = 1 - sum(np.logical_and(e>=0, e<=0.5))/len(e)
    quals.loc[i,'time_not_flatline'] = 1 - (sum(np.diff(e)==0)/len(np.diff(e)))
    quals.loc[i,'99th_quantile'] = np.quantile(e,0.99)
    quals.loc[i,'gini_coeff'] = gini(e)

winners = pd.DataFrame([quals[x].rank() for x in quals.columns]).T
idx = winners.sum(axis=1).idxmax()
f_name = f_chans[idx]
f_dat = fch_dat.data[0][idx]
print(f'Best Frontal channel based on auto-QC is: {f_name}')

# Occipital
o_chans = [x for x in chans if 'O1' in x or 'O2' in x or 'Oz' in x]
och_dat = d.read_data(chan=o_chans)
quals = pd.DataFrame(data=None, dtype = float, 
                     columns=['integral',
                              'stdev',
                              'time_not_in_low_range',
                              'time_not_flatline',
                              '99th_quantile',
                              'gini_coeff'])
for i, name in enumerate(o_chans):
    print(f'Checking quality of {name}')
    s_freq = d.header['orig']['n_samples_per_record'][chans.index(name)]
    e = moving_average_time(abs(och_dat.data[0][i]),s_freq,5)
    quals.loc[i,'integral'] = sum(e) #intragral (AUC)
    quals.loc[i,'stdev'] = np.std(e) #standard deviation 
    quals.loc[i,'time_not_in_low_range'] = 1 - sum(np.logical_and(e>=0, e<=0.5))/len(e)
    quals.loc[i,'time_not_flatline'] = 1 - (sum(np.diff(e)==0)/len(np.diff(e)))
    quals.loc[i,'99th_quantile'] = np.quantile(e,0.99)
    quals.loc[i,'gini_coeff'] = gini(e)

winners = pd.DataFrame([quals[x].rank() for x in quals.columns]).T
idx = winners.sum(axis=1).idxmax()
o_name = o_chans[idx]
o_dat = och_dat.data[0][idx]
print(f'Best Occipital channel based on auto-QC is: {o_name}')

    

## Choose EMG channel
emg_chans = [x for x in chans if 'EMG' in x or 'emg' in x or 'chin' in x]
emg_dat = d.read_data(chan=emg_chans)
quals = pd.DataFrame(data=None, dtype = float, 
                     columns=['integral',
                              'stdev',
                              'time_not_in_low_range',
                              'time_not_flatline',
                              '99th_quantile',
                              'gini_coeff'])

for i, name in enumerate(emg_chans):
    print(f'Checking quality of {name}')
    s_freq = d.header['orig']['n_samples_per_record'][chans.index(name)]
    e = moving_average_time(abs(emg_dat.data[0][i]),s_freq,5)
    quals.loc[i,'integral'] = sum(e) #intragral (AUC)
    quals.loc[i,'stdev'] = np.std(e) #standard deviation 
    quals.loc[i,'time_not_in_low_range'] = 1 - sum(np.logical_and(e>=0, e<=0.5))/len(e)
    quals.loc[i,'time_not_flatline'] = 1 - (sum(np.diff(e)==0)/len(np.diff(e)))
    quals.loc[i,'99th_quantile'] = np.quantile(e,0.99)
    quals.loc[i,'gini_coeff'] = gini(e)

winners = pd.DataFrame([quals[x].rank() for x in quals.columns]).T
idx = winners.sum(axis=1).idxmax()
emg_name = emg_chans[idx]
emg_dat = emg_dat.data[0][idx]
print(f'Best EMG channel based on auto-QC is: {emg_name}')
    
    
## Choose ECG channel
ecg_chans = [x for x in chans if 'ECG' in x or 'ecg' in x]
ecg_dat = d.read_data(chan=ecg_chans)
quals = pd.DataFrame(data=None, columns=['hr_qual',
                                         'stdev',
                                         'time_not_flatline',
                                         'PSD_qual'])

for i, name in enumerate(ecg_chans):
    print(f'Checking quality of {name}')
    s_freq = d.header['orig']['n_samples_per_record'][chans.index(name)]
    e = ecg_dat.data[0][i]
    
    # Check BPM
    beats = detect_heartbeats(e, s_freq)
    hr = len(beats)/(len(e)/(s_freq*60))
    if np.logical_and(hr>40, hr<120):
        print(f'Heart rate ({round(hr)}bpm) seems physiologically plausible. Continuing...')
        quals.loc[i,'hr_qual'] = 2
    elif np.logical_and(hr>120, hr<200):
        quals.loc[i,'hr_qual'] = 1
        print(f'WARNING: Heart rate ({round(hr)}bpm) is physiologically plausible but high. Check...')
    else:
        quals.loc[i,'hr_qual'] = 0
        print(f'WARNING: Heart rate ({round(hr)}bpm) is NOT physiologically plausible. Signal quality is likely poor...')
    
    # Check Std Dev.
    quals.loc[i,'stdev'] = np.std(abs(e)) #standard deviation 
    
    # Check for flatlines
    quals.loc[i,'time_not_flatline'] = 1 - (sum(np.diff(e)==0)/len(np.diff(e)))
    
    # Check PSD        
    f_p, Pxx_spec = periodogram(e, s_freq)
    Pxx_spec_sm = np.concatenate((np.zeros(2000),
                                  moving_average(Pxx_spec, 4000),
                                  np.zeros(1999)))
    top = int((len(e)/s_freq)*20) # corresponds to 20Hz
    peaks = find_peaks(Pxx_spec_sm[:top], prominence=800)[0]
    
    if len(peaks) > 4:
        quals.loc[i,'PSD_qual'] = 2
    elif len(peaks) > 2:
        quals.loc[i,'PSD_qual'] = 1
    else:
        print(f'WARNING: PSD of {name} looks unusual..')
        quals.loc[i,'PSD_qual'] = 0

winners = pd.DataFrame([quals[x].rank() for x in quals.columns]).T
idx = winners.sum(axis=1).idxmax()
ecg_name = ecg_chans[idx]
ecg_dat = ecg_dat.data[0][idx]
print(f'Best ECG channel based on auto-QC is: {ecg_name}')


print('Resampling data..')
f_dat = samplerate.resample(f_dat, 200/s_freq, 'sinc_fastest')
o_dat = samplerate.resample(o_dat, 200/s_freq, 'sinc_fastest')
emg_dat = samplerate.resample(emg_dat, 200/s_freq, 'sinc_fastest')
ecg_dat = samplerate.resample(ecg_dat, 200/s_freq, 'sinc_fastest')

mat = np.vstack((f_dat, o_dat, emg_dat, ecg_dat))

# Resize array to 16m and save
ideal_len = 4096*4096
mat = np.concatenate((mat, np.zeros( (mat.shape[0], ideal_len-mat.shape[1]) )), axis=1)
np.save('/Users/ncro8394/DeepSleep/data/feature_16m/sub-IN003.npy' , mat)


## Format labels 
annot = Annotations('/Users/ncro8394/Documents/projects/seapipe/derivatives/staging_manual/sub-IN003/ses-V1/sub-IN003_ses-V1_eeg.xml')
evts = annot.get_events(name = 'Arou')
labels = np.zeros(mat.shape[1])
for e in evts:
    start = int(e['start']*s_freq)
    end = int(e['end']*s_freq)
    labels[start:end] = 1
    
# Resize array to 16m and save    
labels = np.reshape(labels, (1, len(labels)))    
labels = np.concatenate((labels, np.zeros((labels.shape[0], ideal_len-labels.shape[1])) ), axis=1)
np.save('/Users/ncro8394/DeepSleep/data/label_16m/sub-IN003.npy' , labels)



