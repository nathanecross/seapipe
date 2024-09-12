#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 31 16:17:18 2024

@author: m.ela
"""

#%%
" Opening a recording in Wonambi."

# 1. Annotations file

# a. load function
from wonambi.attr import Annotations 
import shutil
xml = '/Users/m.ela/Desktop/DATA_temp/QC_ROCK_temp/ROCK25/STIM1/ROCK25_E_STIM1_scores_LB_AP.xml'

new_xml = '/Users/m.ela/Desktop/DATA_temp/QC_ROCK_temp/OUT/marker/ROCK25_E_STIM1_scores_test.xml'

shutil.copy(xml, new_xml)

# b. load annotations file
annot = Annotations(new_xml)

events = annot.get_events()  ## get all events in a list
#events = annot.get_events(name = 'Staresina2015', chan = 'Cz') #get event named X on channel X


# 2. EDF
# a. load function
from wonambi import Dataset
edf_file = '/Users/m.ela/Desktop/DATA_temp/QC_ROCK_temp/ROCK25/STIM1/ROCK25_E_STIM1.edf'

# b. Load edf in Wonambi
dset = Dataset(edf_file)
#header = dset.header #to read header with info on .edf
s_freq = dset.header['s_freq']

# c. VERSION1: Read data 
    # i. (between specific timepoints)
data_ref = dset.read_data(chan=['AUX'], begtime=0, endtime=3000) #extract data from channel of interest from time x to timey
# it will create de Chan * Time datatype with a weird array within array 
#dir(data_ref) #to see the header 
#data_ref.data #to see the actually data (512pt per second)
data = data_ref.data[0][0] # 2 brackes to have data in correct shape (that can be visualize)
    #ii. visualize raw
import matplotlib.pyplot as plt
plt.plot(data)



#d. VERSION 2: load in segments
from wonambi.trans import fetch
cat = (1,1,1,1)
chans = ['AUX']
ref_chan = ['M1','M2']
stage = ['NREM2','NREM3']
chan_full = ['AUX (aux)']
evt_type = None
segments = fetch(dset, annot, cat=cat, chan_full=chan_full, 
                 cycle=None, evt_type=evt_type, stage=stage,
                 buffer=3)
    
segments.read_data(chan=chans) # this makes segments into a list of segments

seg = segments[0] # call the first segment, which is a dictionary of metadata and data
seg['data'].data[0][0] # Just like in VERSION1, # 2 brackes to have data in correct shape (that can be visualize)



# 3. Filtering
from numpy import roll
from scipy.signal import firwin, kaiserord, lfilter

#dat = seg['data'].data[0][0].T*-1 # in case you want to invert the signal
#dat = seg['data'].data[0][0]
dat = data

# Filter (w FIR)
bandwidth = (0.25, 0.5)
dat_orig = dat
nyq_rate = s_freq / 2.0
width = 5.0/nyq_rate
N, beta = kaiserord(60.0, width)
cutoff = bandwidth[1]       # for bandpass, if you want lowpass, set cutoff = bandwidth[1]
b = firwin(N, cutoff/nyq_rate, window=('kaiser', beta))
fir = lfilter(b, 1.0, dat_orig)
delay = int(-0.5 * (N-1))
dat_filt = roll(fir, delay)


# 3. Visualise data
plt.plot(dat)
plt.plot(dat_filt)


from scipy.signal import find_peaks

peaks, _ = find_peaks(dat_filt, height=550)


evt_name = 'bed'

evts = []
for i, pk in enumerate(peaks):
    if pk - peaks[i-1] > s_freq:
        pk_sec = pk/s_freq
        evt = {'chan':['AUX (aux)'],
               'cycle':'',
               'end': pk_sec+0.1,
               'name': evt_name,
               'quality': 'Good',
               'stage': '',
               'start': pk_sec-0.1
               }
        evts.append(evt)

from wonambi.graphoelement import Graphoelement

grapho = Graphoelement()
grapho.events = evts
grapho.to_annot(annot, evt_name)