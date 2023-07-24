# -*- coding: utf-8 -*-
"""
Sleep Events Analysis pipeline (SEA)
to perform spectral analyses and brain oscillation detection (spindle, SO, CFC)

@author: Nathan Cross
This script calls seapipe (Nathan Cross) and wonambi (Jordan O'Byrne) functions

Attention: The data must be organized in BIDS format 

"""
from numpy import asarray
from .stats.sleepstats import sleepstats
from .spindle.whales import whale_it, whale_farm
from .sw.slowwave import slow_it
from .cfc.mean_amps import (pac_it, pac_it_2, generate_adap_bands, 
                           cfc_grouplevel, watson_williams)
from .cfc.synchrony import event_sync
from .cfc.plots import (plot_prefphase, plot_mean_amps, plot_prefphase_group, 
                        plot_meanamps_group)


#%% 
"""  ################################################################ """
"""           0.       DEFINE GENERAL INPUTS                          """
"""  ################################################################ """

####### 1) Define root directory (where is the data)

## Decomment the one you need
#tu run from laptop (spyder)
#root_dir = r'/Users/m.ela/Desktop/DATA_temp/Workshop/' #AP laptop
root_dir = '/Users/nathancross/Library/CloudStorage/GoogleDrive-nathan.cross.90@gmail.com/My Drive/Projects/test/' #NC laptop
#root_dir = r'/Volumes/Research/Projects/90000100/MISC/AURORE/HAMAC_poor/' #on server 

#or to run on cluster (user-defined)
#root_dir = sys.argv[1] #/NAS/Projects/90000100/MISC/APNC/SPINDLE/
#part = [sys.argv[2]] #if need to run to every subject - comment for whale_farm


####### 2) Define the basic info of your data

rater=None #rater the to analyse: if one rater put None or if several rater on the .xml put name ex:'AP'

s_freq = 512 #sampling rate (ex: 512Hz)





####### 3) Define the subject you want to analyze

part =  'all' #'all' or ['SP006'] #comment if running sleepstat, whale_it, slow_it, Cfc
visit = 'all' #'all' or ['V1'] 



#%%
"""  ######################  Sleep Stats ############################ """
"""             1.     Calculate Sleep Statistics                     """
"""  ################################################################ """
#Extract macro-sleep statistic and arousal stat
#.edf and .xml in in_dir  

####### 1.a Set general parameters
in_dir = root_dir + 'DATA/' # folder inside root_dir where BIDS data is stored
xml_dir = root_dir + 'DATA/'
out_dir = root_dir + 'OUT' # folder inside root_dir for output files
lights = 'light_times_hamac.xlsx' #xlsx file inside root_dir with lights on and off times 
                                  ## attention: needs one column for ID, one column for visit, 
                                  ## one column for loff and lon
times = root_dir + lights
cycle_idx = None #  [1,2,3] or enter None to run analyses ignoring cycles
stage = ['NREM1', 'NREM2', 'NREM3', 'REM', 'Wake'] # stages: ['NREM1', 'NREM2', 'NREM3', 'REM', 'Wake']
evt_type = ['Arousal'] #name of arousal event on .xml ex: ['Arou']


# ####### 1.b Call wozombie function: sleepstats
sleepstats(in_dir,xml_dir,out_dir,rater,evt_type,stage,times,part=part,visit=visit)



#%%
"""  #######################  Whale It! ############################# """
"""               2.       Detect spindles                            """
"""  ################################################################ """
#Run automatic pre-determined spindle algorithms
#.edf in rec_dir and .xml can be in same folder or in another if needed

####### 2.a Set paths
#rec_dir = root_dir + 'DATA/' # folder inside root_dir where BIDS eeg recording data is stored
rec_dir = root_dir + 'DATA/'
xml_dir = root_dir + 'DATA/'# folder inside root_dir where BIDS xml is stored
out_dir = root_dir + 'OUT/spindle/' # folder inside root_dir for output files

####### 2.b set general parameters 
method = [ 'Moelle2011']
#### Options include:
##[Ferrarelli2007, Nir2011, Martin2013, Moelle2011, Wamsley2012, Ray2015, 
## Moelle2011, Lacourse2018, FASST, FASST2, Concordia or UCSD]
chan = ['Fz'] #channel of interest
ref_chan = ['M1', 'M2'] # re-referencing to mastoid : enter [] if no re-referencing needed
grp_name = 'eeg' #if no group name, put None
stage = ['NREM2']
duration=(0.5, 3)
frequency = (10, 16) 
cycle_idx = None #[1,2,3] or None

# Adapted bands
# NOTE: IF USING ADAPTED BANDS, you need to specify a string path to a csv or xlsx sheet 
#                               containing number of columns = 1 + number of channels: 
#                               Column (1) subject id; (2) peak of chan1; (3) peak of chan2 etc.
#                               Each peak column needs to be titled THE SAME as in the variable
#                               chan = [] above.

adap_bands = True # in order of (phase, amplitude)
peaks = root_dir + 'OUT/fooof/fooof_n2_sigmapeaks.csv'
width = 2 # width in Hz either side of the peak that you want to include in bands
if adap_bands:
    frequency = generate_adap_bands(peaks,width,chan)

###### 2.c Concatenation parameters (LEAVE AS DEFAULT)
## Concatenation Setup
# How would you like to concatenate the data where segmentation has been incurred
# e.g. to combine all segments of one sleep stage (e.g. NREM2) across the night 
# e.g. to remove artefacts and stitch back discontinuous data together
concat_cycles = False # If True, will concatenate cycles together (eg cycles 1 + 2 + 3); 
                     # if False, will return ouput per sleep cycle
concat_stages = False # If True, will concatenate different stages together (eg N2 + N3); 
                      # if False, will return ouput per stage
concat_discontinuous = True # If True, will concatenate signal of same type that is 
                            # interrupted in time (eg two bouts of N2 separated by artefact)
concat_evttypes = False # concatenate different event types together 
                        # (Leave False; True is a rare use case)
cat = asarray([concat_cycles, concat_stages, concat_discontinuous, concat_evttypes], 
              dtype=int)

####### 2.d Call wozombie function: whale_it()
whale_it(rec_dir, xml_dir, out_dir, method, chan, ref_chan, rater, cat, stage, grp_name, 
         cycle_idx=cycle_idx, frequency=frequency, adap_bands=adap_bands, duration=duration, 
          part=part, visit=visit)

#%%
"""  #######################  Slow It! ############################## """
"""          3.      Detect slow oscillation                          """
"""  ################################################################ """
# #Run automatic pre-determined slow oscillation algorithms
# #.edf in rec_dir and .xml can be in same folder or in another if needed (ex: could be on previous spindle detection 
#to have all the event together)

####### 3.a Set paths
rec_dir = root_dir  # folder inside root_dir where BIDS eeg recording data is stored
xml_dir = root_dir + 'OUT/spindle/' # folder inside root_dir where BIDS xml is stored (can be same than rec dir or dif)
out_dir = root_dir + 'OUT/spso/' # folder inside root_dir for output files

####### 3.b Set general parameters  
method = ['Staresina2015'] # Method to run - choose one
#### Options include:
##['Ngo2015','Staresina2015','Massimini2004'] - Attention: you can't run all at the same time

chan = ['Cz'] #channel of interest
ref_chan = ['A1', 'A2'] # re-referencing to mastoid : enter [] if no re-referencing needed
grp_name = 'eeg'
stage = ['NREM2', 'NREM3'] #stage of interest
event_name = 'SO' #name of the event in .xml
duration=(0.8, 2)
cycle_idx = None #[1,2,3] or None
average_channels = False #only for Ngo (True)
invert = False # put True to invert polarity

####### 3.c Concatenation parameters
## Concatenation Setup
# How would you like to concatenate the data where segmentation has been incurred
# e.g. to combine all segments of one sleep stage (e.g. NREM2) across the night 
# e.g. to remove artefacts and stitch back discontinuous data together
concat_cycles = True # If True, will concatenate cycles together (eg cycles 1 + 2 + 3); 
                     # if False, will return ouput per sleep cycle
concat_stages = True # If True, will concatenate different stages together (eg N2 + N3); 
                      # if False, will return ouput per stage
concat_discontinuous = True # If True, will concatenate signal of same type that is 
                            # interrupted in time (eg two bouts of N2 separated by artefact)
concat_evttypes = False # concatenate different event types together 
                        # (Leave False; True is a rare use case)
cat = asarray([concat_cycles, concat_stages, concat_discontinuous, concat_evttypes], 
              dtype=int)

####### 3.d Call wozombie function: slow_it()
slow_it(rec_dir, xml_dir, out_dir, method, chan, rater, cat, stage, ref_chan, grp_name, 
            event_name=event_name, cycle_idx=cycle_idx, duration= duration, part=part, 
            visit=visit, average_channels = average_channels, invert = invert)




#%%
"""  #######################  Whale Farm ############################ """
"""           4.    Extract dataset (spindle/SO)                      """
"""  ################################################################ """
# Extract spindle or SO dataset 
#.edf in rec_dir and .xml in folder OUT
# #extract parameters (count, amplitude etc etc) of each event type
# #do it one by one

####### 4.a Set paths
rec_dir = root_dir + 'DATA/' #folder inside root_dir where BIDS eeg recording data is stored
xml_dir = root_dir + 'OUT/so/' #folder inside root_dir where xml with detection are
out_dir = root_dir + 'OUT/' # folder inside root_dir for output files

####### 4.b Set general parameters 
chan = ['C3','Cz', 'C4'] #channel of interest
ref_chan = ['M1', 'M2'] # re-referencing to mastoid : enter [] if no re-referencing needed
grp_name = 'eeg'
cycle_idx = None #  [1,2,3] or enter None to run detection ignoring cycles
stage = ['NREM2', 'NREM3']
segments=None #[('N2_ON','N2_OFF'), ('N3_ON','N3_OFF')] #or None if whole night
evt_name = 'SO' # name of the spindle detection or name of the event (ex: 'SO' or 'Ferrarelli2007')
frequency=(0.1,3.5)
Ngo = {'run':False, 'chan':['av']}
reject_artf=['Artefact', 'Arou', 'Arousal'] #or None if you dont want to reject

####### 3.c Call wozombie function: whale_farm()
whale_farm(rec_dir, xml_dir, out_dir, chan, grp_name, rater, stage=stage, keyword = None,
            ref_chan=ref_chan, evt_name=evt_name, segs=segments, cycle_idx=cycle_idx,frequency=frequency,
            part=part, visit=visit, param_keys=None, exclude_poor=False, 
            reject_artf=reject_artf, epoch_dur=30, n_fft_sec=4, Ngo=Ngo)


#%%
"""  ######################              ############################ """
"""                 5.     EVENT SYNCHRONY                           """
"""  ################################################################ """
#The script takes a target event (e.g. slow oscillations) and calculates the 
#synchrony with a probe event (e.g. spindles) - 
#the stats for the agreement between events (true positives, true negatives, F1 score etc.) 
#are all stored in an output dataset, stats_file

###### 5.a Set Paths
root_dir = r'/Volumes/GoogleDrive/My Drive/Projects/healthy/'
rec_dir = root_dir + 'EDFS/'
xml_dir = root_dir + 'OUT/spso/' 
out_dir = root_dir + 'OUT/complexes/' # for exporting xml with so+sp, so-sp evts

###### 5.b Set general paramters
rater = None # if None, uses first rater
rater_new = None # if None, uses <rater>
evttype_target = 'SO' # target events
evttype_probe = 'Moelle2011_adap' # probe events
iu_thresh = 0.1 # value between 0 and 1, intersection / union
evttype_tp_target = 'SO+mol' # Name for target events that co-occurred with probe
evttype_fn = 'SO-mol' # Name for target events that did not co-occur with probe
grp = 'eeg'

###### 5.c Concatenation Parameters
## Concatenation Setup
# How would you like to concatenate the data where segmentation has been incurred
# e.g. to combine all segments of one sleep stage (e.g. NREM2) across the night 
# e.g. to remove artefacts and stitch back discontinuous data together
concat_cycles = True # If True, will concatenate cycles together (eg cycles 1 + 2 + 3); 
                     # if False, will return ouput per sleep cycle
concat_stages = True # If True, will concatenate different stages together (eg N2 + N3); 
                      # if False, will return ouput per stage
concat_discontinuous = True # (----> For PAC concat_discontinuous should = False)
concat_evttypes = False # concatenate different event types together 
                        # (Leave False; True is a rare use case)
cat = asarray([concat_cycles, concat_stages, concat_discontinuous, concat_evttypes], dtype=int)

####### 5.d Call wozombie function: event_sync()
# Channels must be run 1 at a a time
# Fz
chan = ['Fz']
event_sync(rec_dir, xml_dir, out_dir, part, visit, cat, evttype_target, evttype_probe, 
               iu_thresh, evttype_tp_target, evttype_fn, chan=chan, stage=stage, grp=grp, rater=None)
# Cz
chan = ['Cz']
event_sync(rec_dir, xml_dir, out_dir, part, visit, cat, evttype_target, evttype_probe, 
               iu_thresh, evttype_tp_target, evttype_fn, chan=chan, stage=stage, grp=grp, rater=None)
# Pz
chan = ['_REF']
event_sync(rec_dir, xml_dir, out_dir, part, visit, cat, evttype_target, evttype_probe, 
               iu_thresh, evttype_tp_target, evttype_fn, chan=chan, stage=stage, grp=grp, rater=None)



#%%
"""  ######################              ############################ """
"""          6.       Cross Frequency Coupling   (1 event type)       """
"""  ################################################################ """

####### 6.a Set paths
#root_dir = r'/Volumes/GoogleDrive/My Drive/Projects/healthy/'
rec_dir = root_dir + 'DATA/'
xml_dir = root_dir + 'OUT/spso/'

####### 6.b Set general paramters 
chan = ['Fz','Cz','_REF'] #channel of interest
ref_chan = ['M1', 'M2'] # re-referencing to mastoid : enter [] if no re-referencing needed
grp_name = 'eeg'
rater = None
cycle_idx = None #  [1,2,3] or enter None to run detection ignoring cycles 
stage = ['NREM2','NREM3']
polar = 'normal' # polarity of recording, set to either 'normal' or 'opposite' 
                  # or a list of length (#participants) e.g. ['normal', 'normal', 'opposite']
band_pairs = 'so_sigma'
evt_type = ['star']

###### 6.c Concatenation Setup 
# Concatenation Setup
# How would you like to concatenate the data where segmentation has been incurred
# e.g. to combine all segments of one sleep stage (e.g. NREM2) across the night 
# e.g. to remove artefacts and stitch back discontinuous data together
concat_cycles = True # If True, will concatenate cycles together (eg cycles 1 + 2 + 3); 
                     # if False, will return ouput per sleep cycle
concat_stages = True # If True, will concatenate different stages together (eg N2 + N3); 
                      # if False, will return ouput per stage
concat_discontinuous = False # (----> For PAC concat_discontinuous should = False)
concat_evttypes = False # concatenate different event types together 
                        # (Leave False; True is a rare use case)
cat = asarray([concat_cycles, concat_stages, concat_discontinuous, concat_evttypes], dtype=int)

###### 6.d Filtering  
######### i. Laplacian filtering
laplacian = True # To run surface Laplacian
oREF = 'Pz' # What the online reference was
lapchan = ['F3', 'Fz', 'F4', 'C3', 'C4', 
           'P3', 'P4', 'O1', 'O2', 'F7', 'T5', 'T6', 'T3', 'T4',
           'F8', 'M1','M2','Cz','_REF'] #channels to include in laplacian #channels to include in laplacian
laplacian_rename=False # (Do you need to rename your channels to match the MNE format for laplacian filtering?)
# Format for renaming: {'MNE_name':'chan_name_in_your_edf'}
renames = {'F3':'F3-(M1+M2)', 'Fz':'Fz-(M1+M2)', 'F4':'F4-(M1+M2)', 'C3':'C3-(M1+M2)', 'C4':'C4-(M1+M2)', 
           'P3':'P3-(M1+M2)', 'P4':'P4-(M1+M2)', 'O1':'O1-(M1+M2)', 'O2':'O2-(M1+M2)', 'F7':'F7-(M1+M2)', 
           'F8':'F8-(M1+M2)', 'M1':'M1', 'M2':'M2', 'Pz':'_REF', 'Cz':'Cz'}

######### ii. other filtering
# Apply notch filter?
chan_rename = False  # (Do you need to rename your channels to match the MNE format for notch filtering?)
notch = False
notch_freq = 60 #in Hz
notch_harmonics = False


###### 6.e Advanced parameters 
idpac = (2, 3, 4)       # Choose the combination of methods to use in order to extract PAC 
                        # (see Tensorpac documentation for further info)
peaks = root_dir + 'OUT/fooof/fooof_n2_sigmapeaks.csv' # for adapted peaks
width = 2 # width in Hz either side of the peak that you want to include in bands
nbins = 18              # number of phase bins to calculate
dcomplex = 'hilbert' 
filtcycle = (3, 6)      # for fir1 filter (pha, amp)
min_dur = 0.3           # for minimum duration of events to run PAC across 
buffer = 3              # buffer (length of time) around events to extract to prevent edge artefacts
norm = True

filter_opts={'notch':notch,'notch_harmonics':notch_harmonics, 'notch_freq':notch_freq,
             'laplacian':laplacian, 'lapchan':lapchan,'laplacian_rename':laplacian_rename, 
             'oREF':oREF,'chan_rename':chan_rename,'renames':renames}

###### 6.e Call wozombie function: pac_it() 
#### Calculate mean amplitudes data per subject (saves to pickle .py and .csv files)

# SO - Sigma (fixed)
out_dir = root_dir + 'OUT/cfc_nolap/' # you can set specific output directories, or a general one
adap_bands=(False,False)
fpha = (0.5,1.25)    # Frequency bands for phase frequency
famp = (11,16)      # Frequency bands for amplitude frequency
pac_it(rec_dir, xml_dir, out_dir, part, visit, cycle_idx, chan, rater, stage,
               polar, grp_name, cat, evt_type, buffer, ref_chan, nbins, idpac, 
               fpha, famp, dcomplex, filtcycle, width, min_dur, band_pairs,
               adap_bands=adap_bands, filter_opts=filter_opts,progress=False)

# SO - Sigma (adapted bands)
out_dir = root_dir + 'OUT/cfc_adap_nolap/'
fpha = (0.5,1.25)    # Frequency bands for phase frequency
adap_bands=(False,True) # in order of (phase, amplitude)
famp = generate_adap_bands(peaks,width,chan)
pac_it(rec_dir, xml_dir, out_dir, part, visit, cycle_idx, chan, rater, stage,
               polar, grp_name, cat, evt_type, buffer, ref_chan, nbins, idpac, 
               fpha, famp, dcomplex, filtcycle, width, min_dur, band_pairs,
               adap_bands=adap_bands, filter_opts=filter_opts,progress=False)


#%%
"""  ######################              ############################ """
"""          7.       Cross Frequency Coupling   (2 event types)      """
"""  ################################################################ """
# see wozombie.cfc.mean_amps for more information on the difference 
# between this function and the previous.

####### 7.a Set general parameters 
rec_dir = root_dir + 'EDFS/'
xml_dir = root_dir + 'OUT/complexes/'

# For adapted peaks
peaks = root_dir + 'OUT/fooof/fooof_n2_sigmapeaks.csv'
width = 2 # width in Hz either side of the peak that you want to include in bands

####### 7.b Specify channels and stages 
chan = ['Fz','Cz','_REF'] #channel of interest
oREF = 'Pz'
ref_chan = ['M1', 'M2'] # re-referencing to mastoid : enter [] if no re-referencing needed
grp_name = 'eeg'
rater = None
cycle_idx = None #  [1,2,3] or enter None to run detection ignoring cycles 
stage = ['NREM2','NREM3']
polar = 'normal' # polarity of recording, set to either 'normal' or 'opposite' 
                  # or a list of length (#participants) e.g. ['normal', 'normal', 'opposite']

####### 7.c Concatenation paramters
concat_cycles = True # If True, will concatenate cycles together (eg cycles 1 + 2 + 3); 
                     # if False, will return ouput per sleep cycle
concat_stages = True # If True, will concatenate different stages together (eg N2 + N3); 
                      # if False, will return ouput per stage
concat_discontinuous = False # (----> For PAC concat_discontinuous should = False)
concat_evttypes = False # concatenate different event types together 
                        # (Leave False; True is a rare use case)
cat = asarray([concat_cycles, concat_stages, concat_discontinuous, concat_evttypes], dtype=int)

####### 7.d Advanced parameters 
nbins = 18
idpac = (2,3,4)
dcomplex = 'hilbert' 
filtcycle = (3, 6)      # for fir1 filter (pha, amp)
min_dur = 0.3           # for minimum duration of events to run PAC across 
buffer = 3              # buffer (length of time) around events to extract to prevent edge artefacts
norm = True

filter_opts={'notch':False,'notch_harmonics':False, 'notch_freq':60,
             'laplacian':False, 'lapchan':None,'laplacian_rename':False, 
             'oREF':oREF,'chan_rename':chan_rename,'renames':None}


####### 7.e Call wozombie function: pac_it_2() 
#### Calculate mean amplitudes data per subject (saves to pickle .py and .csv files)

# SO - Moelle2011 (adapted)
out_dir = root_dir + 'OUT/evt_cfc_adap/' # you can set specific output directories, or a general one
target = 'SO+mol'
probe = 'Moelle2011_adap'
band_pairs = 'so_sigma_mol'
adap_bands=(False,True) # in order of (phase, amplitude)
fpha = (0.1,1.25)
famp = generate_adap_bands(peaks,width,chan)
pac_it_2(rec_dir, xml_dir, out_dir, part, visit, cycle_idx, chan, rater, stage,
               polar, grp_name, cat, target, probe, buffer, ref_chan, nbins, idpac, 
               fpha, famp, dcomplex, filtcycle, width, min_dur, band_pairs,
               adap_bands=adap_bands,filter_opts=filter_opts,progress=False)

# SO - Moelle2011 (fixed)
out_dir = root_dir + 'OUT/evt_cfc/'
target = 'SO+ray'
probe = 'Ray2015_adap'
band_pairs = 'so_sigma_ray'
adap_bands=(False,False) # in order of (phase, amplitude)
fpha = (0.1,1.25)
famp = (11,16)
pac_it_2(rec_dir, xml_dir, out_dir, part, visit, cycle_idx, chan, rater, stage,
               polar, grp_name, cat, target, probe, buffer, ref_chan, nbins, idpac, 
               fpha, famp, dcomplex, filtcycle, width, min_dur, band_pairs,
               adap_bands=adap_bands,filter_opts=filter_opts,progress=False)


#%% 
"""  ######################              ############################ """
"""          8.       Cross Frequency Coupling   (plotting)           """
"""  ################################################################ """

### Plots mean amplitudes and preferred phase for all participants and visits.
### Figures are automatically saved to {out_dir}

## Plotting parameters 
colors = ['darkblue','darkred'] # colors for visits
band_pairs = ['so_sigma']

plot_mean_amps(in_dir, out_dir, part, visit, chan, stage, band_pairs, colors, 
               nbins)
plot_prefphase(in_dir, out_dir, part, visit, chan, stage, band_pairs, colors, 
               nbins)


### Plots a figure of the group-level preferred phase of different comparisons.
'''
How to set:
comps = the comparisons you would like to plot, in the format
        of [(participants, visit)] 
        
        e.g. [('all', ['visit1']), 
              ('all', ['visit2'])]
        
        or   [(['HC001','HC002'], ['visit1']),
              (['PT001','PT002'], ['visit1'])]
        
'''  
comps = [('all',['V1']), ('all',['V2'])]
colors = ['darkblue','darkred'] # colors for comps (must be same no. as no. of comps)

plot_prefphase_group(in_dir, out_dir, band_pairs, chan, cycle_idx, stage, nbins, 
                     norm, colors = colors, comps = comps, layout=None)
plot_meanamps_group(in_dir, out_dir, band_pairs, chan, cycle_idx, stage, nbins, 
                    norm, colors = colors, comps = comps, layout=None)   


                        
