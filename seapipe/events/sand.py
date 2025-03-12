#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 18:13:45 2024

@author: ncro8394
"""

from numpy import append, diff, insert, percentile, where 
from os import listdir, mkdir, path
from wonambi import Dataset, graphoelement
from wonambi.attr import Annotations, create_empty_annotations
from wonambi.detect.spindle import transform_signal
from wonambi.trans import fetch
import mne
from scipy.signal import find_peaks, peak_widths
import yasa
from numpy import (array, concatenate, convolve, cumsum, full, multiply, nan, 
                   ones, repeat, roll, where, zeros)
from copy import deepcopy
from ..utils.logs import create_logger
from ..utils.load import load_channels, load_sessions, rename_channels
from ..utils.misc import merge_events, reconstruct_stitches, remove_duplicate_evts


def detect_above_zero_regions(signal):
    # Find transitions
    starts = where(diff(signal.astype(int)) >0)[0] + 1
    ends = where(diff(signal.astype(int)) <0)[0] + 1

    # Edge case: If signal starts above 0
    if signal[0]:
        starts = insert(starts, 0, 0)
    
    # Edge case: If signal ends above 0
    if signal[-1]:
        ends = append(ends, len(signal))

    return list(zip(starts, ends))    
    


class SAND:
    
    """ Seapipe Artefact and Noise Detection (S.A.N.D)

        This module runs automated artefact detection with the option of using
        previously published staging algorithms:
            1. YASA (standard deviation)
            2. YASA (covariance)
            3. (More to come..)
        
    """   
    
    def __init__(self, rec_dir, xml_dir, out_dir, eeg_chan, ref_chan,
                 rater = None, grp_name = 'eeg', 
                 subs='all', sessions='all', tracking = None):
        
        self.rec_dir = rec_dir
        self.xml_dir = xml_dir
        self.out_dir = out_dir
        self.eeg_chan = eeg_chan
        self.ref_chan = ref_chan
        self.rater = rater
        self.grp_name = grp_name
        
        self.subs = subs
        self.sessions = sessions
        
        if tracking == None:
            tracking = {}
        self.tracking = tracking


    def detect_artefacts(self, method, label = "allchans", win_size = 5,
                               filetype = '.edf', 
                               stage = ['NREM1', 'NREM2', 'NREM3', 'REM'],
                               logger = create_logger('Detect artefacts')):
        
        ''' Automatically detects artefacts.
        
            Creates a new annotations file if one doesn't already exist.
        
        INPUTS:
            
            method      ->   str of name of automated detection algorithm to 
                             detect staging with. 
                             Current methods supported: 
                                 1. 'Vallat2021' (https://doi.org/10.7554/eLife.70092)
                                 2. 'Cross2025' 
                             
            qual_thresh ->   Quality threshold. Any stages with a confidence of 
                             prediction lower than this threshold will be set 
                             to 'Undefined' for futher manual review.
   
        
        '''
        
        ### 0.a Set up logging
        flag = 0
        tracking = self.tracking
        if label == "allchans":
            chan_msg = "All channels at once."
        else:
            chan_msg = "Each channel individually."
            
        logger.info('')
        logger.debug(rf"""Commencing artefact detection... 
                     
                                             ____
                                      /^\   / -- )
                                     / | \ (____/
                                    / | | \ / /
                                   /_|_|_|_/ /
                                    |     / /
                     __    __    __ |    / /__    __    __
                    [  ]__[  ]__[  ].   / /[  ]__[  ]__[  ]     ......
                    |__            ____/ /___           __|    .......
                       |          / .------  )         |     ..........
                       |         / /        /          |    ............
                       |        / /        / _         |  ...............
                   ~._..-~._,….-ˆ‘ˆ˝\_,~._;––' \_.~.~._.~'\................  
                       
            
                    Seapipe Artefact and Noise Detection
                    (S.A.N.D)

                    Method: {method}
                    
                    Applying to: {chan_msg} 
                    
                                                    """,)
        ### 1. First we check the directories
        # a. Check for output folder, if doesn't exist, create
        if path.exists(self.out_dir):
                logger.debug("Output directory: " + self.out_dir + " exists")
        else:
            mkdir(self.out_dir)
        
        # b. Check input list
        subs = self.subs
        if isinstance(subs, list):
            None
        elif subs == 'all':
                subs = listdir(self.rec_dir)
                subs = [p for p in subs if not '.' in p]
        else:
            logger.error("'subs' must either be an array of subject ids or = 'all' ")  
            return
        
        ### 2. Begin loop through dataset
       
        # a. Begin loop through participants
        subs.sort()
        for i, sub in enumerate(subs):
            tracking[f'{sub}'] = {}
            # b. Begin loop through sessions
            flag, sessions = load_sessions(sub, self.sessions, self.rec_dir, flag, 
                                     logger, verbose=2)
            for v, ses in enumerate(sessions):
                logger.info('')
                logger.debug(f'Commencing {sub}, {ses}')
                tracking[f'{sub}'][f'{ses}'] = {'slowosc':{}} 
                
                # Define recording
                rdir = f'{self.rec_dir}/{sub}/{ses}/eeg/'
                try:
                    edf_file = [x for x in listdir(rdir) if x.endswith(filetype)][0]
                except:
                    logger.warning(f'No input {filetype} file in {rdir}')
                    flag += 1
                    break
                
                ## f. Channel setup 
                pflag = deepcopy(flag)
                flag, chanset = load_channels(sub, ses, self.eeg_chan, 
                                              self.ref_chan, flag, logger)
                if flag - pflag > 0:
                    logger.warning(f'Skipping {sub}, {ses}...')
                    break
                
                newchans = rename_channels(sub, ses, self.eeg_chan, logger) 
                
                # Check if applying to all channels or chan-by-chan
                if label == "allchans":
                    # Check if references are the same for each channel
                    ref_chans = {tuple(val) for val in chanset.values()}
                    # If not, setup for running per channel
                    if len(ref_chans) > 1:
                        logger.warning("Channel setup 'all_chans' was set 'True', but "
                                       f"Channel:Reference pairings are unique for {sub}, "
                                       f"{ses}. Therefore, we'll detect artefacts PER CHANNEL.")
                        flag += 1
                        label = "individual"

                # d. Load/create for annotations file
                if not path.exists(self.xml_dir + '/' + sub):
                    mkdir(self.xml_dir + '/' + sub)
                if not path.exists(self.xml_dir + '/' + sub + '/' + ses):
                     mkdir(self.xml_dir + '/' + sub + '/' + ses)
                xdir = self.xml_dir + '/' + sub + '/' + ses
                
                xml_file = [x for x in listdir(f'{xdir}') if '.xml' in x]
                if len(xml_file) > 1:
                    logger.warning(f'More than 1 annotations file found for '
                                   f'{sub}, {ses} in {xdir}. Skipping...')
                    continue
                if len(xml_file) < 1:
                    logger.warning(f'No annotations file was found for '
                                   f'{sub}, {ses} in {xdir}. Skipping...')
                else:
                    xml_file = f'{xdir}/{xml_file[0]}'
                    
                if not path.exists(xml_file):
                    dset = Dataset(rdir + edf_file)
                    create_empty_annotations(xml_file, dset)
                    logger.warning(f"No annotations file exists. Creating " 
                                   f"annotations file for {sub}, {ses} and" 
                                   "detecting Artefacts WITHOUT hypnogram.")
                    annot = Annotations(xml_file)
                    hypno = None
                else:
                    logger.debug(f'Annotations file exists for {sub}, {ses},'
                                 'staging will be used for Artefact detection.')

                    # Extract hypnogram
                    annot = Annotations(xml_file)
                    hypno = [x['stage'] for x in annot.get_epochs()]
                    stage_key = {'Wake':0,
                                 'NREM1':1,
                                 'NREM2':2,
                                 'NREM3':3,
                                 'REM':4,
                                 'Undefined':0,
                                 'Unknown':0,
                                 'Artefact':0}
                    
                    hypno = array([int(stage_key[x]) for x in hypno])
                
                
                for chan in chanset:
                    if newchans:
                        logger.debug(f'Detecting for {newchans[chan]} : {chanset[chan]}')
                    else:
                        logger.debug(f'Detecting for {chan} : {chanset[chan]}')

                    if 'yasa' in method:     
                        
                        ## c. Load recording
                        try:
                            raw = mne.io.read_raw_edf(rdir + edf_file, 
                                                      include = chan + chanset[chan],
                                                      preload=True, verbose = False)
                            
                            mne.set_eeg_reference(raw, ref_channels = chanset[chan])
                                                  
                            s_freq = raw.info["sfreq"]
                        except Exception as e:
                            logger.warning(f'Error loading {filetype} file in {rdir}, {repr(e)}')
                            flag += 1
                            break
                        
                        yasa_meth = 'covar' if 'covar' in method else 'std'
                            
                        # Convert raw data to array    
                        data = raw.to_data_frame()
                        inds = [x for x in data if x in chan]
                        data = data[inds].T
                        data = data.to_numpy()
                        
                        # Upsample hypnogram to match raw data
                        hypno_up = yasa.hypno_upsample_to_data(hypno, 1/30, data, 
                                                               sf_data=s_freq)
                        
                        # Detect artefacts
                        n_chan_reject = 1 if data.shape[0] == 1 else 2
                        art, zscores = yasa.art_detect(data, s_freq, 
                                                       window = win_size, 
                                                       hypno = hypno_up, 
                                                       include = (1, 2, 3, 4), 
                                                       method = yasa_meth, 
                                                       threshold = 3, 
                                                       n_chan_reject = n_chan_reject, 
                                                       verbose = False)
                        
                        # Upsample artefacts to match raw data
                        art = multiply(art, 1)
                        sf_art = 1/win_size
                        art_up = yasa.hypno_upsample_to_data(art, sf_art, data, 
                                                             s_freq)
                        
                        # Find start/end times of artefacts
                        peaks = find_peaks(art_up)
                        properties = peak_widths(art_up, peaks[0])
                        times = [x for x in zip(properties[2],properties[3])]
        
                    elif 'seapipe' in method:
                        
                        logger.debug('Loading Data..')
                        
                        # Get data 
                        dset = Dataset(rdir + edf_file)
                        s_freq = int(dset.header['s_freq'])
                        segments = fetch(dset, annot, cat = (1,1,1,1), 
                                         stage=stage)
                        
                        # Save stitches (to recompute times later)
                        stitches = segments[0]['times']
                        
                        # Read data
                        segments.read_data(chan, ref_chan = chanset[chan], 
                                           grp_name=self.grp_name)
                        data = segments[0]['data'].data[0]
                        
                        # Filter data above 40Hz
                        dat = transform_signal(data[0], s_freq, 'high_butter', 
                                         method_opt={'freq':40,
                                                     'order':3})
                        

                        ## Step 1. Detect flatlines in data
                        logger.debug("Detecting flatlines...")
                        # Align the filtered data back to the raw EEG recording
                        dat_full = reconstruct_stitches(dat, stitches, s_freq,
                                                        replacement=nan)
                        
                        # Let's setup a function to detect flatlines
                        def detect_flatlines(ts, s_freq, tolerance=1e-5, min_length=5):
                            diffs = abs(diff(ts))
                            is_flat = diffs < tolerance
                        
                            # Find indices where flatlines start and stop
                            change_points = where(diff(concatenate(([False], 
                                                                     is_flat, 
                                                                     [False]))))[0]
                            starts, ends = change_points[::2], change_points[1::2]
                        
                            # Filter out short flatlines and adjust for s_freq offset
                            flatlines = [(max(0, start - int(s_freq / 4)), min(len(ts), 
                                                                               end + int(s_freq / 4))) 
                                         for start, end in zip(starts, ends) if 
                                         (end - start) >= min_length]
                        
                            return flatlines 
                        flatlines = detect_flatlines(dat_full, s_freq)
                        
                        
                        ## Step 2. Find large shifts in fast frequencies (>40Hz)
                        logger.debug("Detecting movement artefacts...")
                        # Mask out negative values
                        dat[dat<0] = 0
                        
                        # Calculate sliding RMS
                        dat = transform_signal(dat, s_freq, 'moving_rms',
                                         method_opt = {'dur':win_size,
                                                       'step':1})
                        
                        # Put back all the samples that were averaged within the window 
                        dat = repeat(dat, s_freq)

                        # Filter (smooth) the RMS
                        dat = transform_signal(dat, s_freq, 'low_butter', 
                                         method_opt={'freq':0.5,
                                                     'order':3})

                        # Align this filtered data back to the raw EEG recording 
                        dat_reconstructed = reconstruct_stitches(dat, 
                                                                  stitches, 
                                                                  s_freq)
                        
                        # Threshold and detect artefact 'events'
                        threshold = percentile(dat, 95)
                        dat_reconstructed[dat_reconstructed<threshold] = 0
                        dat_reconstructed[dat_reconstructed>threshold] = 1
                        movements = detect_above_zero_regions(dat_reconstructed)
                        
                        
                        ## TODO: Step 3. Let's look for sharp slow waves in REM
                        # if 'REM' in stage:
                        #     logger.debug("Detecting eye movement artefacts in REM...")
                        #     segments = fetch(dset, annot, cat = (1,1,1,1), 
                        #                      stage=['REM'])
                            
                        #     # Save stitches (to recompute times later)
                        #     stitches = segments[0]['times']
                        #     # Read data
                        #     segments.read_data(chan, ref_chan=ref, 
                        #                        grp_name=self.grp_name)
                        #     data = segments[0]['data'].data[0]
                            
                        #     # Filter (smooth) the RMS
                        #     dat = transform_signal(data[3], s_freq, 'low_butter', 
                        #                      method_opt={'freq':5,
                        #                                  'order':3})
                            
                        #     dat = transform_signal(dat, s_freq, 'moving_sd',
                        #                      method_opt = {'dur':3,
                        #                                    'step':1})
                            
                        #     # Put back all the samples that were averaged within the window 
                        #     dat = repeat(dat, s_freq)
                            
                            
                        #     # Align this filtered data back to the raw EEG recording
                        #     dat_reconstructed = reconstruct_stitches(dat, 
                        #                                               stitches, 
                        #                                               s_freq)
                            
                            
                        #     # Threshold and detect artefact 'events'
                        #     threshold = 50
                        #     dat_reconstructed[dat_reconstructed<50] = 0
                        #     dat_reconstructed[dat_reconstructed>50] = 1
                        #     REMS = detect_above_zero_regions(dat_reconstructed)
                            

                        ## Step 4. Combine all artefact types
                        times = flatlines + movements
                        
                        times = merge_segments(times)
                        
                        # Convert times back into seconds for annotations
                        times = [(x / s_freq, y / s_freq) for x, y in times]

                        
                        # Convert to wonambi annotations format

                    else:
                        logger.critical("Currently the only methods that are" 
                                        " functioning include:"
                                        "\n'yasa_std', 'yasa_covar', or 'seapipe'."
                                        ) 
                        return
                    
                    # Save events to Annotations file
                    if label == "allchans":
                        chan_name = ['']
                    else:
                        chan_name = [chan]
                    evts = []
                    for x in times:
                        evts.append({'name':'Artefact',
                              'start':float(x[0]),
                              'end':float(x[1]),
                              'chan':chan_name,
                              'stage':'',
                              'quality':'Good',
                              'cycle':''})
                        
                    # Add to annotations file
                    grapho = graphoelement.Graphoelement()
                    grapho.events = evts          
                    grapho.to_annot(annot)
                    
                merge_events(annot, 'Artefact')
                remove_duplicate_evts(annot, 'Artefact')

        return
    
    
    
def merge_segments(segments):
    """
    Merges overlapping or adjacent segments in a list of (start, end) tuples.

    Parameters:
      segments (list of tuples): List of (start, end) index pairs.

    Returns:
      list of tuples: Merged (start, end) segments.
    """
    if not segments:
        return []

    # Sort segments by start time
    segments.sort()

    # Initialize merged list with the first segment
    merged = [segments[0]]

    for start, end in segments[1:]:
        prev_start, prev_end = merged[-1]

        if start <= prev_end:  # Overlapping or adjacent segments
            merged[-1] = (prev_start, max(prev_end, end))  # Merge them
        else:
            merged.append((start, end))  # No overlap, add as new segment

    return merged
    
    
    