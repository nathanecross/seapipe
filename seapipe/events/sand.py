#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 18:13:45 2024

@author: ncro8394
"""

from numpy import (append, array, asarray, concatenate, cumsum, diff, inf, 
                   insert, isnan, median, multiply, nan, nanmean, nanstd, pad, 
                   percentile, repeat, sqrt, squeeze, where)
from os import listdir, mkdir, path
from wonambi import Dataset, graphoelement
from wonambi.attr import Annotations, create_empty_annotations
from wonambi.detect.spindle import transform_signal
from wonambi.trans import fetch
import mne
from scipy.signal import butter, correlate, filtfilt, find_peaks, peak_widths
import yasa
from copy import deepcopy
from ..utils.logs import create_logger
from ..utils.load import load_channels, load_sessions, rename_channels
from ..utils.misc import (merge_events, reconstruct_stitches, remove_duplicate_evts)


class SAND:
    
    """ Seapipe Artefact and Noise Detection (S.A.N.D)

        This module runs automated artefact detection with the option of using
        previously published staging algorithms:
            1. YASA (standard deviation)
            2. YASA (covariance)
            3. Seapipe artefact detection
        
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
                                 2. 'seapipe' 
                             
            qual_thresh ->   Quality threshold. Any stages with a confidence of 
                             prediction lower than this threshold will be set 
                             to 'Undefined' for futher manual review.
   
        
        '''
        
        ### 0.a Set up logging
        flag = 0
        tracking = self.tracking
        if label == "allchans":
            chan_msg = "all channels at once."
        else:
            chan_msg = "each channel individually."
            
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
                    continue
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
                
                evts = []
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
                        dat_full = reconstruct_stitches(dat, stitches, s_freq)
                        
                        # Threshold and detect artefact 'events'
                        threshold = percentile(dat, 95)
                        dat_full[dat_full<threshold] = 0
                        dat_full[dat_full>threshold] = 1
                        movements = detect_above_zero_regions(dat_full)
                        
                        
                        # Step 3. Large deflections in signal
                        logger.debug("Detecting deflection artefacts...")
                        
                        dat = append(0, squeeze(diff(data)))
                        dat_full = reconstruct_stitches(dat, stitches, s_freq,
                                                         replacement=nan)
                        
                        threshold = 10 * nanstd(dat_full)
                        dat_full[dat_full<threshold] = 0
                        dat_full[dat_full>=threshold] = 1
                        deflections = detect_above_zero_regions(dat_full)
                        deflections = [(x[0] - (s_freq/2), 
                                        x[1] + (s_freq/2)) for x in deflections]
                        
                        # dat = transform_signal(data[0], s_freq, 'low_butter', 
                        #                  method_opt={'freq':2,
                        #                              'order':3})
                        # dat = transform_signal(dat, s_freq, 'moving_rms', 
                        #                  method_opt = {'dur':win_size,
                        #                                'step':1})
                        # # Put back all the samples that were averaged within the window 
                        # dat = repeat(dat, s_freq)
                        
                        # # Align the filtered data back to the raw EEG recording
                        # dat_full = reconstruct_stitches(dat, stitches, s_freq,
                        #                                 replacement=nan)
                        
                        # #threshold = percentile(dat, 95)
                        # threshold = 2*nanstd(dat_full)
                        # dat_full[dat_full<threshold] = 0
                        # dat_full[dat_full>=threshold] = 1
                        # deflections = detect_above_zero_regions(dat_full)
                        
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
                        
                        merged = flatlines + deflections +  movements
                        
                        merged = merge_segments(merged)
                        
                        # Convert times back into seconds for annotations
                        times = [(x / s_freq, y / s_freq) for x, y in merged]

                        # Convert to wonambi annotations format

                    else:
                        logger.critical("Currently the only methods that are" 
                                        " functioning include:"
                                        "\n'yasa_std', 'yasa_covar', or 'seapipe'."
                                        ) 
                        return
                    
                    # Save events to Annotations file
                    if label == "allchans":
                        channame = ['']
                    else:
                        channame = f'{chan} ({self.grp_name})'
                        
                    for x in times:
                        evts.append({'name':'Artefact',
                                     'start':float(x[0]),
                                     'end':float(x[1]),
                                     'chan':channame,
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
    
    
    
def detect_ecg_artifacts(eeg_signal, sfreq, bandpass=(1, 40), template_window=0.6,
                         peak_percentile=95, match_percentile=99, min_rr_interval=0.5):
    """
    Automatically detects ECG artifacts in a single EEG channel using template matching.

    Parameters:
    - eeg_signal: 1D numpy array of EEG data
    - sfreq: Sampling frequency in Hz
    - bandpass: Tuple of (low, high) bandpass filter cutoff in Hz
    - template_window: Length of ECG template window in seconds
    - peak_percentile: Percentile threshold for peak detection
    - match_percentile: Correlation percentile threshold for match detection
    - min_rr_interval: Minimum interval between matches in seconds

    Returns:
    - clean_matches: Sample indices of detected ECG artifacts
    - template: The derived ECG template used for matching
    - corr: Cross-correlation trace
    """
    
    
    # 1. Bandpass filter
    def bandpass_filter(data, low, high, fs, order=4):
        nyq = 0.5 * fs
        b, a = butter(order, [low / nyq, high / nyq], btype='band')
        return filtfilt(b, a, data)

    filtered = bandpass_filter(eeg_signal, *bandpass, fs=sfreq)

    # 2. Detect peaks likely to be ECG R-peaks
    min_distance = int(min_rr_interval * sfreq)
    height_thresh = percentile(filtered, peak_percentile)
    peaks, _ = find_peaks(filtered, 
                          height = height_thresh, 
                          distance = min_distance)

    # 3. Extract candidate segments around each peak
    win_len = int(template_window * sfreq)
    half_win = win_len // 2
    segments = []

    for peak in peaks:
        if peak - half_win >= 0 and peak + half_win < len(filtered):
            segment = filtered[peak - half_win:peak + half_win]
            segments.append(segment)

    if not segments:
        raise ValueError("No valid ECG segments found. Try lowering peak_percentile or check EEG signal.")

    segments = array(segments)

    # 4. Create the ECG template
    template = median(segments, axis=0)

    # 5. Cross-correlate the template with the EEG signal
    corr = fast_normalized_cross_correlation(filtered, template)
    threshold = percentile(corr[~isnan(corr)], match_percentile)
    matches = where(corr > threshold)[0]

    # 6. Prune close matches
    pruned = []
    last = -inf
    for m in matches:
        if m - last > min_distance:
            pruned.append(m)
            last = m

    clean_matches = array(pruned)

    return clean_matches, template, corr    
    

def fast_normalized_cross_correlation(signal, template):
    """
    Efficient normalized cross-correlation using cumulative sums.

    Parameters:
    - signal: 1D array
    - template: 1D array (shorter than signal)

    Returns:
    - norm_corr: Normalized cross-correlation (same length as signal)
    """
    signal = asarray(signal)
    template = asarray(template)
    template_len = len(template)
    template = (template - nanmean(template)) / nanstd(template)

    # Raw cross-correlation (same length as signal)
    raw_corr = correlate(signal, template, mode='same')

    # Compute running mean and std of signal (same window length as template)
    window = template_len
    runmean = cumsum(insert(signal, 0, 0))
    runmean2 = cumsum(insert(signal**2, 0, 0))

    sum_window = runmean[window:] - runmean[:-window]
    sum2_window = runmean2[window:] - runmean2[:-window]

    mean_window = sum_window / window
    std_window = sqrt(sum2_window / (window - mean_window**2))

    # Pad to match length of signal
    std_padded = pad(std_window, (window // 2, len(signal) - len(std_window) - 
                                  window // 2), mode='edge')

    # Avoid division by zero
    std_padded[std_padded == 0] = 1e-10

    norm_corr = raw_corr / (std_padded * window)

    return norm_corr


def detect_above_zero_regions(signal):
    """
    Detects regions where the signal is above zero.

    Parameters:
    - signal: 1D array
        Input signal where the regions of interest (above zero) are detected.

    Returns:
    - regions: list of tuples
        Each tuple contains the start and end indices of regions where the signal is above zero.
    """
    signal = asarray(signal)
    valid = ~isnan(signal)

    # Compute diffs and mask invalid transitions
    diffs = diff(signal.astype(float))
    valid_diffs = valid[:-1] & valid[1:]

    raw_starts = where((diffs > 0) & valid_diffs)[0] + 1
    raw_ends = where((diffs < 0) & valid_diffs)[0] + 1

    # Edge case: signal starts above 0 and is valid
    if valid[0] and signal[0] > 0:
        raw_starts = insert(raw_starts, 0, 0)

    # Edge case: signal ends above 0 and is valid
    if valid[-1] and signal[-1] > 0:
        raw_ends = append(raw_ends, len(signal))

    # Now, match starts to ends
    starts = []
    ends = []
    s_idx = 0
    e_idx = 0

    while s_idx < len(raw_starts) and e_idx < len(raw_ends):
        if raw_starts[s_idx] < raw_ends[e_idx]:
            starts.append(raw_starts[s_idx])
            # find first end after this start
            while e_idx < len(raw_ends) and raw_ends[e_idx] < raw_starts[s_idx]:
                e_idx += 1
            if e_idx < len(raw_ends):
                ends.append(raw_ends[e_idx])
                e_idx += 1
            s_idx += 1
        else:
            # skip any unmatched end that comes before the first start
            e_idx += 1

    return list(zip(starts, ends))

    

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
