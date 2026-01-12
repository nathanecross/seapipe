# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 10:34:53 2021

@author: nathancross
"""
from copy import deepcopy
from datetime import datetime, date
from itertools import product
from os import listdir, mkdir, path
import shutil
from pathlib import Path
import csv

import numpy as np
from wonambi import Dataset
from wonambi.attr import Annotations
from wonambi.detect import consensus, DetectSpindle
from wonambi.trans import fetch
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
from ..utils.logs import create_logger, create_logger_outfile
from ..utils.load import (load_channels, load_adap_bands, load_sessions, 
                          rename_channels, read_manual_peaks)
from ..utils.misc import remove_duplicate_evts


class whales:
    
    """ Wonambi Heuristic Approach to Locating Elementary Spindles (WHALES)

        This module runs a consensus approach to detecting sleep spindles. While we hope
        to improve detection, and remove biases that occur based on the use of any one 
        spindle detector, this is not a perfect solution.
        The pipeline runs in three stages:
            1. whale_it: Detect spindles with multiple published algorithms 
                (see Documentation).
            2. whales: Assign 'true' spindle events based upon a pre-set 
                agreement threshold, using a consensus of the events detected 
                independently from step 1. This creates a new event called
                'spindle' in the annotations file.
               
    """   
    
    def __init__(self, rootpath, rec_dir, xml_dir, out_dir, chan, ref_chan, 
                 grp_name, stage, frequency, rater = None, subs = 'all', 
                 sessions = 'all', reject_artf = ['Artefact', 'Arou', 'Arousal'], 
                 tracking = None):
        
        self.rootpath = rootpath
        self.rec_dir = rec_dir
        self.xml_dir = xml_dir
        self.out_dir = out_dir
        self.chan = chan
        self.ref_chan = ref_chan
        self.grp_name = grp_name
        self.stage = stage
        self.frequency = frequency
        self.rater = rater
        self.reject = reject_artf
        
        self.subs = subs
        self.sessions = sessions
        
        if tracking == None:
            tracking = {'spindle':{}}
        self.tracking = tracking

    def whale_it(self, method, cat, cycle_idx = None, adap_bands = False, 
                       adap_bw = 4, duration = (0.5, 3), 
                       filetype = '.edf', logger = create_logger('Detect spindles')):
        
        '''
        Runs one (or multiple) automatic spindle detection algorithms and saves the 
        detected events to a new annotations file.
        
        INPUTS:
            method ->    List of names of automated detection algorithms to detect 
                         events with. e.g. ['Lacourse2018','Moelle2011']
            cat ->       Tuple of 4 digits of either 0 or 1 e.g. (0,1,1,0) 
                         This variable sets the concatenation type when reading in 
                         the data. 0 means no concatenation, 1 means concatenation
                         #position 1: cycle concatenation
                         #position 2: stage concatenation
                         #position 3: discontinuous signal concatenation
                         #position 4: event type concatenation (does not apply here)
                         Set this based on whether you would like to detect across the
                         entire recording e.g. (1,1,1,1) or separately for each cycle
                         e.g. (0,1,1,1) or separately for each stage e.g. (1,0,1,1)
            cycle_idx->  List of indices corresponding to sleep cycle numbers.
            duration ->  Tuple of 2 digits, for the minimum and maximum duration of any 
                         detected events.     
        
        '''
        
        ### 0.a Set up logging
        flag = 0
        logger.info('')
        logger.debug(r"""Commencing spindle detection... 
                     
                                
                                          .
                                      .  • •  .. 
                                  .  • • • • • •  .   
                              •. • • • • • • • • • •        
                           .•  • • • • • • • • • • •   .•..           .•.
           . .      .•.   •    • • • • • • • • • • .  .    •.  .    .•   •. .       
              •. .•    •.•     • • • • • • • • • • •.•       .• •..•         
                               • . • • • • • . •.•             
                                •  • . • .  •            
                                    •   •

                                
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
        
        ### 2. Begin loop through dataset
       
        # a. Begin loop through participants
        subs.sort()
        for i, sub in enumerate(subs):
            if not sub in self.tracking['spindle'].keys():
                self.tracking['spindle'][sub] = {}
            
            # b. Begin loop through sessions
            flag, sessions = load_sessions(sub, self.sessions, self.rec_dir, flag, 
                                     logger, verbose=2)   
            for v, ses in enumerate(sessions):
                logger.info('')
                logger.debug(f'Commencing {sub}, {ses}')
                if not ses in self.tracking['spindle'][sub].keys():
                    self.tracking['spindle'][sub][ses] = {} 
    
                ## c. Load recording
                rdir = self.rec_dir + '/' + sub + '/' + ses + '/eeg/'
                try:
                    edf_file = [x for x in listdir(rdir) if x.endswith(filetype)]
                    dset = Dataset(rdir + edf_file[0])
                except:
                    logger.warning(f' No input {filetype} file in {rdir}')
                    flag+=1
                    break
                
                ## d. Load annotations
                xdir = self.xml_dir + '/' + sub + '/' + ses + '/'
                xml_file = [x for x in listdir(xdir) if x.endswith('.xml')][0]
                # Copy annotations file before beginning
                if not path.exists(self.out_dir):
                    mkdir(self.out_dir)
                if not path.exists(self.out_dir + '/' + sub):
                    mkdir(self.out_dir + '/' + sub)
                if not path.exists(self.out_dir + '/' + sub + '/' + ses):
                    mkdir(self.out_dir + '/' + sub + '/' + ses)
                backup = self.out_dir + '/' + sub + '/' + ses + '/'
                backup_file = (f'{backup}{sub}_{ses}.xml')
                if not path.exists(backup_file):
                    shutil.copy(xdir + xml_file, backup_file)
                else:
                    logger.warning(f'Annotations file already exists for {sub}, {ses}, any previously detected events will be overwritten.')
                    flag += 1
                # Read annotations file
                annot = Annotations(backup_file, rater_name=self.rater)
                
                ## e. Get sleep cycles (if any)
                if cycle_idx is not None:
                    all_cycles = annot.get_cycles()
                    cycle = [all_cycles[y - 1] for y in cycle_idx if y <= len(all_cycles)]
                else:
                    cycle = None
                
                ## f. Channel setup 
                pflag = deepcopy(flag)
                flag, chanset = load_channels(sub, ses, self.chan, self.ref_chan,
                                              flag, logger)
                if flag - pflag > 0:
                    logger.warning(f'Skipping {sub}, {ses}...')
                    flag += 1
                    continue
                
                newchans = rename_channels(sub, ses, self.chan, logger)
                
                for c, ch in enumerate(chanset):
                    
                    # 5.b Rename channel for output file (if required)
                    if newchans:
                        fnamechan = newchans[ch]
                    else:
                        fnamechan = ch
                        
                    # g. Check for adapted bands
                    if adap_bands == 'Fixed':
                        freq = self.frequency
                    elif adap_bands == 'Manual':
                        freq = read_manual_peaks(self.rootpath, sub, ses, ch, 
                                                 adap_bw, logger)
                    elif adap_bands == 'Auto':
                        stagename = '-'.join(self.stage)
                        band_limits = f'{self.frequency[0]}-{self.frequency[1]}Hz'
                        freq = load_adap_bands(self.tracking['fooof'], sub, ses,
                                               fnamechan, stagename, band_limits, 
                                               adap_bw, logger)
                    if not freq:
                        logger.warning('Will use fixed frequency bands instead.')
                        freq = self.frequency
                        flag += 1
                    if not chanset[ch]:
                        logchan = ['(no re-refrencing)']
                    else:
                        logchan = chanset[ch]
                        
                    logger.debug(f"Running detection using frequency bands: "
                                 f"{round(freq[0],2)}-{round(freq[1],2)} Hz for "
                                 f"{sub}, {ses}, {str(ch)}:{'-'.join(logchan)}")    
                    
                    # h. Read data
                    logger.debug(f"Reading EEG data for {sub}, {ses}, {str(ch)}:{'-'.join(logchan)}")
                    try:
                        segments = fetch(dset, annot, cat=cat, chan_full = [ch],
                                         stage=self.stage, 
                                         cycle=cycle, reject_epoch=True, 
                                         reject_artf=self.reject)
                        segments.read_data([ch], ref_chan=chanset[ch], grp_name=self.grp_name)
                    except Exception as error:
                        logger.error(error.args[0])
                        logger.warning(f"Skipping {sub}, {ses}, channel {str(ch)} ... ")
                        flag+=1
                        continue
                    
                    if len(segments) < 1:
                        logger.warning(f"No valid data found for {sub}, {ses}, {str(ch)}:{'-'.join(logchan)}"
                                       f"Skipping ... ")
                        continue
                    
                    ## i. Loop through methods (i.e. Whale it!)
                    for m, meth in enumerate(method):
                        logger.debug(f'Using method: {meth}')
                        
                        
                        # j. Define detection
                        detection = DetectSpindle(meth, frequency=freq, duration=duration)

                        ## k. Run detection and save to Annotations file
                        for s, seg in enumerate(segments):
                            if cat[0] == 1:
                                stage_cycle = f"Stage {self.stage[s]}"  
                            else:
                                stage_cycle = f"Cycle {seg['cycle'][2]}, Stage {seg['stage']}"
                            logger.debug(f'Detecting events in {stage_cycle}')
                            spindle = detection(seg['data']) # detect spindles
                            if adap_bands == 'Fixed':
                                evt_name = meth 
                            else:
                                evt_name = f'{meth}_adap'
                            spindle.to_annot(annot, evt_name) # write spindles to annotations file
                            if len(spindle.events) == 0:
                                logger.warning(f'No events detected by {meth} for {sub}, {ses}, {stage_cycle}')    
                                flag += 1
                        # l. Remove any duplicate detected spindles on channel 
                        remove_duplicate_evts(annot, evt_name=evt_name, chan=f'{ch} ({self.grp_name})')
                        
        ### 3. Check completion status and print
        if flag == 0:
            logger.debug('Spindle detection finished without error.')  
        else:
            logger.warning(f'Spindle detection finished with {flag} WARNINGS. See log for details.')
        
        #self.tracking = tracking   ## TO UPDATE - FIX TRACKING
        
        return 
    

    
    def whales(self, method, merge_type, chan, rater, stage, ref_chan, grp_name, keyword,
                     cs_thresh, s_freq = None, 
                     duration= (0.5, 3), evt_out = 'spindle', weights = None,  
                     filetype = '.edf', 
                     logger = create_logger('Detect spindles (WHALES)')):
        

        flag = 0
        logger.info('')
        logger.debug(rf""" Whaling it... 
                     
                                   
                                 'spindles'                 
                                    ":"
                                  ___:____     |"\/"|
                                ,'        `.    \  /
                               |  O        \___/  |
                            ~^~^~^~^~^~^~^~^~^~^~^~^~ 
                                                    
                        
                        Wonambi Heuristic Approach to Locating Elementary Spindles
                        (W.H.A.L.E.S)
                        
                        
                        Combining events using method: {merge_type}
                        Consesus threshold = {cs_thresh}
                        
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
        
        ### 2. Begin loop through dataset
       
        # a. Begin loop through participants
        subs.sort()
        for i, sub in enumerate(subs):
            if not sub in self.tracking['spindle'].keys():
                self.tracking['spindle'][sub] = {}
            # b. Begin loop through sessions
            flag, sessions = load_sessions(sub, self.sessions, self.rec_dir, flag, 
                                     logger, verbose=2)  
            
            for v, ses in enumerate(sessions):
                logger.info('')
                logger.debug(f'Commencing {sub}, {ses}')
                if not ses in self.tracking['spindle'][sub].keys():
                    self.tracking['spindle'][sub][ses] = {} 
                     
                # Get sampling frequency
                if not s_freq:
                    rdir = self.rec_dir + '/' + sub + '/' + ses + '/eeg/'
                    try:
                        edf_file = [x for x in listdir(rdir) if x.endswith(filetype)]
                        s_freq = Dataset(rdir + edf_file[0]).header['s_freq']
                    except:
                        logger.warning(f'No input {filetype} file in {rdir}, cannot obtain sampling frequency for {sub}, {ses}. Skipping...')
                        flag +=1
                        continue
                    
                ## c. Load annotations
                xdir = self.xml_dir + '/' + sub + '/' + ses + '/'
                try:
                    xml_file = [x for x in listdir(xdir) if x.endswith('.xml')]
                    # Copy annotations file before beginning
                    if not path.exists(self.out_dir):
                        mkdir(self.out_dir)
                    if not path.exists(self.out_dir + '/' + sub):
                        mkdir(self.out_dir + '/' + sub)
                    if not path.exists(self.out_dir + '/' + sub + '/' + ses):
                        mkdir(self.out_dir + '/' + sub + '/' + ses)
                    backup = self.out_dir + '/' + sub + '/' + ses + '/'
                    backup_file = (f'{backup}{sub}_{ses}.xml')
                    if not path.exists(backup_file):
                        shutil.copy(xdir + xml_file[0], backup_file)
                    else:
                        logger.warning(f'Annotations file already exists for {sub}, {ses}, new events will be written into the same annotations file. Any existing events labelled {evt_out} will be overwritten.')
                        flag += 1
                except:
                    logger.warning(f' No input annotations file in {xdir}')
                    flag += 1
                    break
                
                # Read annotations file
                annot = Annotations(backup_file, rater_name=self.rater)
                
                ## f. Channel setup 
                pflag = deepcopy(flag)
                flag, chanset = load_channels(sub, ses, self.chan, self.ref_chan,
                                              flag, logger)
                if flag - pflag > 0:
                    logger.warning(f'Skipping {sub}, {ses}...')
                    break
                
                for c, ch in enumerate(chanset):
                    all_events = []
                    for m in method:
                        evts = annot.get_events(name = m, chan = f'{ch} ({grp_name})')
                        if len(evts) == 0:
                            logger.warning(f'No events: {m} found for {sub}, {ses} on {ch}')
                            flag += 1
                        all_events.append(sorted(evts, key=lambda d: d['end']))
                    
                    if all(len(arr) == 0 for arr in all_events):
                        logger.warning(f'No events to merge for {sub}, {ses} on {ch}, skipping...')
                        continue
                    logger.debug(f'Combining events for {sub}, {ses}, {ch}')
                    cons = consensus(all_events, cs_thresh, s_freq, 
                                     min_duration = duration[0],
                                     weights = weights)

                    cons.to_annot(annot, evt_out, chan= f'{ch} ({grp_name})')
                    remove_duplicate_evts(annot, evt_name = evt_out, 
                                          chan=f'{ch} ({self.grp_name})')
                    
        ### 3. Check completion status and print
        if flag == 0:
            logger.debug('Spindle merging (WHALES) finished without error.')  
        else:
            logger.warning(f'Spindle merging (WHALES) finished with {flag} WARNINGS. See log for details.')
        
        return 
    

def _is_number(text):
    try:
        float(text)
        return True
    except ValueError:
        return False


def _load_peak_frequencies(path, target_column="Peak power frequency (Hz)", cluster_column="Peak cluster"):
    header = None
    freq_idx = -1
    cluster_idx = None
    freqs = []
    rows = []
    count_value = None

    with path.open(newline="") as f:
        reader = csv.reader(f)
        for row in reader:
            if not row:
                continue

            if header is None:
                if row[0].strip() == "Count" and len(row) > 1 and _is_number(row[1]):
                    count_value = float(row[1])
                    rows.append(row)
                    continue
                if target_column in row:
                    header = row
                    freq_idx = row.index(target_column)
                    if cluster_column in row:
                        cluster_idx = row.index(cluster_column)
                    rows.append(row)
                else:
                    rows.append(row)
                continue

            rows.append(row)
            if len(row) <= freq_idx:
                continue
            if not row[0].strip().isdigit():
                continue
            freq_text = row[freq_idx].strip()
            if not _is_number(freq_text):
                continue
            freqs.append(float(freq_text))

    if header is None or freq_idx < 0:
        raise ValueError(f"{path} is missing the '{target_column}' column.")
    if not freqs:
        raise ValueError(f"No usable frequency values found in {path}.")

    return header, freq_idx, cluster_idx, freqs, rows, count_value


def _two_peak_threshold(freqs, bins=60, smooth_sigma=1.5):
    freqs_arr = np.asarray(freqs, dtype=float)
    counts, edges = np.histogram(freqs_arr, bins=bins)
    centers = 0.5 * (edges[:-1] + edges[1:])
    smoothed = gaussian_filter1d(counts.astype(float), sigma=smooth_sigma)

    peak_indices, _ = find_peaks(smoothed)
    if len(peak_indices) >= 2:
        top_two = peak_indices[np.argsort(smoothed[peak_indices])[-2:]]
        top_two.sort()
        peak_positions = centers[top_two]
    elif len(peak_indices) == 1:
        peak_positions = np.array([centers[peak_indices[0]], np.median(freqs_arr)])
    else:
        peak_positions = np.array([freqs_arr.min(), freqs_arr.max()])

    low_peak, high_peak = np.sort(peak_positions)
    return float(0.5 * (low_peak + high_peak))


def _append_cluster_column(rows, header_row_idx, freq_idx, cluster_idx, threshold, cluster_column):
    updated = []
    for idx, row in enumerate(rows):
        if idx < header_row_idx:
            updated.append(row)
            continue
        if idx == header_row_idx:
            if cluster_idx is None:
                updated.append(row + [cluster_column])
                cluster_idx = len(row)
            else:
                updated.append(row)
            continue

        if (
            idx > header_row_idx
            and len(row) > freq_idx
            and row
            and row[0].strip().isdigit()
            and _is_number(row[freq_idx])
        ):
            freq_value = float(row[freq_idx])
            cluster = "low_peak" if freq_value < threshold else "high_peak"
            if cluster_idx is None:
                updated.append(row + [cluster])
            else:
                if len(row) <= cluster_idx:
                    row = row + [""] * (cluster_idx - len(row) + 1)
                row[cluster_idx] = cluster
                updated.append(row)
        else:
            if cluster_idx is None:
                updated.append(row + [""])
            else:
                if len(row) <= cluster_idx:
                    row = row + [""] * (cluster_idx - len(row) + 1)
                updated.append(row)
    return updated


def cluster_peaks(csv_dir, target_column="Peak power frequency (Hz)", cluster_column="Peak cluster", bins=60, smooth_sigma=1.5):
    """
    Add a 'Peak cluster' column to every CSV in the provided directory tree by
    locating two peaks in the target column's distribution and labeling events
    as low/high relative to the midpoint between the peaks.
    """
    logger = create_logger("Cluster peaks")
    csv_dir = Path(csv_dir)
    if not csv_dir.exists():
        raise FileNotFoundError(f"{csv_dir} does not exist")

    csv_paths = [p for p in csv_dir.rglob("*.csv") if p.is_file()]
    if not csv_paths:
        logger.warning(f"No CSV files found under {csv_dir}")
        return []

    results = []
    for csv_path in csv_paths:
        try:
            header, freq_idx, cluster_idx, freqs, rows, _ = _load_peak_frequencies(
                csv_path, target_column=target_column, cluster_column=cluster_column
            )
            threshold = _two_peak_threshold(freqs, bins=bins, smooth_sigma=smooth_sigma)
            header_row_idx = rows.index(header)
            updated_rows = _append_cluster_column(rows, header_row_idx, freq_idx, cluster_idx, threshold, cluster_column)
            with csv_path.open("w", newline="") as f:
                writer = csv.writer(f)
                writer.writerows(updated_rows)
            results.append((csv_path, threshold))
            logger.debug(f"Clustered {csv_path} at threshold={threshold:.4f} Hz")
        except Exception as exc:  # noqa: BLE001
            logger.warning(f"Failed to cluster {csv_path}: {exc}")
    return results
      
