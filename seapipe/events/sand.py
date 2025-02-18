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
from numpy import array, multiply
from copy import deepcopy
from ..utils.logs import create_logger
from ..utils.load import load_channels, load_sessions

class SAND:
    
    """ Seapipe Artefact and Noise Detection (S.A.N.D)

        This module runs automated artefact detection with the option of using
        previously published staging algorithms:
            1. YASA (standard deviation)
            2. YASA (covariance)
            3. (More to come..)
        
    """   
    
    def __init__(self, rec_dir, xml_dir, out_dir, eeg_chan, ref_chan,
                 eog_chan, emg_chan, rater = None, grp_name = 'eeg', 
                 subs='all', sessions='all', tracking = None):
        
        self.rec_dir = rec_dir
        self.xml_dir = xml_dir
        self.out_dir = out_dir
        self.eeg_chan = eeg_chan
        self.ref_chan = ref_chan
        self.eog_chan = eog_chan
        self.emg_chan = emg_chan
        self.rater = rater
        self.grp_name = grp_name
        
        self.subs = subs
        self.sessions = sessions
        
        if tracking == None:
            tracking = {}
        self.tracking = tracking


    def detect_artefacts(self, method, filetype = '.edf', 
                               win_size = 5, 
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
                
                # Check if references are the same for each channel
                chans = [[x for x in chanset]]
                ref_chans = {tuple(val) for val in chanset.values()}
                # If not, setup for running per channel
                if len(ref_chans) > 1:
                    logger.debug(f'Channel:reference pairings are unique for {sub}, {ses}. Detecting artefact per channel.')
                    ref_chans = [list(tuple(val)) for val in chanset.values()]
                    chans = chans[0]
                else:
                    ref_chans = [list(tup) for tup in ref_chans]
                

                # d. Load/create for annotations file
                if not path.exists(self.xml_dir + '/' + sub):
                    mkdir(self.xml_dir + '/' + sub)
                if not path.exists(self.xml_dir + '/' + sub + '/' + ses):
                     mkdir(self.xml_dir + '/' + sub + '/' + ses)
                xdir = self.xml_dir + '/' + sub + '/' + ses
                xml_file = f'{xdir}/{sub}_{ses}_eeg.xml'
                if not path.exists(xml_file):
                    dset = Dataset(rdir + edf_file)
                    create_empty_annotations(xml_file, dset)
                    logger.warning(f'No annotations file exists. Creating annotations file for {sub}, {ses} and detecting Artefacts WITHOUT hypnogram.')
                    annot = Annotations(xml_file)
                    hypno = None
                else:
                    logger.debug(f'Annotations file exists for {sub}, {ses}, staging will be used for Artefact detection.')

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
                
                
                for chan, ref in zip(chans, ref_chans):

                    if 'yasa' in method:     
                        
                        ## c. Load recording
                        try:
                            raw = mne.io.read_raw_edf(rdir + edf_file, 
                                                      include = chan + ref,
                                                      preload=True, verbose = False)
                            
                            mne.set_eeg_reference(raw, ref_channels = ref)
                                                  
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
        
                        # Convert to wonambi annotations format
                        evts = []
                        for x in times:
                            evts.append({'name':'Artefact_covar',
                                  'start':x[0]/s_freq,
                                  'end':x[1]/s_freq,
                                  'chan':[''],
                                  'stage':'',
                                  'quality':'Good',
                                  'cycle':''})
                            
                        # Add to annotations file
                        grapho = graphoelement.Graphoelement()
                        grapho.events = evts          
                        grapho.to_annot(annot)
    
    
    
                    elif 'Cross2025' in method:
                        
                        
                        # Convert raw data to array 
                        dset = Dataset(rdir + edf_file)
                        segments = fetch(dset, annot, cat = (1,1,1,1), 
                                         stage=stage)
                        segments.read_data(chan, ref_chan=ref, 
                                           grp_name=self.grp_name)
                        
                        # Create mask
                        
                        
                        # Convert raw data to array    
                        # data = raw.to_data_frame()
                        # inds = [x for x in data if x in chan]
                        # data = data[inds].T
                        # chan_ind = data.index.to_list()
                        
                        
                        
                        dat = transform_signal(data, s_freq, 'high_butter', 
                                         method_opt={'freq':40,
                                                     'order':3})
                        
                        dat2 = transform_signal(dat[0], s_freq, 'moving_covar', 
                                         method_opt = {'dur':5,'step':3},
                                         dat2 = dat[1])
                        
                        threshold = percentile(dat[0], 95)
                        
                        dat2[dat2<threshold] = 0
                        dat2[dat2>threshold] = 1
                        
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
                                                
                        times = detect_above_zero_regions(dat2)
                        
                        # Convert to wonambi annotations format
                        evts = []
                        for x in times:
                            evts.append({'name':'Artefact_Cross',
                                  'start':float(x[0]*3),
                                  'end':float(x[1]*3),
                                  'chan':[''],
                                  'stage':'',
                                  'quality':'Good',
                                  'cycle':''})
                            
                        # Add to annotations file
                        grapho = graphoelement.Graphoelement()
                        grapho.events = evts          
                        grapho.to_annot(annot)
                            
                    else:
                        logger.error("Currently the only method that is functioning is 'yasa_std' or 'yasa_covar.")
                    
                    # ### get cycles
                    # if self.cycle_idx is not None:
                    #     all_cycles = annot.get_cycles()
                    #     cycle = [all_cycles[i - 1] for i in self.cycle_idx if i <= len(all_cycles)]
                    # else:
                    #     cycle = None
                    
                    # ### if event channel only, specify event channels
                    # # 4.d. Channel setup
                    # flag, chanset = load_channels(sub, ses, self.chan, 
                    #                               self.ref_chan, flag, logger)
                    # if not chanset:
                    #     flag+=1
                    #     break
                    # newchans = rename_channels(sub, ses, self.chan, logger)
    
                    # # get segments
                    # for c, ch in enumerate(chanset):
                    #     logger.debug(f"Reading data for {ch}:{'/'.join(chanset[ch])}")
                    #     segments = fetch(dset, annot, cat = cat,  
                    #                      stage = self.stage, cycle=cycle,  
                    #                      epoch = epoch_opts['epoch'], 
                    #                      epoch_dur = epoch_opts['epoch_dur'], 
                    #                      epoch_overlap = epoch_opts['epoch_overlap'], 
                    #                      epoch_step = epoch_opts['epoch_step'], 
                    #                      reject_epoch = epoch_opts['reject_epoch'], 
                    #                      reject_artf = epoch_opts['reject_artf'],
                    #                      min_dur = epoch_opts['min_dur'])
                    
        return