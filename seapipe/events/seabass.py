#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 10 18:40:35 2024

@author: ncro8394
"""


from os import listdir, mkdir, path
from wonambi import Dataset 
from wonambi.attr import Annotations, create_empty_annotations
import mne
import yasa
from copy import deepcopy
from ..utils.logs import create_logger
from ..utils.load import load_sessions, load_stagechan, load_emg, load_eog
from ..utils.misc import infer_eeg, infer_eog, infer_emg

class seabass:
    
    """ Sleep Events Analysis Basic Automated Sleep Staging (S.E.A.B.A.S.S)

        This module runs automated sleep staging with the option of using
        previously published staging algorithms:
            1. Vallat et al. (2020) - YASA
            2. 
        
    """   
    
    def __init__(self, rec_dir, out_dir, eeg_chan, ref_chan,
                 eog_chan, emg_chan, rater = None, subs='all', sessions='all', 
                 tracking = None):
        
        self.rec_dir = rec_dir
        self.out_dir = out_dir
        self.eeg_chan = eeg_chan
        self.ref_chan = ref_chan
        self.eog_chan = eog_chan
        self.emg_chan = emg_chan
        self.rater = rater
        
        self.subs = subs
        self.sessions = sessions
        
        if tracking == None:
            tracking = {}
        self.tracking = tracking

    def detect_stages(self, method, qual_thresh = 0.5, invert = False, 
                            filetype = '.edf', 
                            logger = create_logger('Detect sleep stages')):
        
        ''' Automatically detects sleep stages by applying a published 
            prediction algorithm.
        
            Creates a new annotations file if one doesn't already exist.
        
            Parameters
            ----------
            
            method      ->   str of name of automated detection algorithm to 
                             detect staging with. Currently only 'Vallat2021' 
                             is supported. 
                             (https://doi.org/10.7554/eLife.70092)
                             
            qual_thresh ->   Quality threshold. Any stages with a confidence of 
                             prediction lower than this threshold will be set 
                             to 'Undefined' for futher manual review.
   
        
        '''
        
        flag = 0
        tracking = self.tracking
        
        logger.info('')
        logger.debug(rf"""Commencing sleep stage detection... 
                     
                     
                                  /`·.¸
                                 /¸...;..¸¸:·
                             ¸.·´  ¸       `'·.¸.·´)
                            : © ):´;          ¸    )
                             `·.¸ `·      ¸.·\ ´`·¸)
                                 `\\``''´´\¸.'
                                
                                
                    Sleep Events Analysis Basic Automated Sleep Staging 
                    (S.E.A.B.A.S.S.)
                    
                    Using method: {method}
                    
                                                    """,)
        ### 1. First we check the directories
        # a. Check for output folder, if doesn't exist, create
        if path.exists(self.out_dir):
                logger.debug(f"Output directory: {self.out_dir} exists")
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
                
                # Load EEG
                dset = Dataset(rdir + edf_file)
                pflag = deepcopy(flag)
                flag, chanset = load_stagechan(sub, ses, self.eeg_chan, self.ref_chan,
                                              flag, logger)
                if flag - pflag > 0:
                    logger.debug('Inferring EEG from recording instead...')
                    eeg_chan = infer_eeg(dset, logger)
                    ref_chan = []
                    if not eeg_chan:
                        logger.warning(f'Skipping {sub}, {ses}...')
                        break
                else:    
                    eeg_chan = [x for x in chanset]
                    ref_chan = chanset[eeg_chan[0]]
                
                # Load EMG
                pflag = deepcopy(flag)
                flag, emg_chan = load_emg(sub, ses, self.emg_chan, 
                                    flag, logger)
                if flag - pflag > 0:
                    logger.debug('Inferring EMG from recording instead...')
                    emg_chan = infer_emg(dset, logger)
                    if not emg_chan:
                        logger.warning(f'Skipping {sub}, {ses}...')
                        break
                
                # Load EOG
                pflag = deepcopy(flag)
                flag, eog_chan = load_eog(sub, ses, self.eog_chan,
                                    flag, logger)
                if flag - pflag > 0:
                    logger.debug('Inferring EOG from recording instead...')
                    eog_chan = infer_eog(dset, logger)
                    if not eog_chan:
                        logger.warning(f'Skipping {sub}, {ses}...')
                        break
                              
                if not isinstance(eeg_chan, list):
                    eeg_chan = [eeg_chan]
                if not isinstance(ref_chan, list):
                    ref_chan = [ref_chan]
                if not isinstance(eog_chan, list):
                    eog_chan = [eog_chan]
                if not isinstance(emg_chan, list):
                    emg_chan = [emg_chan]
                
                ## c. Load recording
                try:
                    chans = eeg_chan + ref_chan + eog_chan + emg_chan
                    chans = [x for x in chans if x]
                    raw = mne.io.read_raw_edf(rdir + edf_file, 
                                              include = chans,
                                              preload=True, verbose = False)
                except Exception as e:
                    logger.warning(f'Error loading {filetype} file in {rdir}, {repr(e)}')
                    flag += 1
                    break
                
                # d. Load/create for annotations file
                if not path.exists(f'{self.out_dir}/{sub}'):
                    mkdir(f'{self.out_dir}/{sub}')
                if not path.exists(f'{self.out_dir}/{sub}/{ses}'):
                     mkdir(f'{self.out_dir}/{sub}/{ses}')
                xdir = f'{self.out_dir}/{sub}/{ses}'
                xml_file = f'{xdir}/{sub}_{ses}_eeg.xml'
                if not path.exists(xml_file):
                    create_empty_annotations(xml_file, dset)
                    logger.debug(f'Creating annotations file for {sub}, {ses}')
                else:
                    logger.warning(f'Annotations file exists for {sub}, {ses}, staging will be overwritten.')
                    flag += 1
                annot = Annotations(xml_file)
                
                
                if method == 'Vallat2021':
                    logger.debug(f'Predicting sleep stages file for {sub}, {ses}')
                    if len(eeg_chan) > 1:
                        logger.warning(f'Method: {method} only takes 1 eeg channel, but {len(eeg_chan)} were given. Skipping {sub}, {ses}...')
                        break
                    epoch_length = 30
                    stage_key = {'W': 'Wake',
                                 'N1': 'NREM1',
                                 'N2': 'NREM2',
                                 'N3': 'NREM3',
                                 'R': 'REM'}
                    if len([x for x in ref_chan if x]) > 0:
                        raw.set_eeg_reference(ref_channels=ref_chan, 
                                          verbose = False)
                    if len(emg_chan) < 1:
                        emg_chan = [None]
                    if len(eog_chan) < 1:
                        eog_chan = [None]    
                    sls = yasa.SleepStaging(raw, 
                                            eeg_name=eeg_chan[0], 
                                            eog_name=eog_chan[0],
                                            emg_name=emg_chan[0])
                    hypno = sls.predict()
                    proba = sls.predict_proba()
                    
                else:
                    logger.critical("Currently 'Vallat2021' is the only supported method.")
                    return
                
                # Save staging to annotations
                if method not in annot.raters:
                    annot.add_rater(method)

                idx_epoch = 0
                for i, key in enumerate(hypno):
                    epoch_beg = 0 + (idx_epoch * epoch_length)
                    one_stage = stage_key[key]
                    annot.set_stage_for_epoch(epoch_beg, one_stage,
                                             attr='stage',
                                             save=False)
                    
                    if proba[key][i] < qual_thresh:
                        annot.set_stage_for_epoch(epoch_beg, 'Undefined',
                                                 attr='stage',
                                                 save=False)
                    idx_epoch += 1

                annot.save()
        
        ### 3. Check completion status and print
        if flag == 0:
            logger.debug('Sleep stage detection finished without error.')  
        else:
            logger.warning(f'Sleep stage detection finished with {flag} WARNINGS. See log for details.')
        
        #self.tracking = tracking   ## TO UPDATE - FIX TRACKING
        
        return 
    
    
    
    
                    
                    
                    