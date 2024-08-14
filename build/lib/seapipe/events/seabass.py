#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 10 18:40:35 2024

@author: ncro8394
"""


from os import listdir, mkdir, path
import shutil
from wonambi import Dataset 
from wonambi.attr import Annotations, create_empty_annotations
from wonambi.detect import DetectSlowWave
from wonambi.trans import fetch
import mne
import yasa
from xml.etree.ElementTree import Element, SubElement, tostring, parse

from copy import deepcopy
from datetime import datetime, date
from pandas import DataFrame
from ..utils.logs import create_logger, create_logger_outfile
from ..utils.load import (load_channels, read_inversion, rename_channels)
from ..utils.misc import remove_duplicate_evts
 

class seabass:
    
    """ Sleep Events Analysis Basic Automated Sleep Staging (SEABASS)

        This module runs automated sleep staging with the option of using
        previously published SO detectors:
            1. Vallat et al. (2020) - YASA
            2. 
        
    """   
    
    def __init__(self, rec_dir, xml_dir, out_dir, log_dir, eeg_chan, ref_chan,
                 eog_chan, emg_chan, rater = None, subs='all', sessions='all', 
                 tracking = None, reject_artf = ['Artefact', 'Arou', 'Arousal']):
        
        self.rec_dir = rec_dir
        self.xml_dir = xml_dir
        self.out_dir = out_dir
        self.log_dir = log_dir
        
        self.eeg_chan = eeg_chan
        self.ref_chan = ref_chan
        self.eog_chan = eog_chan
        self.emg_chan = emg_chan
        self.rater = rater
        self.reject = reject_artf
        
        self.subs = subs
        self.sessions = sessions
        
        if tracking == None:
            tracking = {}
        self.tracking = tracking

    def detect_stages(self, method, qual_thresh = 0.5, invert = False, 
                            filetype = '.edf', 
                            outfile = 'auto_sleep_staging_log.txt'):
        
        ''' Automatically detects sleep stages by applying a published 
            prediction algorithm.
        
            Creates a new annotations file if one doesn't already exist.
        
        INPUTS:
            
            method      ->   str of name of automated detection algorithm to 
                             detect staging with. Currently only 'Vallat2021' 
                             is supported. 
                             (https://doi.org/10.7554/eLife.70092)
                             
            qual_thresh ->   Quality threshold. Any stages with a confidence of 
                             prediction lower than this threshold will be set 
                             to 'Undefined' for futher manual review.
   
        
        '''
        
        ### 0.a Set up logging
        flag = 0
        tracking = self.tracking
        if outfile == True:
            evt_out = '_'.join(method)
            today = date.today().strftime("%Y%m%d")
            now = datetime.now().strftime("%H:%M:%S")
            logfile = f'{self.log_dir}/detect_slowosc_{evt_out}_{today}_log.txt'
            logger = create_logger_outfile(logfile=logfile, name='Detect sleep stages')
            logger.info('')
            logger.info(f"-------------- New call of 'Detect slow oscillations' evoked at {now} --------------")
        elif outfile:
            logfile = f'{self.log_dir}/{outfile}'
            logger = create_logger_outfile(logfile=logfile, name='Detect sleep stages')
        else:
            logger = create_logger('Detect sleep stages')
        
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
            tracking[f'{sub}'] = {}
            # b. Begin loop through sessions
            sessions = self.sessions
            if sessions == 'all':
                sessions = listdir(self.rec_dir + '/' + sub)
                sessions = [x for x in sessions if not '.' in x]   
            
            for v, ses in enumerate(sessions):
                logger.info('')
                logger.debug(f'Commencing {sub}, {ses}')
                tracking[f'{sub}'][f'{ses}'] = {'slowosc':{}} 
    
                ## c. Load recording
                rdir = self.rec_dir + '/' + sub + '/' + ses + '/eeg/'
                try:
                    edf_file = [x for x in listdir(rdir) if x.endswith(filetype)]
                    raw = mne.io.read_raw_edf(rdir + edf_file[0], 
                                              include = self.eeg_chan + 
                                                        self.ref_chan + 
                                                        self.eog_chan + 
                                                        self.emg_chan,
                                              preload=True, verbose = False)
                except:
                    logger.warning(f' No input {filetype} file in {rdir}')
                    break
                
                # d. Load/create for annotations file
                if not path.exists(self.xml_dir + '/' + sub):
                    mkdir(self.xml_dir + '/' + sub)
                if not path.exists(self.xml_dir + '/' + sub + '/' + ses):
                     mkdir(self.xml_dir + '/' + sub + '/' + ses)
                xdir = self.xml_dir + '/' + sub + '/' + ses
                xml_file = f'{xdir}/{sub}_{ses}_eeg.xml'
                if not path.exists(xml_file):
                    dset = Dataset(rdir + edf_file[0])
                    create_empty_annotations(xml_file, dset)
                    logger.debug(f'Creating annotations file for {sub}, {ses}')
                else:
                    logger.warning(f'Annotations file exists for {sub}, {ses}, staging will be overwritten.')
                annot = Annotations(xml_file)
                
                
                if method == 'Vallat2021':
                    logger.debug(f'Predicting sleep stages file for {sub}, {ses}')
                    epoch_length = 30
                    stage_key = {'W': 'Wake',
                                 'N1': 'NREM1',
                                 'N2': 'NREM2',
                                 'N3': 'NREM3',
                                 'R': 'REM'}
                    raw.set_eeg_reference(ref_channels=self.ref_chan, 
                                          verbose = False)
                    sls = yasa.SleepStaging(raw, 
                                            eeg_name=self.eeg_chan[0], 
                                            eog_name=self.eog_chan[0],
                                            emg_name=self.emg_chan[0])
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
        return
    
    
    
    def detect_artefacts(self, method, qual_thresh = 0.5, invert = False, 
                               cat = (1,1,1,1), filetype = '.edf', 
                               outfile = 'artefact_detection_log.txt'):
        
        ''' Automatically detects sleep stages by applying a published 
            prediction algorithm.
        
            Creates a new annotations file if one doesn't already exist.
        
        INPUTS:
            
            method      ->   str of name of automated detection algorithm to 
                             detect staging with. Currently only 'Vallat2021' 
                             is supported. 
                             (https://doi.org/10.7554/eLife.70092)
                             
            qual_thresh ->   Quality threshold. Any stages with a confidence of 
                             prediction lower than this threshold will be set 
                             to 'Undefined' for futher manual review.
   
        
        '''
        
        ### 0.a Set up logging
        flag = 0
        tracking = self.tracking
        if outfile == True:
            evt_out = '_'.join(method)
            today = date.today().strftime("%Y%m%d")
            now = datetime.now().strftime("%H:%M:%S")
            logfile = f'{self.log_dir}/detect_slowosc_{evt_out}_{today}_log.txt'
            logger = create_logger_outfile(logfile=logfile, name='Detect artefacts')
            logger.info('')
            logger.info(f"-------------- New call of 'Detect slow oscillations' evoked at {now} --------------")
        elif outfile:
            logfile = f'{self.log_dir}/{outfile}'
            logger = create_logger_outfile(logfile=logfile, name='Detect artefacts')
        else:
            logger = create_logger('Detect artefacts')
        
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
            tracking[f'{sub}'] = {}
            # b. Begin loop through sessions
            sessions = self.sessions
            if sessions == 'all':
                sessions = listdir(self.rec_dir + '/' + sub)
                sessions = [x for x in sessions if not '.' in x]   
            
            for v, ses in enumerate(sessions):
                logger.info('')
                logger.debug(f'Commencing {sub}, {ses}')
                tracking[f'{sub}'][f'{ses}'] = {'slowosc':{}} 
    
                ## c. Load recording
                rdir = self.rec_dir + '/' + sub + '/' + ses + '/eeg/'
                try:
                    edf_file = [x for x in listdir(rdir) if x.endswith(filetype)]
                    raw = mne.io.read_raw_edf(rdir + edf_file[0], 
                                              include = self.eeg_chan + 
                                                        self.ref_chan + 
                                                        self.eog_chan + 
                                                        self.emg_chan,
                                              preload=True, verbose = False)
                except:
                    logger.warning(f' No input {filetype} file in {rdir}')
                    break
                
                # d. Load/create for annotations file
                if not path.exists(self.xml_dir + '/' + sub):
                    mkdir(self.xml_dir + '/' + sub)
                if not path.exists(self.xml_dir + '/' + sub + '/' + ses):
                     mkdir(self.xml_dir + '/' + sub + '/' + ses)
                xdir = self.xml_dir + '/' + sub + '/' + ses
                xml_file = f'{xdir}/{sub}_{ses}_eeg.xml'
                if not path.exists(xml_file):
                    dset = Dataset(rdir + edf_file[0])
                    create_empty_annotations(xml_file, dset)
                    logger.debug(f'Creating annotations file for {sub}, {ses}')
                else:
                    logger.warning(f'Annotations file exists for {sub}, {ses}, staging will be overwritten.')
                annot = Annotations(xml_file)
                
            
                ### get cycles
                if self.cycle_idx is not None:
                    all_cycles = annot.get_cycles()
                    cycle = [all_cycles[i - 1] for i in self.cycle_idx if i <= len(all_cycles)]
                else:
                    cycle = None
                
                ### if event channel only, specify event channels
                # 4.d. Channel setup
                flag, chanset = load_channels(sub, ses, self.chan, 
                                              self.ref_chan, flag, logger)
                if not chanset:
                    flag+=1
                    break
                newchans = rename_channels(sub, ses, self.chan, logger)

                # get segments
                for c, ch in enumerate(chanset):
                    logger.debug(f"Reading data for {ch}:{'/'.join(chanset[ch])}")
                    segments = fetch(dset, annot, cat = cat,  
                                     stage = self.stage, cycle=cycle,  
                                     epoch = epoch_opts['epoch'], 
                                     epoch_dur = epoch_opts['epoch_dur'], 
                                     epoch_overlap = epoch_opts['epoch_overlap'], 
                                     epoch_step = epoch_opts['epoch_step'], 
                                     reject_epoch = epoch_opts['reject_epoch'], 
                                     reject_artf = epoch_opts['reject_artf'],
                                     min_dur = epoch_opts['min_dur'])
                    