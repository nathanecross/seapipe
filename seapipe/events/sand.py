#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 18:13:45 2024

@author: ncro8394
"""


from os import listdir, mkdir, path
import shutil
from wonambi import Dataset, graphoelement
from wonambi.attr import Annotations, create_empty_annotations
from wonambi.detect import DetectSlowWave
from wonambi.trans import fetch
import mne
from scipy.signal import find_peaks, peak_widths
import yasa
from xml.etree.ElementTree import Element, SubElement, tostring, parse
from numpy import array, multiply, nan
from copy import deepcopy
from datetime import datetime, date
from pandas import DataFrame
from ..utils.logs import create_logger, create_logger_outfile
from ..utils.load import (load_channels, load_sessions, load_stagechan, load_emg, 
                          load_eog, read_inversion, rename_channels)
from ..utils.misc import remove_duplicate_evts
 

class SAND:
    
    """ Seapipe Artefact and Noise Detection (S.A.N.D)

        This module runs automated artefact detection with the option of using
        previously published staging algorithms:
            1. YASA (standard deviation)
            2. YASA (covariance)
            3. (More to come..)
        
    """   
    
    def __init__(self, rec_dir, xml_dir, out_dir, log_dir, eeg_chan, ref_chan,
                 eog_chan, emg_chan, rater = None, subs='all', sessions='all', 
                 tracking = None):
        
        self.rec_dir = rec_dir
        self.xml_dir = xml_dir
        self.out_dir = out_dir
        self.log_dir = log_dir
        
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


    def detect_artefacts(self, method, invert = False, filetype = '.edf', 
                               win_size = 5,
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
            logger.info(f"-------------- New call of 'Detect artefacts' evoked at {now} --------------")
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
                
                # Load chans
                pflag = deepcopy(flag)
                flag, chanset = load_channels(sub, ses, self.eeg_chan, self.ref_chan,
                                              flag, logger)
                if flag - pflag > 0:
                    logger.warning(f'Skipping {sub}, {ses}...')
                    break
                eeg_chan = [x for x in chanset]
                ref_chan = chanset[eeg_chan[0]]
                flag, emg_chan = load_emg(sub, ses, self.eeg_chan, 
                                    flag, logger)
                flag, eog_chan = load_eog(sub, ses, self.eeg_chan,
                                    flag, logger)
                if not isinstance(eeg_chan, list):
                    eeg_chan = [eeg_chan]
                if not isinstance(ref_chan, list):
                    ref_chan = [ref_chan]
                if not isinstance(eog_chan, list):
                    eog_chan = [eog_chan]
                if not isinstance(emg_chan, list):
                    emg_chan = [emg_chan]
                
                ## c. Load recording
                rdir = self.rec_dir + '/' + sub + '/' + ses + '/eeg/'
                try:
                    edf_file = [x for x in listdir(rdir) if x.endswith(filetype)]
                    chans = eeg_chan + ref_chan + eog_chan + emg_chan
                    chans = [x for x in chans if x]
                    raw = mne.io.read_raw_edf(rdir + edf_file[0], 
                                              include = chans,
                                              preload=True, verbose = False)
                except:
                    logger.warning(f' No input {filetype} file in {rdir}')
                    flag += 1
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
                                 'Unknown':0}
                    
                    hypno = array([int(stage_key[x]) for x in hypno])
                    sf = raw.info["sfreq"]
                
                if 'yasa' in method:
                    if 'covar' in method:
                        yasa_meth = 'covar'
                    else:
                        yasa_meth = 'std'
                        
                    # Convert raw data to array    
                    data = raw.to_data_frame()
                    inds = [x for x in data if x in eeg_chan]
                    data = data[inds].T
                    data = data.to_numpy()
                    
                    # Upsample hypnogram to match raw data
                    hypno_up = yasa.hypno_upsample_to_data(hypno, 1/30, data, 
                                                           sf_data=sf)
                    
                    # Detect artefacts
                    art, zscores = yasa.art_detect(data, sf, window = win_size, 
                                                           hypno = hypno_up, 
                                                           include = (1, 2, 3, 4), 
                                                           method = yasa_meth, 
                                                           threshold = 3, 
                                                           n_chan_reject = 2, 
                                                           verbose = False)
                    
                    # Upsample artefacts to match raw data
                    art = multiply(art, 1)
                    sf_art = 1/win_size
                    art_up = yasa.hypno_upsample_to_data(art, sf_art, data, sf)
                    
                    # Find start/end times of artefacts
                    peaks = find_peaks(art_up)
                    properties = peak_widths(art_up, peaks[0])
                    times = [x for x in zip(properties[2],properties[3])]
    
                    # Convert to wonambi annotations format
                    evts = []
                    for x in times:
                        evts.append({'name':'Artefact',
                              'start':x[0]/sf,
                              'end':x[1]/sf,
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