#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 11:23:20 2025

@author: ncro8394
"""


from os import listdir, mkdir, path
from pandas import concat, DataFrame
from wonambi import Dataset, graphoelement
from wonambi.attr import Annotations
from wonambi.trans import fetch
import shutil
from yasa import rem_detect
from copy import deepcopy
from ..utils.logs import create_logger
from ..utils.load import load_sessions, load_stagechan, load_emg, load_eog
from ..utils.misc import infer_eog

class remora:
    
    """ Rapid Eye Movement Oscillation Recognition Algorithm (R.E.M.O.R.A)

        This module runs automated sleep staging with the option of using
        previously published staging algorithms:
            1. Vallat et al. (2020) - YASA
            2. 
        
    """   
    
    def __init__(self, rec_dir, xml_dir, out_dir, eog_chan, 
                 ref_chan = None, rater = None, grp_name = 'eeg',
                 reject_artf = ['Artefact', 'Arou', 'Arousal'],
                 subs='all', sessions='all', 
                 tracking = None):
        
        self.rec_dir = rec_dir
        self.xml_dir = xml_dir
        self.out_dir = out_dir
        self.eog_chan = eog_chan
        self.ref_chan = ref_chan
        self.rater = rater
        self.grp_name = grp_name
        self.reject = reject_artf
        self.subs = subs
        self.sessions = sessions
        
        if tracking == None:
            tracking = {}
        self.tracking = tracking
        
        
    def detect_rems(self, method, amplitude = (50, 325), duration = (0.3, 1.5), 
                          stage = ['REM'], cycle_idx = None, filetype = '.edf', 
                          logger = create_logger('Detect sleep stages')):
         
         ''' Automatically detects sleep stages by applying a published 
             prediction algorithm.
         
             Creates a new annotations file if one doesn't already exist.
         
             Parameters
             ----------
             
             method      ->   str of name of automated detection algorithm to 
                              detect staging with. Currently only 'YASA' 
                              is supported. 
                              (https://raphaelvallat.com/yasa/generated/yasa.rem_detect.html)
                              
             qual_thresh ->   Quality threshold. Any stages with a confidence of 
                              prediction lower than this threshold will be set 
                              to 'Undefined' for futher manual review.
    
         
         '''
         
         flag = 0
         tracking = self.tracking
         
         logger.info('')
         logger.debug(rf"""Commencing rapid eye movements detection... 
                      
                      
                                        
                                       /`·.¸       ¸..¸
                     ¸.···´"·........."..···'····'´    '·....¸¸ 
                    >  © ):´;                         ¸.¸¸¸    `````·.¸¸¸·´\ 
                     `·.¸ `·                   ¸.···´´     ````````````. ¸ )     
                         ```''\¸.'´´´´´´´´´´´´                         \  /     
                                                                        '  
                                 
                     Rapid Eye Movement Oscillation Recognition Algorithm 
                     (R.E.M.O.R.A)
                     
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
             return
         
         ### 2. Begin loop through dataset
        
         # a. Begin loop through participants
         subs.sort()
         for i, sub in enumerate(subs):
             # b. Begin loop through sessions
             flag, sessions = load_sessions(sub, self.sessions, self.rec_dir, flag, 
                                      logger, verbose=2)
             for v, ses in enumerate(sessions):
                 logger.info('')
                 logger.debug(f'Commencing {sub}, {ses}')


                 ## c. Load recording
                 rdir = f'{self.rec_dir}/{sub}/{ses}/eeg/'
                 try:
                     edf_file = [x for x in listdir(rdir) if x.endswith(filetype)][0]
                     dset = Dataset(f'{rdir}/{edf_file}/')
                 except:
                     logger.warning(f'No input {filetype} file in {rdir}')
                     break
                 
                 # d. Load/create/read for annotations file
                 xdir = f'{self.xml_dir}/{sub}/{ses}'
                 xml_file = [x for x in listdir(xdir) if x.endswith('.xml')][0]
                 if not path.exists(f'{self.out_dir}/{sub}'):
                     mkdir(f'{self.out_dir}/{sub}')
                 if not path.exists(f'{self.out_dir}/{sub}/{ses}'):
                      mkdir(f'{self.out_dir}/{sub}/{ses}')
                 backup = f'{self.out_dir}/{sub}/{ses}/'
                 backup_file = (f'{backup}{sub}_{ses}.xml')
                 if not path.exists(backup_file):
                     shutil.copy(f'{xdir}/{xml_file}', backup_file)
                 else:
                     logger.warning(f'Annotations file already exists for {sub}, {ses}, any previously detected events will be overwritten.')

                 annot = Annotations(backup_file, rater_name = self.rater)
                 
                 ## e. Get sleep cycles (if any)
                 if cycle_idx is not None:
                     all_cycles = annot.get_cycles()
                     cycle = [all_cycles[y - 1] for y in cycle_idx if y <= len(all_cycles)]
                 else:
                     cycle = None
                 
                 # f. Load chans
                 pflag = deepcopy(flag)
                 flag, eog_chan = load_eog(sub, ses, self.eog_chan,
                                      flag, logger)
                 if flag - pflag > 0 or len(eog_chan) != 2:
                     logger.debug("Unable to obtain EOG from tracking sheet. Ensure there are 2 EOG channels listed under 'eog'.")
                     logger.debug("Inferring EOG directly from file instead.")
                     eog_chan = infer_eog(dset, logger)
                     if not eog_chan:
                          logger.warning(f'Unable to infer EOG. Skipping {sub}, {ses}...')
                          break
                 if not isinstance(eog_chan, list):
                      eog_chan = [eog_chan]
                      
                 # Load ref
                 flag, chanset = load_stagechan(sub, ses, eog_chan, self.ref_chan,
                                               flag, logger)
                 ref_chan = chanset[eog_chan[0]]
                 
                 # h. Read data
                 logger.debug(f"Reading EOG data for {sub}, {ses}, {','.join(eog_chan)}")
                 try:
                    segments = fetch(dset, annot, cat = (0,0,1,1), stage = stage, 
                                     cycle = cycle, reject_epoch = True, 
                                     reject_artf = self.reject)
                    segments.read_data(eog_chan, ref_chan, 
                                       grp_name = self.grp_name)
                 except Exception as error:
                    logger.error(error.args[0])
                    logger.warning(f"Problem loading channel {str(eog_chan)} for {sub}, {ses} ... ")
                    flag+=1
                    continue
                 
                 events_all = DataFrame(data=None)
                 for seg in segments:
                     loc = seg['data'].data[0][0]
                     roc = seg['data'].data[0][1]
                     rem = rem_detect(roc, loc, sf = dset.header['s_freq'], 
                                      hypno = None, include = 4, amplitude = amplitude, 
                                      duration = duration, freq_rem = (0.5, 5), 
                                      remove_outliers = False, verbose = False)
                     
                     events = rem.summary()
                     
                     # Get start of REM segment 
                     seg_start = seg['data'].axis['time'][0][0]
                     
                     rem_evts = []
                     
                     for e in events.index:
                         rem_evts.append({'name':'rem',
                                          'start': seg_start + events.loc[e]['Start'],
                                          'end': seg_start + events.loc[e]['End'],
                                          'chan': [''],
                                          'stage': '', 
                                          'quality': 'Good', 
                                          'cycle': ''})
                         
                         # Update summary with real start and end times too.
                         events.loc[e, 'Start'] = seg_start + events.loc[e]['Start']
                         events.loc[e, 'End'] = seg_start + events.loc[e]['End']
                         
                     grapho = graphoelement.Graphoelement()
                     grapho.events = rem_evts       
                     grapho.to_annot(annot, 'rem')
                     concat((events_all, events))
                     
                 #Save summary to dataset
                 events_all.to_csv(f'{backup}/{sub}_{ses}_rems_summary.csv')
                 
                 print('HOLD')



