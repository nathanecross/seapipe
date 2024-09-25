# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 11:00:11 2021

@author: labottdv
"""

from os import listdir, mkdir, path
import shutil
from wonambi import Dataset 
from wonambi.attr import Annotations
from wonambi.detect import DetectSlowWave
from wonambi.trans import fetch

from copy import deepcopy
from datetime import datetime, date
from pandas import DataFrame
from ..utils.logs import create_logger, create_logger_outfile
from ..utils.load import (load_channels, read_inversion, load_sessions)
from ..utils.misc import remove_duplicate_evts


class seasnakes:
    
    """ Sleep Events Analysis of Slow Neocortical oscillations And K-complexES (SEASNAKES)

        This module runs a slow oscillation detection with the option of using
        previously published SO detectors:
            1. Massimini et al. (2007)
            2. Massimini et al. (2007) [Adapted to AASM criteria]
            3. Ngo et al. (2015)
            4. Staresina et al. (2015)
        
    """   
    
    def __init__(self, rec_dir, xml_dir, out_dir, log_dir, chan, ref_chan, 
                 grp_name, stage, rater = None, subs='all', 
                 sessions='all', tracking = None,
                 reject_artf = ['Artefact', 'Arou', 'Arousal']):
        
        self.rec_dir = rec_dir
        self.xml_dir = xml_dir
        self.out_dir = out_dir
        self.log_dir = log_dir
        
        self.chan = chan
        self.ref_chan = ref_chan
        self.grp_name = grp_name
        self.stage = stage
        self.rater = rater
        self.reject = reject_artf
        
        self.subs = subs
        self.sessions = sessions
        
        if tracking == None:
            tracking = {}
        self.tracking = tracking

    def detect_slowosc(self, method, cat, cycle_idx = None, duration = (0.5, 3), 
                       average_channels = False, invert = False, filetype = '.edf', 
                       outfile = 'detect_slowosc_log.txt'):
        
        ''' Detect slow oscillations
        
        Runs one (or multiple) slow oscillation detection algorithms and saves the 
        detected events to a new annotations file.
        
        INPUTS:
            method ->    List of names of automated detection algorithms to detect 
                         events with. e.g. ['Massimini2004','AASM/Massimini2004',
                                            'Staresina2015', 'Ngo2015']
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
        tracking = self.tracking
        if outfile == True:
            evt_out = '_'.join(method)
            today = date.today().strftime("%Y%m%d")
            now = datetime.now().strftime("%H:%M:%S")
            logfile = f'{self.log_dir}/detect_slowosc_{evt_out}_{today}_log.txt'
            logger = create_logger_outfile(logfile=logfile, name='Detect slow oscillations')
            logger.info('')
            logger.info(f"-------------- New call of 'Detect slow oscillations' evoked at {now} --------------")
        elif outfile:
            logfile = f'{self.log_dir}/{outfile}'
            logger = create_logger_outfile(logfile=logfile, name='Detect slow oscillations')
        else:
            logger = create_logger('Detect slow oscillations')
        
        logger.info('')
        logger.debug(r"""Commencing slow oscillation detection... 
                     
                     
                                           .•. 
                                          •   •.
                              •.         •      •       
                           .•  •        •        •.       .•.
           . .      .•.   •    •       •           •.   .•   •. .       
              •. .•    •.•     •      •              •..         
                               •     •                  
                               •    •                
                               •   •                  
                               •  •           
                               •.•  
                                
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
            flag, sessions = load_sessions(sub, self.sessions, self.rec_dir, flag, 
                                     logger, verbose=2)
            
            for v, ses in enumerate(sessions):
                logger.info('')
                logger.debug(f'Commencing {sub}, {ses}')
                tracking[f'{sub}'][f'{ses}'] = {'slowosc':{}} 
    
                ## c. Load recording
                rdir = self.rec_dir + '/' + sub + '/' + ses + '/eeg/'
                try:
                    edf_file = [x for x in listdir(rdir) if x.endswith(filetype)]
                    dset = Dataset(rdir + edf_file[0])
                except:
                    logger.warning(f' No input {filetype} file in {rdir}')
                    break
                
                ## d. Load annotations
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
                    backup_file = (f'{backup}{sub}_{ses}_slowosc.xml')
                    if not path.exists(backup_file):
                        shutil.copy(xdir + xml_file[0], backup_file)
                    else:
                        logger.warning(f'Annotations file already exists for {sub}, {ses}, any previously detected events will be overwritten.')
                except:
                    logger.warning(f' No input annotations file in {xdir}')
                    break
                
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
                    break
                
                for c, ch in enumerate(chanset):
                    
                    # g. Check if channel needs to be inverted for detection
                    if type(invert) == type(DataFrame()):
                        inversion = read_inversion(sub, ses, invert, ch, logger)
                        if not inversion:
                            logger.warning(f"NO inversion will be applied to channel {ch} prior to detection for {sub}, {ses}. To turn off this warning, select invert = 'False'")
                        else: 
                            logger.debug(f'Inverting channel {ch} prior to detection for {sub}, {ses}')
                    elif type(invert) == bool:
                        inversion = invert
                    else:
                        logger.critical(f"The argument 'invert' must be set to either: 'True', 'False' or 'None'; but it was set as {invert}.")
                        logger.info('Check documentation for how to set up staging data:')
                        logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
                        logger.info('-' * 10)
                        return
                        
                    # h. Read data
                    logger.debug(f"Reading EEG data for {sub}, {ses}, {str(ch)}:{'-'.join(chanset[ch])}")
                    try:
                        segments = fetch(dset, annot, cat=cat, stage=self.stage, 
                                         cycle=cycle, reject_epoch=True, 
                                         reject_artf=self.reject)
                        segments.read_data([ch], ref_chan=chanset[ch], grp_name=self.grp_name,
                                           average_channels=average_channels)
                    except Exception as error:
                        logger.error(error.args[0])
                        logger.warning(f"Skipping {sub}, {ses}, channel {str(ch)} ... ")
                        flag+=1
                        continue
    
                    ## i. Loop through methods (i.e. Whale it!)
                    for m, meth in enumerate(method):
                        logger.debug(f'Using method: {meth}')
                        
                        # j. define detection
                        detection = DetectSlowWave(meth, duration=duration)
                        detection.invert = inversion

                        ## k. Run detection and save to Annotations file
                        if cat[0] == 1 and cat[1] == 0:
                            for s, seg in enumerate(segments):
                                logger.debug(f'Detecting events in stage {self.stage[s]}')
                                evt_name = meth 
                                event = detection(seg['data']) # detect events
                                event.to_annot(annot, evt_name) # write events to annotations file
                                if len(event.events) == 0:
                                    logger.warning(f'No events detected by {meth} for {sub}, {ses}')    
                            now = datetime.now().strftime("%m-%d-%Y, %H:%M:%S")
                            tracking[f'{sub}'][f'{ses}']['slowosc'][f'{ch}'] = {'Method':meth,
                                                                              'Stage':self.stage,
                                                                              'Cycle':'All',
                                                                              'File':backup_file,
                                                                              'Updated':now}
                        else:
                            for s, seg in enumerate(segments):
                                logger.debug('Detecting events in cycle {} of {}, stages: {}'.format(s + 1, 
                                      len(segments),self.stage))
                                evt_name = meth 
                                event = detection(seg['data']) # detect events
                                event.to_annot(annot, evt_name) # write events to annotations file
                                if len(event.events) == 0:
                                    logger.warning(f'No events detected by {meth} for {sub}, {ses}')
                            now = datetime.now().strftime("%m-%d-%Y, %H:%M:%S")
                            tracking[f'{sub}'][f'{ses}']['slowosc'][f'{ch}'] = {'Method':meth,
                                                                              'Stage':self.stage,
                                                                              'Cycle':list(range(1,len(segments))),
                                                                              'File':backup_file,
                                                                              'Updated':now}
                        
                        # l. Remove any duplicate detected slow oscillations on channel 
                        remove_duplicate_evts(annot, evt_name=meth, chan=f'{ch} ({self.grp_name})')
                        
        ### 3. Check completion status and print
        if flag == 0:
            logger.info('')
            logger.debug('Slow oscillation detection finished without ERROR.')  
        else:
            logger.info('')
            logger.warning('Slow oscillation finished with WARNINGS. See log for details.')
        
        #self.tracking = tracking   ## TO UPDATE - FIX TRACKING
        
        return 
    
    
    
    
    ## FOR FUTURE:
    
    # def kcomplex():
        
        
    #     return

    
class swordfish:
        
        '''
            Slow Wave Online Realtime/Rapid Detection For Intermittent Stimulation H_ (SWORDFISH)
        '''
        def __init__(self, rec_dir, xml_dir, out_dir, log_dir, chan, ref_chan, 
                 grp_name, stage, frequency=(11,16), rater = None, subs='all', 
                 sessions='all', tracking = {}):
        
            self.rec_dir = rec_dir
            self.xml_dir = xml_dir
            self.out_dir = out_dir
            self.log_dir = log_dir
            
            self.chan = chan
            self.ref_chan = ref_chan
            self.grp_name = grp_name
            self.stage = stage
            self.frequency = frequency
            self.rater = rater
            
            self.subs = subs
            self.sessions = sessions
            
            self.tracking = tracking
            
            return
    
    
    
    