# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 10:34:53 2021

@author: nathancross
"""
from copy import deepcopy
from datetime import datetime, date
from itertools import product
from numpy import (asarray, ndarray, sum)
from os import listdir, mkdir, path, walk
from pandas import DataFrame, read_csv
import shutil
from wonambi import Dataset
from wonambi.attr import Annotations
from wonambi.detect import consensus, DetectSpindle
from wonambi.trans import (fetch, get_times)
from wonambi.trans.analyze import event_params, export_event_params
from ..utils.logs import create_logger, create_logger_outfile, create_logger_empty
from ..utils.load import (load_channels, load_adap_bands, rename_channels, read_manual_peaks)
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

    def whale_it(self, method, cat, cycle_idx=None, adap_bands=False, peaks=None, 
                 adap_bw=4, duration=(0.5, 3),filetype='.edf', 
                 outfile='detect_spindles_log.txt'):
        
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
        tracking = self.tracking
        if outfile == True:
            evt_out = '_'.join(method)
            today = date.today().strftime("%Y%m%d")
            now = datetime.now().strftime("%H:%M:%S")
            logfile = f'{self.log_dir}/detect_spindles_{evt_out}_{today}_log.txt'
            logger = create_logger_outfile(logfile=logfile, name='Detect spindles')
            logger.info('')
            logger.info(f"-------------- New call of 'Detect spindles' evoked at {now} --------------")
        elif outfile:
            logfile = f'{self.log_dir}/{outfile}'
            logger = create_logger_outfile(logfile=logfile, name='Detect spindles')
        else:
            logger = create_logger('Detect spindles')
        
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
            tracking[f'{sub}'] = {}
            # b. Begin loop through sessions
            sessions = self.sessions
            if sessions == 'all':
                sessions = listdir(self.rec_dir + '/' + sub)
                sessions = [x for x in sessions if not '.' in x]   
            
            for v, ses in enumerate(sessions):
                logger.info('')
                logger.debug(f'Commencing {sub}, {ses}')
                tracking[f'{sub}'][f'{ses}'] = {'spindle':{}} 
    
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
                    backup_file = (f'{backup}{sub}_{ses}_spindle.xml')
                    if not path.exists(backup_file):
                        shutil.copy(xdir + xml_file[0], backup_file)
                    else:
                        logger.debug(f'Annotations file already exists for {sub}, {ses}, any previously detected events will be overwritten.')
                    annot = Annotations(backup_file, rater_name=self.rater)
                except:
                    logger.warning(f' No input annotations file in {xdir}')
                    break
                
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
                        freq = read_manual_peaks(sub, ses, peaks, ch, 
                                                 adap_bw, logger)
                    elif adap_bands == 'Auto':
                        stagename = '-'.join(self.stage)
                        bandwidth = f'{self.frequency[0]}-{self.frequency[1]}Hz'
                        freq = load_adap_bands(self.tracking['fooof'], sub, ses,
                                               fnamechan, stagename, bandwidth, 
                                               adap_bw, logger)
                    if not freq:
                        logger.warning('Will use fixed frequency bands instead.')
                        freq = self.frequency
                    logger.debug(f"Running detection using frequency bands: {round(freq[0],2)}-{round(freq[1],2)} Hz for {sub}, {ses}, {str(ch)}:{'-'.join(chanset[ch])}")    
                    
                    # h. Read data
                    logger.debug(f"Reading EEG data for {sub}, {ses}, {str(ch)}:{'-'.join(chanset[ch])}")
                    try:
                        segments = fetch(dset, annot, cat=cat, stage=self.stage, 
                                         cycle=cycle, reject_epoch=True, 
                                         reject_artf=['Artefact', 'Arou', 'Arousal'])
                        segments.read_data([ch], ref_chan=chanset[ch], grp_name=self.grp_name)
                    except Exception as error:
                        logger.error(type(error).__name__, "–", error)
                        flag+=1
                        break
    
                    ## i. Loop through methods (i.e. Whale it!)
                    for m, meth in enumerate(method):
                        logger.debug(f'Using method: {meth}')
                        
                        
                        # j. Define detection
                        detection = DetectSpindle(meth, frequency=freq, duration=duration)

                        ## k. Run detection and save to Annotations file
                        if cat[0] == 1 and cat[1] == 0:
                            for s, seg in enumerate(segments):
                                logger.debug(f'Detecting events in stage {self.stage[s]}')
                                spindle = detection(seg['data']) # detect spindles
                                if adap_bands == 'Fixed':
                                    evt_name = meth 
                                else:
                                    evt_name = f'{meth}_adap'
                                spindle.to_annot(annot, evt_name) # write spindles to annotations file
                                if len(spindle.events) == 0:
                                    logger.warning(f'No events detected by {meth} for {sub}, {ses}')    
                            now = datetime.now().strftime("%m-%d-%Y, %H:%M:%S")
                            tracking[f'{sub}'][f'{ses}']['spindle'][f'{ch}'] = {'Method':meth,
                                                                              'Stage':self.stage,
                                                                              'Cycle':'All',
                                                                              'File':backup_file,
                                                                              'Updated':now}
                        else:
                            for s, seg in enumerate(segments):
                                logger.debug('Detecting events in cycle {} of {}, stages: {}'.format(s + 1, 
                                      len(segments),self.stage))
                                spindle = detection(seg['data'])
                                spindle.to_annot(annot, meth)
                                if len(spindle.events) == 0:
                                    logger.warning(f'No events detected by {meth} for {sub}, {ses}')
                            now = datetime.now().strftime("%m-%d-%Y, %H:%M:%S")
                            tracking[f'{sub}'][f'{ses}']['spindle'][f'{ch}'] = {'Method':meth,
                                                                              'Stage':self.stage,
                                                                              'Cycle':list(range(1,len(segments))),
                                                                              'File':backup_file,
                                                                              'Updated':now}
                        
                        # l. Remove any duplicate detected spindles on channel 
                        remove_duplicate_evts(annot, evt_name=meth, chan=f'{ch} ({self.grp_name})')
                        
        ### 3. Check completion status and print
        if flag == 0:
            logger.debug('Spindle detection finished without ERROR.')  
        else:
            logger.warning('Spindle detection finished with WARNINGS. See log for details.')
        
        #self.tracking = tracking   ## TO UPDATE - FIX TRACKING
        
        return 
    
    
    def whales(self, out_dir, method, chan, rater, cat, stage, ref_chan, grp_name, keyword,
               cs_thresh, s_freq, min_duration, frequency=(11, 16), duration= (0.5, 3),  
                 part='all', visit='all', evt_type='spindle', weights=None,
                 outfile='export_params_log.txt'):
        
        
        ### 0. Set up logging
        logger = create_logger('Detect spindles')
        tracking = self.tracking
        flag = 0
        
        # loop through records
        if isinstance(part, list):
            None
        elif part == 'all':
                part = listdir(out_dir + '/')
                part = [ p for p in part if not '.' in p]
        else:
            print('')
            print("ERROR: 'part' must either be an array of subject ids or = 'all' ")
            print('')
            
        logger.debug(r"""Whaling it... 
                                   
                             'spindles'                 
                                ":"
                             ___:____    |"\/"|
                           ,'        `.   \  /
                          |  O        \___/  |
                        ~^~^~^~^~^~^~^~^~^~^~^~^~ 
                                                    """,)

                
        for i, p in enumerate(part):
            if visit == 'all':
                visit = listdir(out_dir + '/' + p)
                visit = [x for x in visit if not '.' in x]
            for v, vis in enumerate(visit): 
                if not path.exists(out_dir + '/'+ p + '/' + vis + '/'):
                    print(f'WARNING: whale_it has not been run for Subject {p}, skipping..')
                    continue
                elif not path.exists(out_dir + '/' + p + '/' + vis + r'/consensus/'):
                    mkdir(out_dir + '/' + p + '/' + vis + r'/consensus/')
                backup_dir = out_dir + '/'+ p + '/' + vis + r'/consensus'
                xml_file = [x for x in listdir(out_dir + '/'+ p + '/' + vis) if x.endswith('.xml')] 
                for x, file in enumerate(xml_file):
                    pre = file.split(".")[0]
                    ext = file.split(".")[1]
                    backup_file = (f'{backup_dir}/{pre}_{keyword}.{ext}')
                    orig_file = out_dir + '/' + p + '/' + vis + '/' + file
                    shutil.copy(orig_file, backup_file)
            
                    # Create consensus events and export to XML    
                    annot = Annotations(backup_file, rater_name=rater)
                    for c, ch in enumerate(chan):
                        print(ch)
                        all_events = []
                        for m in method:
                            all_events.append(annot.get_events(name=m, chan=ch + ' (' + grp_name + ')'))
                        print('Coming to a consensus for ' + file + ', ' + "channel '" + ch + "'...")
                        cons = consensus(all_events, cs_thresh, s_freq, min_duration=min_duration,
                                         weights=weights)
                        cons.to_annot(annot, evt_type, chan= ch + ' (' + grp_name + ')') # New consensus XML
     
        return 
    


      