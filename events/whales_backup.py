# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 10:34:53 2021

@author: nathancross
"""
from copy import deepcopy
from datetime import datetime
from itertools import product
from numpy import (asarray, concatenate, char, empty, nan, ndarray, sum, zeros)
from os import listdir, mkdir, path, walk
from pandas import concat, DataFrame, read_csv
from pickle import dump
import shutil
import sys
import traceback
from wonambi import Dataset
from wonambi.attr import Annotations
from wonambi.detect import consensus, DetectSpindle
from wonambi.trans import (fetch, get_times)
from wonambi.trans.analyze import event_params, export_event_params
from utils.logs import create_logger, create_logger_outfile, create_logger_empty
from utils.load import load_channels
from utils.misc import remove_duplicate_evts


class Whales:
    
    """ Wonambi Heuristic Approach to Locating Elementary Spindles (WHALES)

        This script runs a consensus approach to detecting sleep spindles. While we hope
        to improve detection, and remove biases that occur based on the use of any one 
        spindle detector, this is not a perfect solution.
        The pipeline runs in three stages:
            1. Whale_it: Detect spindles with pre-determined published algorithms (see Documentation).
            2. Whales: Assign 'true' events based upon a pre-set agreement threshold, using a consensus 
               of the events detected independently from step 1. This creates a new event called
               'spindle' in the annotations file.
            3. Whale_farm: Exports the parameters of these 'spindle' events.
    """   
    
    def __init__(self, rec_dir, xml_dir, out_dir, chan, ref_chan, grp_name, stage, 
                 frequency=(11,16), rater = None, subs='all', sessions='all'):
        
        self.rec_dir = rec_dir
        self.xml_dir = xml_dir
        self.out_dir = out_dir
        
        self.chan = chan
        self.ref_chan = ref_chan
        self.grp_name = grp_name
        self.stage = stage
        self.frequency = frequency
        self.rater = rater
        
        self.subs = subs
        self.sessions = sessions
        
        self.tracking = {}

    def whale_it(self, method, cat, cycle_idx=None, adap_bands=False, 
                 duration=(0.5, 3)):
        
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
        ### 0. Set up logging
        logger = create_logger('Detect spindles')
        tracking = self.tracking
        flag = 0
        
        ### 1. First we set up the output directory
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
        logger.debug(r"""Whaling it... 
                                   
                             'spindles'                 
                                ":"
                             ___:____    |"\/"|
                           ,'        `.   \  /
                          |  O        \___/  |
                                                    """,)
        logger.info("                        ~^~^~^~^~^~^~^~^~^~^~^~^~ ")
       
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
                tracking[f'{sub}'][f'{ses}'] = {'spindles':{}} 
    
                ## c. Load recording
                rdir = self.rec_dir + sub + '/' + ses + '/eeg/'
                try:
                    edf_file = [x for x in listdir(rdir) if x.endswith('.edf') or x.endswith('.rec') or x.endswith('.eeg')]
                    dset = Dataset(rdir + edf_file[0])
                except:
                    logger.warning(f' No input recording file: {edf_file}')
                    break
                
                ## d. Load annotations
                xdir = self.xml_dir + sub + '/' + ses + '/'
                try:
                    xml_file = [x for x in listdir(xdir) if x.endswith('.xml')]
                    # Copy annotations file before beginning
                    if not path.exists(self.out_dir):
                        mkdir(self.out_dir)
                    if not path.exists(self.out_dir + sub):
                        mkdir(self.out_dir + sub)
                    if not path.exists(self.out_dir + sub + '/' + ses):
                        mkdir(self.out_dir + sub + '/' + ses)
                    backup = self.out_dir + sub + '/' + ses + '/'
                    backup_file = (f'{backup}{sub}_{ses}_spindles.xml')
                    if not path.exists(backup_file):
                        shutil.copy(xdir + xml_file[0], backup_file)
                    else:
                        logger.debug(f'Annotations file already exists for {sub}, {ses}, any previously detected events will be overwritten.')
                    annot = Annotations(backup_file, rater_name=self.rater)
                except:
                    logger.warning(f' No input annotations file: {xml_file}')
                    break
                
                ## e. Get sleep cycles (if any)
                if cycle_idx is not None:
                    all_cycles = annot.get_cycles()
                    cycle = [all_cycles[y - 1] for y in cycle_idx if y <= len(all_cycles)]
                else:
                    cycle = None
                
                ## f. Channel setup 
                pflag = deepcopy(flag)
                flag, chanset = load_channels(sub, ses, self.chan, self.ref_chan, flag, logger)
                if flag - pflag > 0:
                    logger.warning(f'Skipping {sub}, {ses}...')
                    break
                
                ### g. Select and read data
                for c, ch in enumerate(chanset):
                    logger.debug(f'Reading EEG data for {sub}, {ses}, ' + str(ch))
                    try:
                        segments = fetch(dset, annot, cat=cat, stage=self.stage, 
                                         cycle=cycle, reject_epoch=True, 
                                         reject_artf=['Artefact', 'Arou', 'Arousal'])
                        segments.read_data([ch], ref_chan=chanset[ch], grp_name=self.grp_name)
                    except Exception as error:
                        logger.error(type(error).__name__, "–", error)
                        flag+=1
                        break
    
                    ## h. Loop through methods (i.e. Whale it!)
                    for m, meth in enumerate(method):
                        logger.debug(f'Using method: {meth}')
                        
                        # i. Check for adapted bands
                        if adap_bands is True:
                            freq = self.frequency[ch][sub + '_' + ses]
                            logger.debug(f'Using adapted bands for {sub}: {round(freq[0],2)}-{round(freq[1],2)} Hz')
                        else:
                            freq = self.frequency
                        
                        # j. Define detection
                        detection = DetectSpindle(meth, frequency=freq, duration=duration)

                        ## k. Run detection and save to Annotations file
                        if cat[0] == 1 and cat[1] == 0:
                            for s, seg in enumerate(segments):
                                logger.debug(f'Detecting events in stage {self.stage[s]}')
                                spindles = detection(seg['data']) # detect spindles
                                spindles.to_annot(annot, meth) # write spindles to annotations file
                                if len(spindles.events) == 0:
                                    logger.warning(f'No events detected by {meth} for {sub}, {ses}')    
                            now = datetime.now().strftime("%m-%d-%Y, %H:%M:%S")
                            tracking[f'{sub}'][f'{ses}']['spindles'][f'{ch}'] = {'Method':meth,
                                                                              'Stage':self.stage,
                                                                              'Cycle':'All',
                                                                              'File':backup_file,
                                                                              'Updated':now}
                        else:
                            for s, seg in enumerate(segments):
                                logger.debug('Detecting events in cycle {} of {}, stages: {}'.format(s + 1, 
                                      len(segments),self.stage))
                                spindles = detection(seg['data'])
                                spindles.to_annot(annot, meth)
                                if len(spindles.events) == 0:
                                    logger.warning(f'No events detected by {meth} for {sub}, {ses}')
                            now = datetime.now().strftime("%m-%d-%Y, %H:%M:%S")
                            tracking[f'{sub}'][f'{ses}']['spindles'][f'{ch}'] = {'Method':meth,
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
            logger.warning('Spindle detection finished with ERRORS. See log for details.')
        
        #self.tracking = tracking   ## TO UPDATE - FIX TRACKING
        
        return 
    
    
    def whales(out_dir, method, chan, rater, cat, stage, ref_chan, grp_name, keyword,
               cs_thresh, s_freq, min_duration, frequency=(11, 16), duration= (0.5, 3),  
                 part='all', visit='all', evt_type='spindle', weights=None):
        
        
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
    
    
    def whale_farm(self, keyword = None, evt_name = 'spindle', segs = None, 
                   cycle_idx = None, adap_bands = False, param_keys = None, 
                   exclude_poor = False, reject_artf = ['Artefact', 'Arou', 'Arousal'], 
                   epoch_dur = 30, n_fft_sec = 4, Ngo = {'run':False}, 
                   outfile='whales_log.txt'):
        
        '''
        Extracts event parameters per participant and session.
        
        segs: to extract parameters between certain markers, these need to be defined
                in the Annotations file first. 
                Format should be a list of tuples, with both tags named
                e.g. [('N2_ON','N2_OFF'), ('N3_ON','N3_OFF')]
        '''
        
        ### 0. Set up logging
        if outfile:
            logfile = f'{self.out_dir}/{outfile}'
            logger = create_logger_outfile(logfile=logfile, name='Export params')
        else:
            logger = create_logger('Export params')
        flag = 0

        try:
            logger.info('') 
       
            # Loop through subjects
            subs = self.subs
            if isinstance(subs, list):
                None
            elif subs == 'all':
                    subs = next(walk(self.xml_dir))[1]
            else:
                logger.error("'subs' must either be an array of participant ids or = 'all' ")
                
            # Set base parameters for output dataframes
            params = {}
            sublist = subs
            header = []
            
            if param_keys is None:
                param_keys = ['count', 'density','dur', 'ptp', 'energy', 'peakef'] # Default spindle metrics to extract
            
            if cycle_idx is not None:
                    for m, param in enumerate(param_keys):
                        params[param] = zeros((len(subs),len(self.chan)*len(cycle_idx)*len(self.stage)+len(self.chan)+len(self.stage)))
                        params[param].fill(nan)
            else:
                    for m, param in enumerate(param_keys):
                        params[param] = zeros((len(subs),len(self.chan)*len(self.stage)+len(self.chan)))
                        params[param].fill(nan)
            
            subs.sort()               
            for i, sub in enumerate(subs):
                dat = []

                # Loop through visits
                if self.sessions == 'all':
                    visits = next(walk(f'{self.xml_dir}/{sub}'))[1]
                else:
                    visits = self.sessions
                    
                # Update size of file based on number of visits(sessions)
                if i == 0:
                    for m, param in enumerate(param_keys):
                            params[param] = zeros((len(subs), len(visits)*len(params[param][0])))
                            params[param].fill(nan)
                
                visits.sort()
                for v, vis in enumerate(visits):  
                    logger.debug(f'Extracting {evt_name} parameters for Subject {sub}, visit {vis}..') 
                    
                    # Define files
                    rdir = self.rec_dir + '/' + sub + '/' + vis + '/eeg/'
                    edf_file = [x for x in listdir(rdir) if x.endswith('.edf') or x.endswith('.rec') 
                                or x.endswith('.eeg') if not x.startswith('.')]

                    if keyword is not None:
                        xml_file = [x for x in listdir(self.xml_dir + '/' + sub + '/' + vis) 
                                if x.endswith('.xml') if not x.startswith('.') if keyword in x] 
                    else:
                        xml_file = [x for x in listdir(self.xml_dir +  '/' + sub + '/' + vis ) 
                                if x.endswith('.xml') if not x.startswith('.')]
                        
                    if len(xml_file) == 0:                
                        logger.warning(f'{evt_name} has not been detected for Subject {sub}, visit {vis} - skipping..')
                    else:
                        xml_file_path = self.xml_dir + '/' + sub + '/' + vis + r'/' + xml_file[0]
                        # Open dataset
                        dataset = Dataset(rdir + edf_file[0])
                
                        # Import Annotations file
                        annot = Annotations(xml_file_path, rater_name=self.rater)
                    
                        # Get sleep cycles (if any)
                        if cycle_idx is not None:
                            all_cycles = annot.get_cycles()
                            cycle = [all_cycles[y - 1] for y in cycle_idx if y <= len(all_cycles)]
                        else:
                            cycle = None
                        
                        ## Channel setup
                        flag, chanset = load_channels(sub, vis, self.chan, self.ref_chan, flag, logger)
                        
                        # Run through channels
                        for ch, channel in enumerate(chanset):

                            if Ngo['run'] == False:
                                if self.grp_name is not None:
                                    chan_ful = [f'{channel} ({self.grp_name})']
                                else:
                                    chan_ful = [channel]
                                    logger.debug(f'channel is {chan_ful}')
                            else:
                                chan_ful = Ngo['chan']
                            
                            if adap_bands is True:
                                freq = self.frequency[channel][sub + '_' + vis]
                                logger.debug(f'Using adapted bands for {sub}: {round(freq[0],2)}-{round(freq[1],2)} Hz')
                            else:
                                freq = self.frequency
                            
                            # Create header for output file
                            for m, param in enumerate(param_keys):
                                if i == 0:
                                    header.append(f'{param}_{channel}_wholenight_visit_{vis}')
                                    
                            if segs is not None:
                                for s, seg in enumerate(segs):
                                    for m, param in enumerate(param_keys):
                                                    if i == 0:
                                                        header.append(f'{param}_{channel}_{seg[0]}_visit_{vis}')  
                            else:
                                for s, st in enumerate(self.stage):
                                    for m, param in enumerate(param_keys):
                                                    if i == 0:
                                                        header.append(f'{param}_{channel}_{st}_visit_{vis}')
                                    if cycle_idx is not None:
                                        for cy, cyc in enumerate(cycle_idx):
                                            for m, param in enumerate(param_keys):
                                                    if i == 0:
                                                        header.append('{param}_{channel}_{st}_cycle_{str(cy+1)}_visit_{vis}')
                        
                            
                            ### WHOLE NIGHT ###
                            # Select and read data
                            logger.debug('Reading data for ' + sub + ', visit ' + vis + ' '+ str(channel))
 
                            segments = fetch(dataset, annot, cat=(0,0,0,0), evt_type=[evt_name], cycle=cycle, 
                                             chan_full=chan_ful, reject_epoch=True, 
                                             reject_artf = reject_artf)
                            segments.read_data([channel], chanset[channel], grp_name=self.grp_name)
                            if type(chan_ful)!=str and len(chan_ful)>1:
                                chan_ful = chan_ful[0]
                            
                            if len(segments) == 0:
                                logger.info('')
                                logger.warning(f"Events haven't been detected for {sub}, {vis} on channel {channel}, skipping...")
                                logger.info('')

                            # Calculate event density (whole night)
                            poi = get_times(annot, stage=self.stage, cycle=cycle, chan=[channel], exclude=exclude_poor)
                            total_dur = sum([x[1] - x[0] for y in poi for x in y['times']])
                            evts = annot.get_events(name=evt_name, time=None, chan = chan_ful, stage = self.stage)
                            count = len(evts)
                            density = len(evts) / (total_dur / epoch_dur)
                            logger.info('')
                            logger.debug('----- WHOLE NIGHT -----')
                            logger.debug(f'No. Segments = {len(segments)}, Total duration (s) = {total_dur}')
                            logger.debug(f'Density = {density} per epoch')
                            logger.info('')

                            dat.append(count)
                            dat.append(density)
                            
                            # Set n_fft
                            n_fft = None
                            if segments and n_fft_sec is not None:
                                s_freq = segments[0]['data'].s_freq
                                n_fft = int(n_fft_sec * s_freq)
                            
                            # Export event parameters (whole night)
                            lg = create_logger_empty()
                            data = event_params(segments, params='all', band=freq, n_fft=n_fft)
                            if data:
                                    if not path.exists(self.out_dir + '/' + sub):
                                        mkdir(self.out_dir + '/' + sub)
                                    if not path.exists(self.out_dir + '/' + sub + '/' + vis):
                                        mkdir(self.out_dir + '/' + sub + '/' + vis)
                                    out_subs = self.out_dir + '/' + sub + '/' + vis
                                    data = sorted(data, key=lambda x: x['start'])
                                    outputfile = out_subs + '/' + sub + '_' + vis + '_' + channel + '_' + evt_name + '.csv'
                                    logger.debug('Writing to ' + outputfile)  
                                    export_event_params(outputfile, data, count=len(evts), 
                                                        density=density)
                            else:
                                logger.warning('No valid data found.')
                                
                            for ev in data:
                                ev['ptp'] = ev['ptp']()[0][0]
                                ev['energy'] = list(ev['energy'].values())[0]
                                ev['peakef'] = list(ev['peakef'].values())[0]
                                #ev['minamp'] = data[0]['minamp'].data[0][0]
                                #ev['maxamp'] = data[0]['maxamp'].data[0][0]
                                #ev['rms'] = data[0]['rms'].data[0][0]
                            
                            for m, param in enumerate(param_keys[2:]):
                                dat.append(asarray([x[param] for x in data]).mean())

                            if segs is not None:
                                 
                            ### PER SEGMENT ###
                                for s,seg in enumerate(segs): 
                                    
                                    # Calculate event density (per segment type)
                                    poi = get_times(annot, evt_type=[seg[0],seg[1]], stage=None, 
                                                    chan=[channel], exclude=exclude_poor)
                                    duos=[]
                                    [duos.extend([(poi[0]['times'][x][0],poi[1]['times'][x][1])]) 
                                                         for x,item in enumerate(poi[0]['times'])]
                                    total_dur = sum([x[1] - x[0] for x in duos])
                                    evts =[]
                                    for d in duos:
                                        evts.extend(annot.get_events(name=evt_name, time=d, 
                                                                     chan = chan_ful, stage = None))
                                    count = len(evts)
                                    density = len(evts) / (total_dur / epoch_dur)
                                    logger.info('')
                                    logger.debug(f'----- Segment {seg} -----')
                                    logger.debug(f'No. Segments = {len(segments)}, Total duration (s) = {total_dur}')
                                    logger.debug(f'Density = {density} per epoch')
                                    logger.info('')


                                    dat.append(count)
                                    dat.append(density)
                                    
                                    
                                    # Set n_fft
                                    n_fft = None
                                    if segments and n_fft_sec is not None:
                                        s_freq = segments[0]['data'].s_freq
                                        n_fft = int(n_fft_sec * s_freq)
                                    
                                    # Export event parameters
                                    data = event_params(segments, params='all', band=freq, n_fft=n_fft)
                                    if data:
                                            data = sorted(data, key=lambda x: x['start'])
                                            outputfile = self.out_dir + '/' + sub + '_' + vis + '_' + channel + '_' + evt_name + '.csv'
                                            logger.debug('Writing to ' + outputfile)
                                            with open(logfile, 'a') as f:
                                                print('', file=f) 
                                                print('Writing to ' + outputfile, file=f) 
                                            export_event_params(outputfile, data, count=len(evts), 
                                                                density=density)
                                    else:
                                        logger.warning('No valid data found.')
                                        with open(logfile, 'a') as f:
                                                print('', file=f) 
                                                print('No valid data found.', file=f) 
                                        
                                    for ev in data:
                                        ev['ptp'] = ev['ptp']()[0][0]
                                        ev['energy'] = list(ev['energy'].values())[0]
                                        ev['peakef'] = list(ev['peakef'].values())[0]
                                        #ev['minamp'] = data[0]['minamp'].data[0][0]
                                        #ev['maxamp'] = data[0]['maxamp'].data[0][0]
                                        #ev['rms'] = data[0]['rms'].data[0][0]
                                    
                                    for m, param in enumerate(param_keys[2:]):
                                        dat.append(asarray([x[param] for x in data]).mean())

                            else:
                                    
                             ### PER STAGE ###
                                for s, st in enumerate(self.stage):
                                    
                                    if not isinstance(chan_ful, ndarray):
                                        chan_ful = [chan_ful]
                                    
                                    segments = fetch(dataset, annot, cat=(0,0,0,0), evt_type=[evt_name], 
                                                     stage = [st], cycle=cycle, 
                                                     chan_full=chan_ful, reject_epoch=True, 
                                                     reject_artf = reject_artf, min_dur=0.5)
                                    segments.read_data([channel], chanset[channel], grp_name=self.grp_name)
                                    if len(chan_ful) > 1:
                                        chan_ful = chan_ful[0]
                
                                    # Calculate event density (per stage)
                                    poi = get_times(annot, stage=[st], cycle=cycle, chan=[channel], exclude=exclude_poor)
                                    total_dur = sum([x[1] - x[0] for y in poi for x in y['times']])
                                    evts = annot.get_events(name=evt_name, time=None, chan = chan_ful, stage = st)
                                    count = len(evts)
                                    density = len(evts) / (total_dur / epoch_dur)
                                    logger.info('')
                                    logger.debug(f'---- STAGE {st} ----')
                                    logger.debug(f'No. Segments = {len(segments)}, Total duration (s) = {total_dur}')
                                    logger.debug(f'Density = {density} per epoch')
                                    logger.info('')
                                    with open(logfile, 'a') as f:
                                            print('', file=f)
                                            print(f'---- STAGE {st} ----', file=f)
                                            print(f'No. Segments = {len(segments)}, Total duration (s) = {total_dur}', file=f)
                                            print(f'Density = {density} per epoch', file=f)
                                    dat.append(count)
                                    dat.append(density)
                                    
                                    # Set n_fft
                                    n_fft = None
                                    if segments and n_fft_sec is not None:
                                        s_freq = segments[0]['data'].s_freq
                                        n_fft = int(n_fft_sec * s_freq)
                                    
                                    # Export event parameters (per stage)
                                    data = event_params(segments, params='all', band=freq, n_fft=n_fft)
                                        
                                    for ev in data:
                                        ev['ptp'] = ev['ptp']()[0][0]
                                        ev['energy'] = list(ev['energy'].values())[0]
                                        ev['peakef'] = list(ev['peakef'].values())[0]
            
                                    
                                    for m, param in enumerate(param_keys[2:]):
                                        dat.append(asarray([x[param] for x in data]).mean())
                                    
                                    ### PER CYCLE ###
                                    if cycle_idx is not None: 
                                        try:
                                            for cy, cycc in enumerate(cycle_idx):
                                                cyc = cycle[cy]
                                                if not isinstance(chan_ful, list):
                                                    chan_ful = [chan_ful]
                                                    if len(chan_ful) > 1:
                                                        chan_ful = chan_ful[0] 
                                                segments = fetch(dataset, annot, cat=(0,0,0,0), evt_type=[evt_name], 
                                                                 stage = [st], cycle=[cyc], 
                                                                 chan_full=chan_ful, reject_epoch=True, 
                                                                 reject_artf = reject_artf, min_dur=0.5)
                
                                                
                                                segments.read_data([channel], chanset[channel], grp_name=self.grp_name)
                                                
                                                if isinstance(chan_ful, ndarray):
                                                    if len(chan_ful) > 1:
                                                        chan_ful = chan_ful[0]
            
                                                if len(chan_ful) > 1:
                                                    chan_ful = chan_ful[0]
            
                                                # Calculate event density (per cycle)
                                                poi = get_times(annot, stage=[st], cycle=[cyc], chan=[channel], exclude=exclude_poor)
                                                total_dur = sum([x[1] - x[0] for y in poi for x in y['times']])
                                                evts = annot.get_events(name=evt_name, time=cycle[cy][0:2], chan = chan_ful, stage = st)
                                                count = len(evts)
                                                density = len(evts) / (total_dur / epoch_dur)
                                                logger.info('')
                                                logger.info(f'---- STAGE {st}, CYCLE {cy+1} ----')
                                                logger.debug(f'No. Segments = {len(segments)}, Total duration (s) = {total_dur}')
                                                logger.debug(f'Density = {density} per epoch')
                                                logger.info('')
                                                with open(logfile, 'a') as f:
                                                    print('', file=f)
                                                    print(f'---- STAGE {st}, CYCLE {cy+1} ----', file=f)
                                                    print(f'No. Segments = {len(segments)}, Total duration (s) = {total_dur}', file=f)
                                                    print(f'Density = {density} per epoch', file=f)
                                               
                                                dat.append(count)
                                                dat.append(density)
                                                # if i == 0:
                                                #         header.append('density_' + channel + '_' + st + '_' + str(cy+1) + '_visit_' + vis)
                                                
                                                # Set n_fft
                                                n_fft = None
                                                if segments and n_fft_sec is not None:
                                                    s_freq = segments[0]['data'].s_freq
                                                    n_fft = int(n_fft_sec * s_freq)
                                                
                                                # Export event parameters (per cycle)
                                                data = event_params(segments, params='all', band=freq, n_fft=n_fft)
                                                    
                                                for ev in data:
                                                    ev['ptp'] = ev['ptp']()[0][0]
                                                    ev['energy'] = list(ev['energy'].values())[0]
                                                    ev['peakef'] = list(ev['peakef'].values())[0]
            
                                                
                                                for m, param in enumerate(param_keys[2:]):
                                                    dat.append(asarray([x[param] for x in data]).mean())
                                        except Exception as e:
                                            logger.warning(f'No STAGE {st} in CYCLE {cy+1}')
                                            dat.append(nan)
                                            for m, param in enumerate(param_keys[1:]):
                                                dat.append(nan)
                                                # if i == 0:
                                                #     header.append(param + '_' + channel + '_' + st + '_cycle' + str(cy+1) + '_visit_' + vis)
                                            continue        
                
                if i == 0:
                    output = DataFrame(dat)
                else:
                    output = concat([output,DataFrame(dat)], axis=1)
                    
            
            filler = empty((len(header)-output.shape[0],output.shape[1]))
            filler[:] = nan
            output = concatenate((output, filler))
            output = DataFrame(output)
            output = DataFrame.transpose(output) 
            output.columns=header  
            output.index=sublist                 
            DataFrame.to_csv(output, f'{self.out_dir}/{evt_name}_dataset.csv', sep=',')            
            with open(f'{self.out_dir}/{evt_name}.p', 'wb') as f:
                dump(params, f)
                
            return output
            
        except Exception as e:
                logger.error(e)
                tb = sys.exc_info()[2]
                tbinfo = traceback.format_tb(tb)[0]
                logger.info(tbinfo)
                with open(logfile, 'a') as log:
                    print(tbinfo, file=log)
                    print(e, file=log)
                    
    def whale_farm_from_csvs(self, chan, evt_name, part='all', visit='all'):
    
        rec_dir = self.xml_dir
        out_dir = self.out_dir        
    
        if not path.exists(out_dir):
                mkdir(out_dir)
        
        # Loop through subjects
        subs = self.subs
        if isinstance(part, list):
            None
        elif part == 'all':
                part = next(walk(self.xml_dir))[1]
        else:
            logger.error("'subs' must either be an array of participant ids or = 'all' ")
         
        # Loop through visits
        part.sort()                
        for i, p in enumerate(part):
            # Loop through visits
            if self.sessions == 'all':
                visit = next(walk(f'{self.xml_dir}/{p}'))[1]
            else:
                visit = self.sessions
        
        # Set variable names and combined with visits
        variables = ['Count','Density','Duration_mean (s)','Duration_stdv (s)',
                   'Min_amplitude_mean (uV)','Min_amplitude_stdv (uV)', 'Max_amplitude_mean (uV)',
                   'Max_amplitude_stdv (uV)','Ptp_amplitude_mean (uV)', 'Ptp_amplitude_stdev (uV)',
                   'Power_mean (uV^2)','Power_stdev (uV^2)', 'Peak_power_frequency_mean (Hz)',
                   'Peak_power_frequency_std (Hz)']
        
        chanvar = []
        for pair in product(chan, variables):
            chanvar.append('_'.join(pair))
        
        columns = []    
        for pair in product(visit, chanvar):
            columns.append('_'.join(pair))
        
        # Create empty output directory
        dataset = zeros((len(part), len(visit)*len(chan)*14))
        
        
        for p, pp in enumerate(part):
            
            for v, vis in enumerate(visit):  
                print(f'Extracting {evt_name} parameters for Subject {pp}, visit {vis}..')
                
                for c, ch in enumerate(chan):
                    
                    data_file = rec_dir + '/' + pp + '/' + vis + f'/{pp}_{vis}_{ch}_{evt_name}.csv' 
                    
                    if path.isfile(data_file):
                        print(f'Extracting from: {data_file}')
                        # Delimiter
                        data_file_delimiter = ','
                        
                        # The max column count a line in the file could have
                        largest_column_count = 0
                        
                        # Loop the data lines
                        with open(data_file, 'r') as temp_f:
                            # Read the lines
                            lines = temp_f.readlines()
                        
                            for l in lines:
                                # Count the column count for the current line
                                column_count = len(l.split(data_file_delimiter)) + 1
                                
                                # Set the new most column count
                                largest_column_count = column_count if largest_column_count < column_count else largest_column_count
                        
                        # Generate column names (will be 0, 1, 2, ..., largest_column_count - 1)
                        column_names = [i for i in range(0, largest_column_count)]
                        
                        # Read csv
                        df = read_csv(data_file, header=None, delimiter=data_file_delimiter, names=column_names,
                                      index_col=0)
                        count = float(df.loc['Count'][1])
                        dens = round(float(df.loc['Density'][1]),3)
                        
                        df = read_csv(data_file, skiprows=3, header=0, delimiter=data_file_delimiter, 
                                       index_col=0)
                        
                        dur = df['Duration (s)'].loc['Mean']
                        minamp = df['Min. amplitude (uV)'].loc['Mean']
                        maxamp = df['Max. amplitude (uV)'].loc['Mean']
                        ptp = df['Peak-to-peak amplitude (uV)'].loc['Mean']
                        power = df['Power (uV^2)'].loc['Mean']
                        peakfreq = df['Peak power frequency (Hz)'].loc['Mean']
                        
                        dur_sd = df['Duration (s)'].loc['SD']
                        minamp_sd = df['Min. amplitude (uV)'].loc['SD']
                        maxamp_sd = df['Max. amplitude (uV)'].loc['SD']
                        ptp_sd = df['Peak-to-peak amplitude (uV)'].loc['SD']
                        power_sd = df['Power (uV^2)'].loc['SD']
                        peakfreq_sd = df['Peak power frequency (Hz)'].loc['SD']
                        
                        
                        dataset[p,(v*(len(chan)*14)+(c*13))] = count
                        dataset[p,(v*(len(chan)*14)+(c*13)+1)] = dens
                        dataset[p,(v*(len(chan)*14)+(c*13)+2)] = dur
                        dataset[p,(v*(len(chan)*14)+(c*13)+3)] = dur_sd
                        dataset[p,(v*(len(chan)*14)+(c*13)+4)] = minamp
                        dataset[p,(v*(len(chan)*14)+(c*13)+5)] = minamp_sd
                        dataset[p,(v*(len(chan)*14)+(c*13)+6)] = maxamp
                        dataset[p,(v*(len(chan)*14)+(c*13)+7)] = maxamp_sd
                        dataset[p,(v*(len(chan)*14)+(c*13)+8)] = ptp
                        dataset[p,(v*(len(chan)*14)+(c*13)+9)] = ptp_sd
                        dataset[p,(v*(len(chan)*14)+(c*13)+10)] = power
                        dataset[p,(v*(len(chan)*14)+(c*13)+11)] = power_sd
                        dataset[p,(v*(len(chan)*14)+(c*13)+12)] = peakfreq
                        dataset[p,(v*(len(chan)*14)+(c*13)+13)] = peakfreq_sd
                        
                    else:
                        print(f'WARNING: {data_file} not found')
        
        dataset = DataFrame(dataset, index=part, columns=columns)                
        dataset.to_csv(out_dir + f"{evt_name}_{'_'.join(chan)}_fromcsvs.csv")
                