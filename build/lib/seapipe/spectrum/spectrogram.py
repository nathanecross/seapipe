#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Mon Mar 27 14:18:54 2023

@author: nathancross
"""

" Spectrogram on events."

from copy import deepcopy
from datetime import datetime, date
from numpy import (append, flip, mean, where, zeros)
from os import listdir, mkdir, path, walk
from pandas import DataFrame
from safepickle import dump, load
from scipy.signal import spectrogram
from scipy.stats import binned_statistic
import sys
from termcolor import colored
from ..utils.misc import laplacian_mne, notch_mne, notch_mne2
from ..utils.logs import create_logger, create_logger_outfile
from ..utils.load import load_channels, rename_channels

from wonambi.attr import Annotations 
from wonambi.trans import fetch
from wonambi import Dataset
import shutil

def event_spectrogram(self, in_dir, xml_dir, out_dir, subs, sessions, stage, 
                      cycle_idx, chan, ref_chan, rater, grp_name, evt_type, buffer, 
                      invert, cat, filter_opts, outfile=False, filetype = '.edf', 
                      progress=True, tracking = None):

        
    ### 0.a Set up logging
    flag = 0
    tracking = self.tracking
    if outfile == True:
        today = date.today().strftime("%Y%m%d")
        now = datetime.now().strftime("%H:%M:%S")
        logfile = f'{self.log_dir}/detect_spindles_Spectrogram_{today}_log.txt'
        logger = create_logger_outfile(logfile=logfile, name='Detect spindles')
        logger.info('')
        logger.info(f"-------------- New call of 'Spectrogram' evoked at {now} --------------")
    elif outfile:
        logfile = f'{self.log_dir}/{outfile}'
        logger = create_logger_outfile(logfile=logfile, name='Detect spindles')
    else:
        logger = create_logger('Detect spindles')
    
    logger.debug(r"""    Computing spectrogram...
          
                 
             |................................................
             |........................ ..  ...................
             |...............'.......    .....................
             |....  ......... .......  .......'...............
             |....   ..'..... .......  '......'......'........
             |...       ..... .......   ......'......'........
             |...    .. ..... ................'......'........
             |...  .... ..... ............... '......'........
             |...  ...; .... :......'........  ...... ........
             |... ....   •.,    •... '.......  '..... ........
             |.......;   '.;  . .... '......    ..... ........
             |......'     '.  . '...   ...       ..,  '.......
             |________________________________________________
                     (Hz)
          
          """)
    
    # Make output directory
    if not path.exists(out_dir):
        mkdir(out_dir)
    
    # b. Check input list
    subs = self.subs
    if isinstance(subs, list):
        None
    elif subs == 'all':
            subs = listdir(self.rec_dir)
            subs = [p for p in subs if not '.' in p]
    else:
        logger.error("'subs' must either be an array of subject ids or = 'all' ") 
    
    
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
            tracking[f'{sub}'][f'{ses}'] = {'spectrogram':{}} 

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
                annot = Annotations(backup_file, rater_name=rater)
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
                
                # 5 Channel stuff
                # a. Rename channel for output file (if required)
                if newchans:
                    fnamechan = newchans[ch]
                else:
                    fnamechan = ch
                # b. Set ref channel name for log output
                if not chanset[ch]:
                    logchan = ['(no re-refrencing)']
                else:
                    logchan = chanset[ch]
                    
                if filter_opts['laplacian']:
                    chans = filter_opts['lapchan']
                else:
                    chans = [ch]
                    ref_chan=chanset[chans]

                # h. Read data
                logger.debug(f"Reading EEG data for {sub}, {ses}, {str(ch)}:{'-'.join(logchan)}")
                logger.debug(f"Using filter settings:\nNotch filter: {filter_opts['notch']}")
                logger.debug(f"Notch harmonics filter: {filter_opts['notch_harmonics']}")
                logger.debug(f"Laplacian filter: {filter_opts['laplacian']}")

                
                try:
                    segments = fetch(dset, annot, cat=cat, stage=self.stage, 
                                     cycle=cycle, evt_type=evt_type,
                                     buffer=buffer, reject_artf=self.reject)

                    segments.read_data(chans, ref_chan, grp_name=self.grp_name)
                except Exception as error:
                    logger.error(type(error).__name__, "–", error)
                    flag+=1
                    break
                
                if len(segments) <1:
                    logger.error(f'No events found for {sub}, {ses}, {str(ch)}, skipping...')
                    continue
                
                nsegs=[]
                # Check concatenation of stages
                if cat[0] == 0 and cat[1] == 0:
                    logger.debug('Splitting cycles and stages')
                elif cat[0] == 0:
                    logger.debug('Splitting cycles')
                    for cy in cycle:
                        segs = [s for s in segments if cy in s['cycle']]
                        nsegs.append(segs)
                elif cat[1] == 0:
                    logger.debug('Splitting stages')
                    for st in stage:
                        segs = [s for s in segments if st in s['stage']]
                        nsegs.append(segs)
                else:
                    nsegs = [segments]
                




                
    # Check input visits
    if isinstance(sessions, list):
        None
    elif sessions == 'all':
        lenvis = set([len(next(walk(in_dir + x))[1]) for x in subs])
        if len(lenvis) > 1:
            print(colored('WARNING |', 'yellow', attrs=['bold']),
                  colored('number of visits are not the same for all subjects.',
                          'white', attrs=['bold']))
            print('')
    else:
        print('')
        print(colored('ERROR |', 'red', attrs=['bold']),
              colored("'visit' must either be an array of subject ids or = 'visit'" ,
                      'cyan', attrs=['bold']))
        print('')
    
    # Ensure correct concatenation setup
    cat = list(cat)
    cat[2] = 0 #enforce non-concatenation of signal
    cat[3] = 0 #enforce non-concatenation of event types
    cat = tuple(cat)
    
    # Loop through participants and visits
    part.sort()
    for i, p in enumerate(part):
        
        if not path.exists(out_dir + '/' + p):
            mkdir(out_dir + '/' + p)
        
        if visit == 'all':
            visit = listdir(xml_dir + '/' + p)
            visit = [x for x in visit if not '.' in x]
        visit.sort()    
        for j, vis in enumerate(visit): 
            if not path.exists(xml_dir + '/' + p + '/' + vis + '/'):
                print(colored('WARNING |', 'yellow', attrs=['bold']),
                      colored(f'input folder missing for Subject {p}, visit {j} ({vis}), skipping...',
                              'white', attrs=['bold']))
                print('')
                continue
            else:
                
                if not path.exists(out_dir + '/' + p + '/' + vis):
                    mkdir(out_dir + '/' + p + '/' + vis)
                rec_file = [s for s in listdir(rec_dir + '/' + p + '/' + vis) if 
                            ('.edf') in s or ('.rec') in s or ('.eeg')  in s if not s.startswith(".")]
                xml_file = [x for x in listdir(xml_dir + '/' + p + '/' + vis) if 
                            x.endswith('.xml') if not x.startswith(".")] 
                
                # Open recording and annotations files (if existing)
                if len(xml_file) == 0:
                    print(colored('WARNING |', 'yellow', attrs=['bold']),
                          colored(f'annotations does not exist for Subject {p}, visit {j} ({vis}) - check this. Skipping...',
                                  'white', attrs=['bold']))
                    print('')
                elif len(xml_file) >1:
                    print(colored('WARNING |', 'yellow', attrs=['bold']),
                          colored(f'multiple annotations files exist for Subject {p}, visit {j} ({vis}) - check this. Skipping...',
                                  'white', attrs=['bold']))
                    print('')
                else:
                    
                    ###########################            DEBUGGING              ###########################
                    with open(out_dir + f'/debug_{p}.txt', 'w') as f:
                        f.write('opening participant edf and xml')
                    ###########################            DEBUGGING              ###########################
    
                # Load data
                dset = Dataset(rec_dir + '/' + p + '/' + vis + '/' + rec_file[0]) 
                s_freq = dset.header['s_freq']
                annot = Annotations(xml_dir + '/' + p + '/' + vis + '/' + xml_file[0], 
                                    rater_name=None)
                
                # Get sleep cycles
                if cycle_idx is not None and cat[0] == 1:
                    all_cycles = annot.get_cycles()
                    scycle = [all_cycles[i - 1] for i in cycle_idx if i <= len(all_cycles)]
                else:
                    scycle = [None]
                
                # Loop through channels
                for k, ch in enumerate(chan):
                    
                    print(f'Reading data for {p}, visit {vis}, channel {ch}')
                    chan_full = ch + ' (' + grp_name + ')'
                    
                    
                    # Loop through sleep cycles 
                    for l, cyc in enumerate(scycle):
                        print('')
                        # Select and read data
                        if cycle_idx is not None:
                            print(f'Analysing, cycle {l+1}')
                        else:
                            print('Analysing, all cycles')
                        print('')
                        print('Using filter settings:')
                        print('')
                        print(colored('Notch filter:','white', attrs=['bold']),
                              colored(f"{filter_opts['notch']}", 'yellow', attrs=['bold']))
                        print(colored('Notch harmonics filter:','white', attrs=['bold']),
                              colored(f"{filter_opts['notch_harmonics']}", 'yellow', attrs=['bold']))
                        print(colored('Laplacian filter:','white', attrs=['bold']),
                              colored(f"{filter_opts['laplacian']}", 'yellow', attrs=['bold']))
                        
                        # Fetch segments from recording
                        segments = fetch(dset, annot, cat=cat, chan_full=[chan_full], 
                                         cycle=cyc, evt_type=evt_type, stage=stage,
                                         buffer=buffer)
                        if filter_opts['laplacian'] or filter_opts['notch'] or filter_opts['notch_harmonics']:
                            chans = filter_opts['lapchan']
                        else:
                            chans = [ch]
                        segments.read_data(chan=chans, ref_chan=ref_chan)
                        
                        if len(segments) <1:
                            print(colored('WARNING |', 'yellow', attrs=['bold']),
                                  colored('No segments found.',
                                          'white', attrs=['bold']))
                            
                        nsegs=[]
                        # Check concatenation of stages
                        if cat[1] == 0:
                            print('Splitting stages')
                            for s, st in enumerate(stage):
                                segs = [s for s in segments if st in s['stage']]
                                nsegs.append(segs)
                        else:
                            nsegs = [segments]
                        
                        for s in range(len(nsegs)):
                            segments = nsegs[s]
                            
                            if cat[1] == 1:
                                stagename = ''.join(stage) 
                            else: 
                                stagename = stage[s]
                                
                            print('')
                            print(f'Stage {stagename}')
                            print('Creating spectrogram')
                            print(f'No. Segments = {len(segments)}')
                           
                            z=0
                            for m, seg in enumerate(segments):
                                
                                # Print out progress
                                if progress:
                                    z +=1
                                    j = z/len(segments)
                                    sys.stdout.write('\r')
                                    sys.stdout.write(f"Progress: [{'=' * int(50 * j):{50}s}] {int(100 * j)}%")
                                    sys.stdout.flush()
                                
                                # Select data from segment
                                data = seg['data']
                                
                                # Find sampling frequency
                                s_freq = data.s_freq
    
                                # Apply filtering (if necessary)
                                if filter_opts['notch']:
                                    selectchans = list(data.chan[0])
                                    data.data[0] = notch_mne(data, oREF=filter_opts['oREF'], 
                                                                channel=selectchans, 
                                                                freq=filter_opts['notch_freq'],
                                                                rename=filter_opts['chan_rename'],
                                                                renames=filter_opts['renames'])
                                
                                if filter_opts['notch_harmonics']: 
                                    selectchans = list(data.chan[0])
                                    data.data[0] = notch_mne2(data, oREF=filter_opts['oREF'], 
                                                              channel=selectchans,
                                                              rename=filter_opts['chan_rename'],
                                                              renames=filter_opts['renames'])
                                
                                if filter_opts['laplacian']:
                                    data = laplacian_mne(data, oREF=filter_opts['oREF'], channel=ch, 
                                                         ref_chan=ref_chan, 
                                                         laplacian_rename=filter_opts['laplacian_rename'], 
                                                         renames=filter_opts['renames'])
                                    dat = data[0]
                                else:
                                    dat = data()[0][0]
                                
                                
                                # Check polarity of recording
                                if isinstance(polar, list):
                                    polarity = polar[i]
                                else:
                                    polarity = polar
                                if polarity == 'opposite':
                                    dat = dat*-1 
                        
                                if m == 0:
                                    dat = dat[int((len(dat)/2)-(s_freq*buffer)):int((len(dat)/2)+(s_freq*buffer))]
                                    out = dat.reshape((1,len(dat)))
                                else:
                                    dat = dat[int((len(dat)/2)-(s_freq*buffer)):int((len(dat)/2)+(s_freq*buffer))]
                                    dat = dat.reshape((1,len(dat)))
                                    out = append(out,dat,axis=0)
                                     
                        
                        # Calaculate spectrogram
                        f,t,spect = spectrogram(out, fs=s_freq, nfft=filter_opts['nfft'], 
                                                noverlap=filter_opts['noverlap'],)
                        
                        # Apply mask to extract only frequencies you want
                        freqs = filter_opts['bandpass']
                        highpass = where(f==freqs[0])[0][0]
                        lowpass = where(f==freqs[1])[0][0]
                        spect = flip(spect[:,highpass:lowpass,:], axis=1)
                        f = f[highpass:lowpass]
                        
                        # Prepare filename
                        if cat[0] == 1:
                            cyclename = 'wholenight'
                        else: 
                            cyclename = f'cycle{l+1}' 
                            
                        # Save average spectrogram to dataframe
                        spectav = mean(spect,axis=0)
                        d = DataFrame(spectav)
                        d.columns = t
                        d.index = flip(f)
                        d.to_csv(path_or_buf=out_dir + '/' + p + '/' + vis + '/' + 
                                 p + '_' + vis + '_' + ch + '_' + stagename +  '_' + cyclename + 
                                 '_spectrogram.csv', sep=',')
                        
                        # Save spectrogram for every event to pickle file
                        spectro = {'spect':spect, 'f':f, 't':t}
                        with open(out_dir + '/' + p + '/' + vis + '/' + 
                                 p + '_' + vis + '_' + ch + '_' + stagename + '_' + cyclename + 
                                  '_spectrogram.p', 'wb') as ff:
                             dump(spectro, ff)
                             
                        # Save raw and averaged event signal to pickle file
                        datav = mean(out, axis=0)
                        signal = {'raw':out, 'average':datav}
                        with open(out_dir + '/' + p + '/' + vis + '/' + 
                                 p + '_' + vis + '_' + ch + '_' + stagename + '_' + cyclename + 
                                  '_signal.p', 'wb') as ff:
                             dump(signal, ff)
    
    print('The function event_spectrogram completed without error.')
    return
                
                
                
                

def event_spectrogram_grouplevel(in_dir, out_dir, part, visit, chan, stage, cat, 
                                 cycle_idx):                        
             
    '''
    This script combines the output from the function event_spectrogram, and formats it
    in a group-level dataframe for statistical analyses.
    The outputs provided by this script will be, for each visit and EEG channel:
        1. A pickle file (per channel and visit) with the night-average spectrogram for
            every subject and the night-average signal of each event.
        2. A .csv file (per channel and visit) containing the group averaged spectrogram.
        
    '''   
    
    # Make output directory
    if not path.exists(out_dir):
        mkdir(out_dir)
    
    ## BIDS CHECKING
    # Check input participants
    if isinstance(part, list):
        None
    elif part == 'all':
            part = listdir(in_dir)
            part = [ p for p in part if not '.' in p]
    else:
        print('')
        print(colored('ERROR |', 'red', attrs=['bold']),
              colored("'part' must either be an array of subject ids or = 'all'" ,
                      'cyan', attrs=['bold']))
        print('')
        
    # Check input visits
    if isinstance(visit, list):
        None
    elif visit == 'all':
        lenvis = set([len(next(walk(in_dir + x))[1]) for x in part])
        if len(lenvis) > 1:
            print(colored('WARNING |', 'yellow', attrs=['bold']),
                  colored('number of visits are not the same for all subjects.',
                          'white', attrs=['bold']))
            print('')
        visit = list(set([y for x in part for y in listdir(in_dir + x)  if '.' not in y]))
    else:
        print('')
        print(colored('ERROR |', 'red', attrs=['bold']),
              colored("'visit' must either be an array of subject ids or = 'visit'" ,
                      'cyan', attrs=['bold']))
        print('')
            
    # Create output dataframe
    if cycle_idx is not None:
        all_ampbin = zeros((7, len(part), 6), dtype='object')
    else:
        all_ampbin = zeros((7, len(part)), dtype='object')
    
    # Check for stage setup
    if cat[1] == 1:
        stage = [''.join(stage) ]

    for st, stagename in enumerate(stage): # Loop through stages
        for k, ch in enumerate(chan):      # Loop through channels
            print('')
            print(f'CHANNEL {ch}')
            
            
            for j, vis in enumerate(visit): 
                z=0
                index=[]
                part.sort()
                spectgrp = zeros((len(part)), dtype='object')
                signalgrp = zeros((len(part)), dtype='object')
                for i, p in enumerate(part):    # Loop through participants
                    index.append(p)
                    if not path.exists(in_dir + '/' + p + '/' + vis + '/'):
                        print(colored('WARNING |', 'yellow', attrs=['bold']),
                              colored(f'input folder missing for Subject {p}, visit {vis}, skipping..',
                                      'white', attrs=['bold']))
                        continue
                    else:
                        
                        # Spectrogram
                        # Define pickle files for spectrogram
                        p_files = [s for s in listdir(in_dir + '/' + p + '/' + vis) 
                                   if '_spectrogram.p' in s] 
                        p_files = [s for s in p_files if ch in s]
                        p_files = [s for s in p_files if '_'+stagename+'_' in s]
                        
                        # Open files containing spectrogram (if existing)
                        if len(p_files) == 0:
                            print(colored('WARNING |', 'yellow', attrs=['bold']),
                                  colored(f'spectrogram does not exist for Subject {p}, visit {vis}, stage {stagename}, channel {ch}. Skipping..',
                                          'white', attrs=['bold']))
                        elif len(p_files)>1:
                            print(colored('WARNING |', 'yellow', attrs=['bold']),
                                  colored(f'multiple spectrogram files exist for Subject {p}, visit {vis}, stage {stagename}, channel {ch}. Skipping..',
                                          'white', attrs=['bold']))
                            print('')
                        else:
                            print(f'Extracting spectrogram for ... Subject {p}, visit {vis}')
                            pfile = in_dir + '/' + p + '/' + vis + '/' + p_files[0]
                            with open(pfile, 'rb') as ff:
                                pick = load(ff)
                            
                            spect = pick['spect']
                            f = pick['f']
                            t = pick['t']
                            spectgrp[i] = mean(spect, axis=0)
                            
                        ## Raw signal  
                        # Define pickle files for raw signal
                        p_files = [s for s in listdir(in_dir + '/' + p + '/' + vis) 
                                   if '_signal.p' in s] 
                        p_files = [s for s in p_files if ch in s]
                        p_files = [s for s in p_files if '_'+stagename+'_' in s]
                        
                        # Open files containing event signal (if existing)
                        if len(p_files) == 0:
                            print(colored('WARNING |', 'yellow', attrs=['bold']),
                                  colored(f'signal file does not exist for Subject {p}, visit {vis}, stage {stagename}, channel {ch}. Skipping..',
                                          'white', attrs=['bold']))
                        elif len(p_files)>1:
                            print(colored('WARNING |', 'yellow', attrs=['bold']),
                                  colored(f'multiple signal files exist for Subject {p}, visit {vis}, stage {stagename}, channel {ch}. Skipping..',
                                          'white', attrs=['bold']))
                            print('')
                        else:
                            print(f'Extracting raw signal for ... Subject {p}, visit {vis}')
                            pfile = in_dir + '/' + p + '/' + vis + '/' + p_files[0]
                            with open(pfile, 'rb') as ff:
                                pick = load(ff)
                            
                            signal = pick['average']
                            signalgrp[i] = binned_statistic(range(0,len(signal)), signal,  
                                                            bins=spectgrp[i].shape[1])[0]
                            
                    
                # Save all spectrograms and event signal to pickle file
                grouped_spectrogram = {'spectrogram':spectgrp, 'signal':signalgrp}
                with open(out_dir + ch + '_' + stagename + '_visit_' + vis +
                          '_group_spectrogram.p', 'wb') as ff:
                     dump(grouped_spectrogram, ff)
                 
                    
                # Save average spectrogram to dataframe
                spectav = mean(spectgrp,axis=0)
                d = DataFrame(spectav)
                d.columns = t
                d.index = flip(f)
                d.to_csv(path_or_buf=out_dir + ch + '_' + stagename +  '_visit_' + vis + 
                         '_spectrogram.csv', sep=',')  
                            
    return
                
                
                
                


