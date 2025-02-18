#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 15 11:33:16 2022

@author: nathancross
"""
from .logs import create_logger
import mne
import numpy as np
from pandas import DataFrame
import copy
from operator import itemgetter
from os import listdir, path, mkdir
from datetime import datetime
from wonambi import Dataset, graphoelement
from wonambi.attr.annotations import Annotations, create_empty_annotations
from scipy.signal import find_peaks, periodogram
import shutil
from sleepecg import detect_heartbeats



def clean_annots(xml_dir, out_dir, rater, keep_evts = None, subs = 'all', 
                 sessions = 'all', logger = create_logger('Clean annotations')):
    
    '''
        Copies the annotations files from inside xml_dir into out_dir and 
        removes all annotations from the files.
        
    '''
    
    logger.debug('')
    logger.debug(f'Time start: {datetime.now()}')
    logger.debug(f'Removing events from files in directory {xml_dir}')  
    logger.debug(f'Keeping events = {keep_evts}')
    logger.debug(f'Saving new annotations files to directory {out_dir}')
    logger.debug('')
    
    # Loop through subjects
    if isinstance(subs, list):
        None
    elif subs == 'all':
            subs = listdir(xml_dir)
            subs = [x for x in subs if not '.' in x]
    else:
        logger.error("'subs' must either be an array of subject ids or = 'all' **CHECK**")
    
    subs.sort()                
    for i, sub in enumerate(subs):
        # Loop through visits
        if sessions == 'all':
            sessions = listdir(xml_dir + '/' + '/' + sub)
            sessions = [x for x in sessions if not '.' in x]
        sessions.sort()
        for v, ses in enumerate(sessions): 
            
            # Define files
            xdir = xml_dir + '/' + sub + '/' + ses + '/'
            xml_files = [x for x in listdir(xdir) if x.endswith('.xml') if not x.startswith('.')] 
            
            if len(xml_files) == 0:                
                logger.warning(f'No xml files found for Subject {sub}, skipping..')
  
            else:
                for xml_file in xml_files:
                    ## Copy annotations file before beginning
                    if not path.exists(out_dir):
                        mkdir(out_dir)
                    if not path.exists(out_dir + sub):
                        mkdir(out_dir + sub)
                    if not path.exists(out_dir + sub + '/' + ses):
                        mkdir(out_dir + sub + '/' + ses)
                    backup = out_dir + sub + '/' + ses + '/'
                    backup_file = (f'{backup}{xml_file[0]}')
                    shutil.copy(xdir + xml_file, backup_file)
                                
                    # Import Annotations file
                    annot = Annotations(backup_file, rater_name=rater)
                    evts = list(set([x['name'] for x in annot.get_events()]))
                    
                    logger.debug(f'Removing events for {sub}, {ses}..')
                    if keep_evts == None:
                        keep_evts = []
                    for ev in evts:
                        if ev not in keep_evts:
                            annot.remove_event_type(name=ev)
                

def remove_event(annot, evt_name, chan = None, stage = None):
    
    ''' Remove event from annotation
        Workaround function because wonambi.attr.annotations.remove_event()
        is not working properly (and I couldn't figure out why).
    '''

    if not chan and not stage:
        evts = [x for x in annot.get_events(name=evt_name)]
    elif not stage:
        evts = [x for x in annot.get_events(name=evt_name) 
                if chan not in x['chan']]
    elif not chan:
        evts = [x for x in annot.get_events(name=evt_name) 
                if stage not in x['stage']]
    
    annot.remove_event_type(name=evt_name)
    grapho = graphoelement.Graphoelement()
    grapho.events = evts       
    grapho.to_annot(annot, evt_name)
    
    return
    

def remove_duplicate_evts(annot, evt_name, chan, stage = None):
    
    ''' Removes any events in annotations that are duplicate.
        Workaround function because wonambi.attr.annotations.remove_event()
        is not working properly (and I couldn't figure out why).
    '''
    evts = annot.get_events(name=evt_name, chan = chan, stage = stage)
    evts_trim = copy.deepcopy(evts)
    for e, event in enumerate(evts[:-1]):
        starttime = event['start']
        i = 0
        for ee, eevent in enumerate(evts_trim):
                if eevent['start'] == starttime:
                    if i == 0:
                        None
                    elif i >0:
                        del(evts_trim[ee])
                    i=+1
    
    evts = [x for x in annot.get_events(name=evt_name) if chan not in x['chan']]
    
    annot.remove_event_type(name=evt_name)
    grapho = graphoelement.Graphoelement()
    grapho.events = evts + evts_trim          
    grapho.to_annot(annot, evt_name)
    
    return


def remove_duplicate_evts_bids(in_dir, out_dir, chan, grp_name, rater, 
                               cat=(0, 0, 0, 0), stage = None, 
                               evt_name = None, part = 'all', visit = 'all', 
                               param_keys = None):

    if path.exists(in_dir):
            print(in_dir + " already exists")
    else:
        mkdir(in_dir)
    
    print('')
    print(f'Time start: {datetime.now()}')
    print(f'Removing duplicates from files in directory {in_dir}')  
    print(f'Event = {evt_name[0]}')
    print('')
    
    if evt_name == None:
        evt_name = ['spindle']
    
    # Loop through subjects
    if isinstance(part, list):
        None
    elif part == 'all':
            part = listdir(in_dir)
            part = [ p for p in part if not '.' in p]
    else:
        print("ERROR: 'part' must either be an array of subject ids or = 'all' **CHECK**")
    
    
    part.sort()                
    for i, p in enumerate(part):
        # Loop through visits
        if visit == 'all':
            visit = listdir(in_dir + '/' + p)
            visit = [x for x in visit if not '.' in x]
        
        print(f'Removing duplicate events for Subject {p}..')
        
        visit.sort()
        for v, vis in enumerate(visit): 
            if not path.exists(out_dir + '/' + p + '/' + vis + '/'):
                mkdir((out_dir + '/' + p + '/' + vis + '/'))
            backup = out_dir + '/' + p + '/' + vis + '/' #backup folder

            # Define files
            rec_dir = in_dir + '/' + p + '/' + vis + '/'
            xml_file = [x for x in listdir(rec_dir) if x.endswith('.xml') if not x.startswith('.')] 
            
            if len(xml_file) == 0:                
                print(f'WARNING: whale_it has not been run for Subject {p}, skipping..')
                
            else:
                # backup file
                backup_file = (backup + xml_file[0])
                shutil.copy(rec_dir + xml_file[0], backup_file) 
                
                # Import Annotations file
                annot = Annotations(backup_file, 
                                    rater_name=rater)
                
                # Run through channels
                for ch, channel in enumerate(chan):
                    chan_ful = channel + ' (' + grp_name + ')'
                                    
                    ### WHOLE NIGHT ###
                    # Select and read data
                    print('Reading data for ' + p + ', visit ' + vis + ' ' + channel)
                    
                    # Retrieve events
                    remove_duplicate_evts(annot, evt_name, chan_ful, stage)

    return


def merge_xmls(in_dirs, out_dir, chan, grp_name, stage = None, evt_name = None, 
               part = 'all', visit = 'all'):

    for pth in in_dirs:
        if not path.exists(pth):
            print(pth + " doesn't exist. Check in_dirs is set correctly.")
            
    else: 
        main_dir = in_dirs[0]
        sec_dir = in_dirs[1]
        
        print('')
        print(f'Time start: {datetime.now()}')
        print(f'Initiating the merging of xml files in directory {main_dir}...') 
        print(f'with the xml files in directory {sec_dir}...') 
        print('')
        
        # Get list of subjects
        if isinstance(part, list):
            None
        elif part == 'all':
                part = listdir(main_dir)
                part = [ p for p in part if not '.' in p]
        else:
            print("ERROR: 'part' must either be an array of subject ids or = 'all' **CHECK**")
        
        
        # Loop through subjects  
        part.sort()              
        for i, p in enumerate(part):
            print(f'Merging xml files for Subject {p}..')
            # Get visits
            if visit == 'all':
                visit = listdir(main_dir + '/' + p)
                visit = [x for x in visit if not '.' in x]
            
            # Loop through visits
            visit.sort()
            for v, vis in enumerate(visit): 
                if not path.exists(main_dir + '/' + p + '/' + vis):
                    print(main_dir + '/' + p + '/' + vis + " doesn't exist. Check data is in BIDS format.")
                
                # Check for output directory and create
                if not path.exists(out_dir + '/'):
                     mkdir(out_dir + '/')
                if not path.exists(out_dir + '/' + p):
                    mkdir(out_dir + '/' + p)
                if not path.exists(out_dir + '/' + p + '/' + vis):
                   mkdir(out_dir + '/' + p + '/' + vis)
                output = out_dir + '/' + p + '/' + vis + '/' #backup folder
    
    
                # Define files
                if main_dir == sec_dir:
                
                    rec_dir = main_dir + '/' + p + '/' + vis + '/'
                    xml_file = [x for x in listdir(rec_dir) if x.endswith('.xml') 
                                if not x.startswith('.')] 
                    xml_file.sort()
                    
                else:
                    rec_dir = main_dir + '/' + p + '/' + vis + '/'
                    xml_file1 = [x for x in listdir(rec_dir) if x.endswith('.xml') 
                                if not x.startswith('.')]
                    rec_dir2 = sec_dir + '/' + p + '/' + vis + '/'
                    xml_file2 = [x for x in listdir(rec_dir2) if x.endswith('.xml') 
                                if not x.startswith('.')]
                    xml_file = xml_file1 + xml_file2
                    
                    
                # Merge files
                if len(xml_file) < 2:                
                    print(f'WARNING: Only 1 xml file found for Subject {p}, skipping..')
                else:
                    # Copy first file to output directory
                    output_file = (output + p + '-merged.xml')
                    shutil.copy(rec_dir + xml_file[0], output_file) 
                    
                    # Import Annotations file
                    annot = Annotations(output_file, 
                                        rater_name=None)
                    
                    for x,xml in enumerate(xml_file[1:]):
                        annot2 = Annotations(rec_dir + xml, 
                                        rater_name=None)
                    
                        # Run through channels
                        for ch, channel in enumerate(chan):
                            chan_ful = [channel + ' (' + grp_name + ')']
                                            
                            print('Merging on channel: {channel}')
                            
                            # Calculate event density (whole night)
                            if evt_name is not None:
                                for e,evt in enumerate(evt_name):
                                    evts = annot2.get_events(name=evt, chan = chan_ful, 
                                                             stage = stage)
                                    annot.add_events(evts)
                            else:
                                evts = annot2.get_events(name=None, chan = chan_ful, 
                                                         stage = stage)
                                annot.add_events(evts)

    return


def rainbow_merge_evts(xml_dir, out_dir, chan, grp_name, rater, segments = None, 
                       events = None, evt_name = None, part = 'all', 
                       visit = 'all'):
    
    # Make output directory
    if not path.exists(out_dir):
        mkdir(out_dir)
    
    if evt_name == None:
        evt_name = ['spindle']
    
    # Loop through subjects
    if isinstance(part, list):
        None
    elif part == 'all':
            part = listdir(xml_dir)
            part = [ p for p in part if not(p.startswith('.'))]
    else:
        print("ERROR: 'part' must either be an array of subject ids or = 'all' **CHECK**")
    
    
    part.sort()                
    for i, p in enumerate(part):
        # Loop through visits
        if visit == 'all':
            visit = listdir(xml_dir + '/' + p)
            visit = [x for x in visit if not(x.startswith('.'))]
        
        if not path.exists(out_dir + '/' + p):
            mkdir(out_dir + '/' + p)
        print(" ")      
        print(f'Merging rainbows for Subject {p}..')
        print(r""" 
                 .##@@&&&@@##.
              ,##@&::%&&%%::&@##.
             #@&:%%000000000%%:&@#
           #@&:%00'         '00%:&@#
          #@&:%0'             '0%:&@#
         #@&:%0                 0%:&@#
        #@&:%0                   0%:&@#
        #@&:%0                   0%:&@#
                """,flush=True)   
        
        
        visit.sort()
        for v, vis in enumerate(visit): 
            if not path.exists(out_dir + '/' + p + '/' + vis + r'/'):
                mkdir(out_dir + '/' + p + '/' + vis + r'/')
            backup = out_dir + '/' + p + '/' + vis + r'/' #backup folder

            # Define files
            xdir = xml_dir + '/' + p + '/' + vis + '/'
            xml_file = [x for x in listdir(xdir) if x.endswith('.xml') if not x.startswith('.')] 
            
            if len(xml_file) == 0:                
                print(f'WARNING: no xml files found for Subject {p}, skipping..')
            elif len(xml_file) > 1:
                print(f'WARNING: multiple xml files found for Subject {p}, visit {v}. Check this! Skipping..')
            else:
                # backup file
                backup_file = (backup + p + '_trunc.xml')
                shutil.copy(xdir + xml_file[0], backup_file) 
                
                # Import Annotations file
                annot = Annotations(backup_file,rater_name=rater)
                if segments is not None:
                    times = []
                    for seg in segments:
                        
                        starts = [x['end'] for x in sorted(annot.get_events(name=seg[0]), 
                                                           key=itemgetter('start')) ]
                        ends = [x['start'] for x in sorted(annot.get_events(name=seg[1]), 
                                                           key=itemgetter('start')) ]
                        for i in range(0,len(starts)):
                            times.append((starts[i],ends[i]))
                else:
                    times = [(0,annot.last_second)]
                
                evts = []
                for t,time in enumerate(times):
                    
                    if events is not None:
                        for e in events:
                            evts.extend(annot.get_events(name=e,time=time)[:])
                    else:
                        evts = annot.get_events(time=time)[:]

            
                names = set([e['name'] for e in evts])  
                newevents = copy.deepcopy(evts)
                d1 = {'name': evt_name[0]}
                [x.update(d1) for x in newevents]
                

                for t,typ in enumerate(names):    
                    annot.remove_event_type(name=typ)
                
                grapho = graphoelement.Graphoelement()
                grapho.events = newevents          
                grapho.to_annot(annot)


    return


def rename_evts(xml_dir, out_dir, part, visit, evts):
    
    '''
       Renames events inside annotations files. The xmls are first copied to out_dir
       and the events are deleted from the copied xmls.
       Inputs are :
           xml_dir, out_dir = input xml_directory and output directory for new xmls
           evts             = dictionary for how events are to be renamed:
                               e.g. evts = {'Moelle2011':'spindle',
                                            'Lacourse2018':'spindle,
                                            'Staresina2015':'so'}
    '''
    if not path.exists(out_dir):
        mkdir(out_dir)
    
    # Get list of subjects
    if isinstance(part, list):
        None
    elif part == 'all':
            part = listdir(xml_dir)
            part = [ p for p in part if not '.' in p]
    else:
        print("ERROR: 'part' must either be an array of subject ids or = 'all' **CHECK**")
    
    
    # Loop through subjects  
    part.sort()              
    for i, p in enumerate(part):
        
        # Get visits
        if visit == 'all':
            visit = listdir(xml_dir + '/' + p)
            visit = [x for x in visit if not '.' in x]
        
        # Loop through visits
        visit.sort()
        for v, vis in enumerate(visit): 
            
            print(f'Subject {p}, visit {vis}..')
            
            # Check for input data
            if not path.exists(xml_dir + '/' + p + '/' + vis):
                print(xml_dir + '/' + p + '/' + vis + " doesn't exist. Check data is in BIDS format.")
            
            # Check for output directory and create
            if not path.exists(out_dir + '/'):
                 mkdir(out_dir + '/')
            if not path.exists(out_dir + '/' + p):
                mkdir(out_dir + '/' + p)
            if not path.exists(out_dir + '/' + p + '/' + vis):
               mkdir(out_dir + '/' + p + '/' + vis)
            output = out_dir + '/' + p + '/' + vis + '/' #backup folder

            # Define files
            xdir = xml_dir + '/' + p + '/' + vis + '/'
            xml_file = [x for x in listdir(xdir) if x.endswith('.xml') if not x.startswith('.')] 
            xml_file.sort()
            
            if len(xml_file) > 1:                
                print(f'WARNING: Multiple xml files for Subject {p}, visit {vis} - check this - skipping..')
                
            else:
                
                # Copy first file to output directory
                output_file = (output + xml_file[0])
                shutil.copy(xdir + xml_file[0], output_file) 
                
                # Import Annotations file
                annot = Annotations(output_file, rater_name=None)
                
                # Rename events
                for e,event in enumerate(evts):
                    old = event
                    new = evts[event]
                    
                    print(f"Renaming event '{old}' to '{new}'")
                    
                    annot.rename_event_type(old, new)
            
                
    return
                

def laplacian_mne(data, oREF, channel, ref_chan, laplacian_rename=False, 
                  renames=None, montage='standard_alphabetic'):
    
    
    ch_names = list(data.axis['chan'][0])
    dig = mne.channels.make_standard_montage(montage)
    
    if oREF:
        dig.rename_channels({oREF:'_REF'}, allow_duplicates=False)
    
    if laplacian_rename:
        dig.rename_channels(renames, allow_duplicates=False)
    
    info = mne.create_info(ch_names, data.s_freq,verbose=None)
    mneobj = mne.io.RawArray(data.data[0],info,verbose=40)
    dic = [{x:'eeg' for x in mneobj.ch_names}]
    mneobj.set_channel_types(dic[0],verbose=40)
    mneobj.info.set_montage(dig,verbose=40)
    a = mneobj.pick(['eeg']).load_data()
    a.set_eeg_reference(ref_chan,verbose=40)
    raw_csd = mne.preprocessing.compute_current_source_density(a, verbose=40)
    data = raw_csd.get_data(picks=channel)
    
    return data



def notch_mne(data, oREF, channel, freq, rename=False,
                 renames=None, montage='standard_alphabetic'):
    
    ch_names = list(data.axis['chan'][0])
    dig = mne.channels.make_standard_montage(montage)
    
    if rename:
        dig.rename_channels(renames, allow_duplicates=False)
        
    if oREF:
        dig.rename_channels({oREF:'_REF'}, allow_duplicates=False)
    
    info = mne.create_info(ch_names, data.s_freq,verbose=40)
    mneobj = mne.io.RawArray(data.data[0],info,verbose=40)
    dic = [{x:'eeg' for x in mneobj.ch_names}]
    mneobj.set_channel_types(dic[0],verbose=40)
    mneobj.info.set_montage(dig,verbose=40)
    a = mneobj.pick(['eeg']).load_data()
    anotch = a.notch_filter(freq,verbose=40)
    data = anotch.get_data(picks=channel)
    
    return data

def notch_mne2(data, oREF, channel, rename=False, renames=None,
               montage='standard_alphabetic'):
    
    ch_names = list(data.axis['chan'][0])
    dig = mne.channels.make_standard_montage(montage)
    
    if rename:
        dig.rename_channels(renames, allow_duplicates=False)
        
    if oREF:
        dig.rename_channels({oREF:'_REF'}, allow_duplicates=False)
        
    info = mne.create_info(ch_names, data.s_freq,verbose=40)
    mneobj = mne.io.RawArray(data.data[0],info,verbose=40)
    dic = [{x:'eeg' for x in mneobj.ch_names}]
    mneobj.set_channel_types(dic[0],verbose=40)
    mneobj.info.set_montage(dig,verbose=40)
    a = mneobj.pick(['eeg']).load_data()
    anotch = a.notch_filter(freqs=None, filter_length='auto', 
                            notch_widths=None, trans_bandwidth=1, 
                            method='spectrum_fit', verbose=None)
    data = anotch.get_data(picks=channel)
    
    return data
    

def bandpass_mne(data, oREF, channel, highpass, lowpass, rename=False,
                 renames=None, montage='standard_alphabetic'):
    
    ch_names = list(data.axis['chan'][0])
    dig = mne.channels.make_standard_montage(montage)
    
    if rename:
        dig.rename_channels(renames, allow_duplicates=False)
        
    if oREF:
        dig.rename_channels({oREF:'_REF'}, allow_duplicates=False)
        
    info = mne.create_info(ch_names, data.s_freq,verbose=40)
    mneobj = mne.io.RawArray(data.data[0],info,verbose=40)
    dic = [{x:'eeg' for x in mneobj.ch_names}]
    mneobj.set_channel_types(dic[0],verbose=40)
    mneobj.info.set_montage(dig,verbose=40)
    a = mneobj.pick(['eeg']).load_data()
    a_filt = a.filter(highpass, lowpass,verbose=40)
    data = a_filt.get_data(picks=channel)
    
    return data

def csv_stage_import(edf_file, xml_file, hypno_file, rater):
    
    ''' This function creates a new annoations file and imports staging from a csv 
        file. The staging should be in the 'Alice' format, for further information 
        see: wonambi.attr.annotations.
    '''

    dataset = Dataset(edf_file)
    create_empty_annotations(xml_file, dataset)
    annot = Annotations(xml_file)
    annot.import_staging(hypno_file, 'alice', rater_name=rater, rec_start=dataset.header['start_time'])
    

def infer_eeg(dset, logger = create_logger('Infer EEG')):
    
    ''' This function reads a recording file (e.g. edf) and searches for likely
        EEG channel names. IF the EEG channels are not easly identifiable from
        the edf, then nothing is returned. IF there are multiple (>2) options,
        then the channels will be QC'd and the channels with the best quality 
        will be returned. #NOTE: Central channels will be given highest priority.
    '''
    
    h = dset.header
    
    # 1. Look for central channels first
    eeg = [x for x in h['chan_name'] if any(l in x for l in ['C1', 'C2', 'C3',
                                                             'C4', 'C5', 'C6',
                                                             'Cz', 'CZ'])]  
    
    # 2. If no central channels, look for frontal:
    if len(eeg) == 0:
        eeg = [x for x in h['chan_name'] if any(l in x for l in ['F1', 'F2', 'F3',
                                                                 'F4', 'F5', 'F6',
                                                                 'Fz', 'FZ'])]          
    # 3. Look for any channels names beginning with E
    if len(eeg) == 0:
        eeg = [x for x in h['chan_name'] if 'E' in x if not 
                                                any(l in x for l in ['EMG', 
                                                                     'ECG',
                                                                     'EOG',
                                                                     'Chin',
                                                                     ])]
               

    # If still nothing or when all electrode names contain 'E' (e.g. EGI nets) 
    if len(eeg) == 0 or len(eeg) > 10:
        logger.error('Unable to determine EEG channels from recording. EEG names must be explicitly specified.')
        logger.info('Check documentation for how to specify channels:')
        logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
        logger.info('-' * 10)
        return None
    else:
        eeg = choose_best_eeg(eeg, dset, 1, logger)
            
    return eeg 
    
    
def infer_eog(dset, logger = create_logger('Infer EOG')):
    
    ''' This function reads a recording file (e.g. edf) and searches for likely
        EOG channel names. IF the EOG channels are not easly identifiable from
        the edf, then nothing is returned. IF there are multiple (>2) options,
        then the channels will be QC'd and the channels with the best quality 
        will be returned.
    '''
    
    h = dset.header
    
    # First and most obvious EOG name
    eyes = [x for x in h['chan_name'] if 'EOG' in x]  
    
    # If EOG are not labelled as 'EOG'
    if len(eyes) == 0:
        eyes = [x for x in h['chan_name'] if 'E' in x if not 'EMG' in x 
                                                      if not 'ECG' in x]    
    
    # If still nothing or when all electrode names contain 'E' (e.g. EGI nets) 
    if len(eyes) == 0 or len(eyes) > 6:
        logger.error('Unable to determine EOG channels from recording. EOG names must be explicitly specified.')
        logger.info('Check documentation for how to specify channels:')
        logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
        logger.info('-' * 10)
        return None

    elif len(eyes) > 2:
        
        loc = [x for x in eyes if any(l in x for l in ['l', 'L'])]
        loc = choose_best_eog(loc, dset, 1, logger)
        
        roc = [x for x in eyes if any(l in x for l in ['r', 'R'])]
        roc = choose_best_eog(roc, dset, 1, logger)
        
        eyes = roc + loc
            
    return eyes    

def infer_emg(dset, logger = create_logger('Infer EMG')):
    
    ''' This function reads a recording file (e.g. edf) and searches for likely
        EOG channel names. IF the EOG channels are not easly identifiable from
        the edf, then nothing is returned. IF there are multiple (>2) options,
        then the channels will be QC'd and the channels with the best quality 
        will be returned.
    '''
    
    h = dset.header
    
    # First and most obvious EMG name
    emg = [x for x in h['chan_name'] if any(l in x for l in ['EMG', 
                                                             'Chin', 
                                                             'CHIN'])]  
     
    
    # If still nothing or when all electrode names contain 'E' (e.g. EGI nets) 
    if len(emg) == 0 or len(emg) > 6:
        logger.error('Unable to determine EMG channels from recording. EOG names must be explicitly specified.')
        logger.info('Check documentation for how to specify channels:')
        logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
        logger.info('-' * 10)
        return None

    elif len(emg) > 1:
        emg = choose_best_emg(emg, dset, 1, logger)
            
    return emg  
        
### Quality Control Checks    
def gini(x, w = None):
    
    ''' Calculates the Gini coefficient:
                    a measure of statistical dispersion calculated by comparing 
                    the actual Lorenz curve to the diagonal line of equality.
    '''

    x = np.asarray(x)
    if w is not None:
        w = np.asarray(w)
        sorted_indices = np.argsort(x)
        sorted_x = x[sorted_indices]
        sorted_w = w[sorted_indices]
        # Force float dtype to avoid overflows
        cumw = np.cumsum(sorted_w, dtype=float)
        cumxw = np.cumsum(sorted_x * sorted_w, dtype=float)
        return (np.sum(cumxw[1:] * cumw[:-1] - cumxw[:-1] * cumw[1:]) / 
                (cumxw[-1] * cumw[-1]))
    else:
        sorted_x = np.sort(x)
        n = len(x)
        cumx = np.cumsum(sorted_x, dtype=float)
        # The above formula, with all weights equal to 1 simplifies to:
        return (n + 1 - 2 * np.sum(cumx) / cumx[-1]) / n

def moving_average_time(a, s_freq, win=5):
    n = win*s_freq*60
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

def moving_average(a, win=100):
    n = win
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n


def choose_best_eeg(chans, dset, num = 2, logger = create_logger('QC EEG')):
    
    quals = DataFrame(data=None, dtype = float, 
                          columns=['integral',
                                  'stdev',
                                  'time_not_in_low_range',
                                  'time_not_flatline',
                                  '99th_quantile',
                                  'gini_coeff'])
    
    ch_dat = dset.read_data(chan = chans)
    for i, name in enumerate(chans):
        logger.debug(f'Checking quality of {name}')
        s_freq = dset.header['orig']['n_samples_per_record'][chans.index(name)]
        e = moving_average_time(abs(ch_dat.data[0][i]),s_freq,5)
        quals.loc[i,'integral'] = sum(e) # (AUC)
        quals.loc[i,'stdev'] = np.std(e)
        quals.loc[i,'time_not_in_low_range'] = 1 - sum(np.logical_and(e>=0, e<=0.5))/len(e)
        quals.loc[i,'time_not_flatline'] = 1 - (sum(np.diff(e)==0)/len(np.diff(e)))
        quals.loc[i,'99th_quantile'] = np.quantile(e,0.99)
        quals.loc[i,'gini_coeff'] = gini(e)
    
    winners = DataFrame([quals[x].rank() for x in quals.columns]).T
    idx = winners.sum(axis = 1).nlargest(num).index.to_list()
    eeg_names = [chans[x] for x in idx]
    logger.debug(f'Best EEG channels based on auto-QC are: {eeg_names}')
    
    return eeg_names



def choose_best_eog(chans, dset, num = 2, logger = create_logger('QC EOG')):
    
    quals = DataFrame(data=None, dtype = float, 
                          columns=['integral',
                                  'stdev',
                                  'time_not_in_low_range',
                                  'time_not_flatline',
                                  '99th_quantile',
                                  'gini_coeff'])
    
    ch_dat = dset.read_data(chan = chans)
    for i, name in enumerate(chans):
        logger.debug(f'Checking quality of {name}')
        s_freq = dset.header['orig']['n_samples_per_record'][chans.index(name)]
        e = moving_average_time(abs(ch_dat.data[0][i]),s_freq,5)
        quals.loc[i,'integral'] = sum(e) # (AUC)
        quals.loc[i,'stdev'] = np.std(e)
        quals.loc[i,'time_not_in_low_range'] = 1 - sum(np.logical_and(e>=0, e<=0.5))/len(e)
        quals.loc[i,'time_not_flatline'] = 1 - (sum(np.diff(e)==0)/len(np.diff(e)))
        quals.loc[i,'99th_quantile'] = np.quantile(e,0.99)
        quals.loc[i,'gini_coeff'] = gini(e)
    
    winners = DataFrame([quals[x].rank() for x in quals.columns]).T
    idx = winners.sum(axis = 1).nlargest(num).index.to_list()
    eog_names = [chans[x] for x in idx]
    
    
    
    logger.debug(f'Best EOG channels based on auto-QC are: {eog_names}')
    
    return eog_names




def choose_best_ecg(chans, dset, num = 1, logger = create_logger('QC ECG')):
    
    quals = DataFrame(data=None, columns=['hr_qual',
                                          'stdev',
                                          'time_not_flatline',
                                          'PSD_qual'])
    ecg_dat = dset.read_data(chan=chans)
    
    for i, name in enumerate(chans):
        logger.debug(f'Checking quality of {name}')
        s_freq = dset.header['orig']['n_samples_per_record'][chans.index(name)]
        e = ecg_dat.data[0][i]
        
        # Check BPM
        beats = detect_heartbeats(e, s_freq)
        hr = len(beats)/(len(e)/(s_freq*60))
        if np.logical_and(hr>40, hr<120):
            logger.debug(f'Heart rate ({round(hr)}bpm) seems physiologically plausible. Continuing...')
            quals.loc[i,'hr_qual'] = 2
        elif np.logical_and(hr>120, hr<200):
            quals.loc[i,'hr_qual'] = 1
            logger.warning(f'Heart rate ({round(hr)}bpm) is physiologically plausible but high. Check...')
        else:
            quals.loc[i,'hr_qual'] = 0
            logger.warning(f'Heart rate ({round(hr)}bpm) is NOT physiologically plausible. Signal quality is likely poor...')
        
        # Check Std Dev.
        quals.loc[i,'stdev'] = np.std(abs(e)) #standard deviation 
        
        # Check for flatlines
        quals.loc[i,'time_not_flatline'] = 1 - (sum(np.diff(e)==0)/len(np.diff(e)))
        
        # Check PSD        
        f_p, Pxx_spec = periodogram(e, s_freq)
        Pxx_spec_sm = np.concatenate((np.zeros(2000),
                                      moving_average(Pxx_spec, 4000),
                                      np.zeros(1999)))
        top = int((len(e)/s_freq)*20) # corresponds to 20Hz
        peaks = find_peaks(Pxx_spec_sm[:top], prominence=800)[0]
        
        if len(peaks) > 4:
            quals.loc[i,'PSD_qual'] = 2
        elif len(peaks) > 2:
            quals.loc[i,'PSD_qual'] = 1
        else:
            print(f'WARNING: PSD of {name} looks unusual..')
            quals.loc[i,'PSD_qual'] = 0

    winners = DataFrame([quals[x].rank() for x in quals.columns]).T
    idx = winners.sum(axis = 1).nlargest(num).index.to_list()
    ecg_name = [chans[x] for x in idx]
    logger.debug(f'Best ECG channel based on auto-QC is: {ecg_name}')

    return ecg_name


def choose_best_emg(chans, dset, num = 1, logger = create_logger('QC EMG')):

    emg_dat = dset.read_data(chan=chans)
    quals = DataFrame(data=None, dtype = float, 
                      columns=['integral',
                               'stdev',
                               'time_not_in_low_range',
                               'time_not_flatline',
                               '99th_quantile',
                               'gini_coeff'])
    
    for i, name in enumerate(chans):
        logger.debug(f'Checking quality of {name}')
        s_freq = dset.header['orig']['n_samples_per_record'][chans.index(name)]
        e = moving_average_time(abs(emg_dat.data[0][i]),s_freq,5)
        quals.loc[i,'integral'] = sum(e) #intragral (AUC)
        quals.loc[i,'stdev'] = np.std(e) #standard deviation 
        quals.loc[i,'time_not_in_low_range'] = 1 - sum(np.logical_and(e>=0, e<=0.5))/len(e)
        quals.loc[i,'time_not_flatline'] = 1 - (sum(np.diff(e)==0)/len(np.diff(e)))
        quals.loc[i,'99th_quantile'] = np.quantile(e,0.99)
        quals.loc[i,'gini_coeff'] = gini(e)
    
    winners = DataFrame([quals[x].rank() for x in quals.columns]).T
    idx = winners.sum(axis = 1).nlargest(num).index.to_list()
    emg_name = [chans[x] for x in idx]
    print(f'Best EMG channel based on auto-QC is: {emg_name}')
    
    return emg_name



    