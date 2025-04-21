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
from wonambi.trans import fetch
from wonambi.detect.spindle import transform_signal
from scipy.signal import butter, filtfilt, find_peaks,  periodogram
import shutil
from .logs import create_logger
from .load import check_adap_bands
from .audit import check_fooof

                
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


def remove_duplicate_evts(annot, evt_name, chan = None, stage = None):
    
    ''' Removes any events in annotations that are duplicate.
        Workaround function because wonambi.attr.annotations.remove_event()
        is not working properly (and I couldn't figure out why).
    '''
    evts = annot.get_events(name = evt_name, chan = chan, stage = stage)
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
    
    evts = [x for x in annot.get_events(name = evt_name) if chan not in x['chan']]
    
    annot.remove_event_type(name = evt_name)
    grapho = graphoelement.Graphoelement()
    grapho.events = evts + evts_trim          
    grapho.to_annot(annot, evt_name)
    
    return


def merge_events(annot, evt_name, chan = None, stage = None, segments = None):
    
    ''' Merges any events in annotations that are overlapping in time.
        # NOTE: This merges events with the SAME NAME ONLY. To merge 2 or more
                event types, use the function merge_2_events()
    '''
    
    evts = annot.get_events(name=evt_name, chan = chan, stage = stage)
    merged_events = []
    
    removed = []
    for e, event in enumerate(evts):
        if e not in removed:
            evts_trim = [item for i, item in enumerate(evts) if i != e]
        
            matches = ( [x for x in evts_trim if 
                         x['start'] <= event['start'] <= x['end']] +
                        [x for x in evts_trim if 
                         x['start'] <= event['end'] <= x['end']]
                      )
            
            #Find start and end times of merged events
            if len(matches) > 0:
                event['start'] = min([x['start'] for x in matches])
                event['end'] = max([x['end'] for x in matches])
                
                # Save to new event list
                merged_events.append(event)
            
            # Flag that these events have already been merged and removed
            for m in matches:
                index = [i for i, evt in enumerate(evts) if
                         evt['start'] == m['start'] if evt['end'] == m['end']]
                removed.append(index)
    
    events_trim = [evts[1] for i, x in enumerate (evts) if not i in removed]
    merged_events = merged_events + events_trim
    
    annot.remove_event_type(name = evt_name)
    grapho = graphoelement.Graphoelement()
    grapho.events = merged_events         
    grapho.to_annot(annot, evt_name)       


    return    

def merge_2_events(annot, events, new_event_name, chan = None, stage = None, 
                   segments = None):
    
    ''' Merges 2 different event types in annotations (irregardless of any 
            overlap). However, this will then also combine any events that 
            have a temporal overlap into the 1 event.
    '''
    # Get each event from annotations file
    evts = []
    for evt_name in events:
        evts += annot.get_events(name = evt_name, chan = chan, stage = stage)

    # Rename both events into a single event name
    evts = [{**d, "name": new_event_name} if "name" in d else d for d in evts]
    
    # Save new events to Annotations file
    grapho = graphoelement.Graphoelement()
    grapho.events = evts         
    grapho.to_annot(annot, new_event_name)       
    
    # Finally, combine any overlapping events
    merge_events(annot, new_event_name, chan, stage, segments)

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


def clean_annots(xml_dir, out_dir, rater, keep_evts = None, subs = 'all', 
                 sessions = 'all', logger = create_logger('Clean annotations')):
    
    '''
        Copies the annotations files from inside xml_dir into out_dir and 
        removes all annotations from the files, keeping only staging info.
        
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


def merge_epochs(epochs):
    
    ''' Function to merge epochs of the same stage together into larger segments,
        to extract the start and end times of periods of sleep stages. An 
        exception is granted for single epochs (30s) splitting a larger segment,
        in this case the single epoch is 'absorbed' by the larger segment.
        
        Parameters
        ----------
        epochs  - A list of dict corresponding to segments of data. Each dict 
                  must contain the following parameters: 'start', 'end', 'stage'.
        
        Returns
        -------
        A list of dict corresponding to segments of data. Each dict contains the 
        following parameters: 'start', 'end', 'stage'.
        
    '''
    

    merged = []
    current = epochs[0].copy()
    skip_next = False  # Flag to track if we are skipping an interrupting element

    for i in range(1, len(epochs)):
        entry = epochs[i]

        if skip_next:
            skip_next = False
            continue

        if entry["stage"] == current["stage"] and entry["start"] <= current["end"]:
            # Merge as usual
            current["end"] = max(current["end"], entry["end"])
        elif (
            i + 1 < len(epochs) and epochs[i + 1]["stage"] == current["stage"]
        ):  
            # If the next element matches the current stage, skip this one
            skip_next = True
            current["end"] = max(current["end"], epochs[i + 1]["end"])
        else:
            # Save the merged segment and move to the next one
            merged.append(current)
            current = entry.copy()

    # Append the last segment
    merged.append(current)
    return merged

                            

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
                

def laplacian_mne(data, channel, ref_chan, oREF = None, laplacian_rename=False, 
                  renames=None, montage='standard_alphabetic', 
                  logger = create_logger('Laplacian filter')):
    
    
    ch_names = list(data.axis['chan'][0])
    dig = mne.channels.make_standard_montage(montage)
    
    if oREF:
        try:
            dig.rename_channels({oREF:'_REF'}, allow_duplicates=False)
        except Exception as error:
            logger.error(error)
            logger.error("This is due to an incorrect naming of an online "
                         "reference 'oREF' (see user guide):")
            logger.info('')
            logger.info("Check documentation for how to set up channel data: "
                        "https://seapipe.readthedocs.io/en/latest/index.html")
            logger.info('-' * 10)
            return data, 1
    
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
    
    return data, 0


def notch_mne(data, channel, freq, oREF = None, rename=False, renames=None, 
              montage='standard_alphabetic', 
              logger = create_logger('Notch filter')):
    
    ch_names = list(data.axis['chan'][0])
    dig = mne.channels.make_standard_montage(montage)
    
    if rename:
        dig.rename_channels(renames, allow_duplicates=False)
        
    if oREF:
        try:
            dig.rename_channels({oREF:'_REF'}, allow_duplicates=False)
        except Exception as error:
            logger.error(error)
            logger.error("This is due to an incorrect naming of an online "
                         "reference 'oREF' (see user guide):")
            logger.info('')
            logger.info("Check documentation for how to set up channel data: "
                        "https://seapipe.readthedocs.io/en/latest/index.html")
            logger.info('-' * 10)
            return data, 1
    
    info = mne.create_info(ch_names, data.s_freq,verbose=40)
    mneobj = mne.io.RawArray(data.data[0],info,verbose=40)
    dic = [{x:'eeg' for x in mneobj.ch_names}]
    mneobj.set_channel_types(dic[0],verbose=40)
    mneobj.info.set_montage(dig,verbose=40)
    a = mneobj.pick(['eeg']).load_data()
    anotch = a.notch_filter(freq,verbose=40)
    data = anotch.get_data(picks=channel)
    
    return data, 0 

def notch_mne2(data, channel, oREF = None, rename=False, renames=None,
               montage = 'standard_alphabetic', 
               logger = create_logger('Notch filter')):
    
    ch_names = list(data.axis['chan'][0])
    dig = mne.channels.make_standard_montage(montage)
    
    if rename:
        dig.rename_channels(renames, allow_duplicates=False)
        
    if oREF:
        try:
            dig.rename_channels({oREF:'_REF'}, allow_duplicates=False)
        except Exception as error:
            logger.error(error)
            logger.error("This is due to an incorrect naming of an online "
                         "reference 'oREF' (see user guide):")
            logger.info('')
            logger.info("Check documentation for how to set up channel data: "
                        "https://seapipe.readthedocs.io/en/latest/index.html")
            logger.info('-' * 10)
            return data, 1
        
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
    
    return data, 0
    

def bandpass_mne(data, channel, highpass, lowpass, oREF = None, rename=False,
                 renames=None, montage='standard_alphabetic', 
                 logger = create_logger('Bandpass filter')):
    
    ch_names = list(data.axis['chan'][0])
    dig = mne.channels.make_standard_montage(montage)
    
    if rename:
        dig.rename_channels(renames, allow_duplicates=False)
        
    if oREF:
        try:
            dig.rename_channels({oREF:'_REF'}, allow_duplicates=False)
        except Exception as error:
            logger.error(error)
            logger.error("This is due to an incorrect naming of an online "
                         "reference 'oREF' (see user guide):")
            logger.info('')
            logger.info("Check documentation for how to set up channel data: "
                        "https://seapipe.readthedocs.io/en/latest/index.html")
            logger.info('-' * 10)
            return data, 1
        
    info = mne.create_info(ch_names, data.s_freq,verbose=40)
    mneobj = mne.io.RawArray(data.data[0],info,verbose=40)
    dic = [{x:'eeg' for x in mneobj.ch_names}]
    mneobj.set_channel_types(dic[0],verbose=40)
    mneobj.info.set_montage(dig,verbose=40)
    a = mneobj.pick(['eeg']).load_data()
    a_filt = a.filter(highpass, lowpass,verbose=40)
    data = a_filt.get_data(picks=channel)
    
    return data, 0

def csv_stage_import(edf_file, xml_file, hypno_file, rater):
    
    ''' This function creates a new annoations file and imports staging from a csv 
        file. The staging should be in the 'Alice' format, for further information 
        see: wonambi.attr.annotations.
    '''

    dataset = Dataset(edf_file)
    create_empty_annotations(xml_file, dataset)
    annot = Annotations(xml_file)
    annot.import_staging(hypno_file, 'alice', rater_name=rater, rec_start=dataset.header['start_time'])
    




def reconstruct_stitches(dat, stitches, s_freq, replacement = 0):
    orig_length = int(stitches[-1][1]*s_freq)  # Length of original array2
    
    # Step 1: Create a zero-filled array2
    dat_reconstructed = np.zeros(orig_length)
    dat_reconstructed.fill(replacement)
    
    # Step 2: Reinsert segments into dat_reconstructed
    idx = 0  # Index to track position in array1
    
    for start, end in stitches:
        start_idx = int(start * s_freq)
        end_idx = int(end * s_freq)
        segment_length = end_idx - start_idx
        
        # Insert the segment from array1
        dat_reconstructed[start_idx:end_idx] = dat[idx:idx + segment_length]
        
        # Update the array1 index tracker
        idx += segment_length
        
    return dat_reconstructed


def adap_bands_setup(self, adap_bands, frequency, subs, sessions, chan, ref_chan,
                     stage, cat, concat_cycle, cycle_idx, logger):

   
    if adap_bands == 'Fixed':
        flag = None
        if not frequency:
            frequency = (11,16)
    elif adap_bands == 'Manual':
        logger.debug(f"Checking for spectral peaks in {self.rootpath}/'tracking.tsv' ")
        flag = check_adap_bands(self.rootpath, subs, sessions, chan, logger)
        if flag == 'error':
            logger.critical('Spindle detection finished with ERRORS. See log for details.')
            return
        elif flag == 'review':
            logger.info('')
            logger.warning("Some spectral peak entries in 'tracking.tsv' are "
                           "inconsistent or missing. In these cases, detection will "
                           f"revert to fixed bands: {frequency[0]}-{frequency[1]}Hz")
            logger.info('')
    elif adap_bands == 'Auto': 
        if not frequency:
            frequency = (9,16)           
        self.track(subs, sessions, step = 'fooof', show = False, log = False)
        if not type(chan) == type(DataFrame()):
            logger.critical("For adap_bands = Auto, the argument 'chan' must "
                            "be 'None' and specfied in 'tracking.csv'")
            flag = 'error'
        else:
            flag, pk_chan, pk_sub, pk_ses = check_fooof(self, frequency, 
                                                              chan, ref_chan, 
                                                              stage, 
                                                              cat,
                                                              cycle_idx, 
                                                              logger)
        if flag == 'error':
            logger.critical('Error in reading channel names, check tracking sheet.')
            logger.info("Check documentation for how to set up channel names in tracking.tsv':")
            logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
            logger.info('-' * 10)
            return
        elif flag == 'review':
            logger.debug("Spectral peaks have not been found for all subs, "
                         "analysing the spectral parameters prior to spindle detection..")
            for (sub,ses) in zip(pk_sub, pk_ses):
                self.detect_spectral_peaks(subs = [sub], 
                                       sessions = [ses], 
                                       chan = pk_chan, 
                                       frequency = frequency,
                                       stage = stage, cycle_idx = cycle_idx,
                                       concat_cycle = concat_cycle, 
                                       concat_stage=True)    
    return frequency, flag


def infer_polarity(dset, annot, chan, ref_chan, cat = (1,1,1,1), evt_type = None, 
                   stage = None, cycle = None, logger = create_logger('Check polarity')):

    segments = fetch(dset, annot, cat, evt_type, stage, cycle)
    segments.read_data(chan, ref_chan)
    
    s_freq = dset.header['s_freq']
    
    X = segments.segments[0]['data'].data[0][0]
    
    XNormed = (X - X.mean())/(X.std())
    
    X_trans = transform_signal(XNormed, s_freq, 'double_sosbutter', 
                               method_opt={'freq':(0.1,4.5), 'order':3})
    
    ## Positive component of signal
    X_pos = copy.deepcopy(X_trans)
    X_pos[X_pos<0] = np.nan
    
    # Positive peaks
    peaks_pos = find_peaks(X_pos, height = np.nanstd(X_pos)) 
    pos_peaks_mean = np.nanmean(X_pos[peaks_pos[0]])
    
    # Positive mean
    pos_mean = np.nanmean(X_pos[X_pos>np.nanpercentile(X_pos, 99)])
    
    # Positive PSD
    result = flip_and_switch(X_pos)
    Sxx = periodogram(result, s_freq)
    pos_Sxx_mean = np.nanmean(Sxx[1][0:2000])
    
    ## Negative component of signal
    X_neg = copy.deepcopy(X_trans)
    X_neg[X_neg>0] = np.nan
    X_neg = X_neg * -1
    
    # Negative peaks
    peaks_neg = find_peaks(X_neg, height = np.nanstd(X_neg)) 
    neg_peaks_mean = np.nanmean(X_neg[peaks_neg[0]])
    
    # Negative mean
    neg_mean = np.nanmean(X_neg[X_neg>np.nanpercentile(X_neg, 99)])
    
    # Negative PSD
    result = flip_and_switch(X_neg)
    Sxx = periodogram(result, s_freq)
    neg_Sxx_mean = np.nanmean(Sxx[1][0:2000])
    
    # Compare Postive to Negative
    rating = ( 
              int(pos_peaks_mean>neg_peaks_mean) + 
              int(pos_Sxx_mean>neg_Sxx_mean) + 
              int(pos_mean > neg_mean) 
              )

    if rating > 0:
        logger.debug('Signal polarity appears to be reversed.')
        invert = True
    else:
        invert = False
        logger.debug('Signal polarity appears correct.')
    return invert


def flip_and_switch(data):
    
    # Step 1: Identify segments between NaNs
    segments = []
    current_segment = []
    
    for i, val in enumerate(data):
        if not np.isnan(val):
            current_segment.append(val)
        elif current_segment:
            segments.append(np.array(current_segment))
            current_segment = []
    
    # Add the last segment if it exists
    if current_segment:
        segments.append(np.array(current_segment))
    
    # Step 2: Flip every even-indexed segment (starting at index 1)
    for i in range(len(segments)):
        if i % 2 == 1:  # 0-based index: 1st, 3rd, etc.
            segments[i] = -segments[i]
    
    # Step 3: Stitch the segments together
    result = np.concatenate(segments)
    
    return result



