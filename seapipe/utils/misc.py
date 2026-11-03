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
from scipy.signal import find_peaks,  periodogram
import shutil
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
                if chan in x['chan']]
    elif not chan:
        evts = [x for x in annot.get_events(name=evt_name) 
                if stage in x['stage']]
    else:
        evts = [x for x in annot.get_events(name=evt_name) 
                if stage in x['stage'] if chan in x['chan']]
    
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
    
    if chan:
        evts = [x for x in annot.get_events(name = evt_name) if chan not in x['chan']]
    else:
        evts = []
    
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
    
    events_trim = [evts[i] for i, x in enumerate (evts) if not i in removed]
    merged_events = merged_events + events_trim
    
    
    if chan:
        evts = [x for x in annot.get_events(name = evt_name) if chan not in x['chan']]
    else:
        evts = []
    
    annot.remove_event_type(name = evt_name)
    grapho = graphoelement.Graphoelement()
    grapho.events = evts + merged_events         
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




def _as_event_name_list(event_names):
    if event_names is None:
        return []
    if isinstance(event_names, str):
        return [event_names]
    return list(event_names)


def _pair_sensible_event_windows(start_events, end_events):
    start_times = sorted(evt['end'] for evt in start_events)
    end_times = sorted(evt['start'] for evt in end_events)

    windows = []
    end_idx = 0
    for start_time in start_times:
        while end_idx < len(end_times) and end_times[end_idx] <= start_time:
            end_idx += 1

        if end_idx >= len(end_times):
            break

        windows.append((start_time, end_times[end_idx]))
        end_idx += 1

    return windows


def _merge_time_windows(windows):
    if not windows:
        return []

    windows = sorted(windows, key=itemgetter(0, 1))
    merged = [list(windows[0])]

    for start, end in windows[1:]:
        last = merged[-1]
        if start <= last[1]:
            last[1] = max(last[1], end)
        else:
            merged.append([start, end])

    return [tuple(win) for win in merged]


def _event_overlap_ratio(event, windows):
    duration = event['end'] - event['start']
    if duration <= 0:
        for start, end in windows:
            if start <= event['start'] <= end:
                return 1.0
        return 0.0

    overlap = 0.0
    for start, end in windows:
        overlap += max(0.0, min(event['end'], end) - max(event['start'], start))

    return min(overlap / duration, 1.0)


def mask_events_annots(xml_dir, out_dir, event_name, segs_start, segs_end,
                       subs='all', sessions='all', rater=None,
                       min_overlap=0.5, logger=create_logger('Mask events')):

    """
        Copies annotations files from xml_dir to out_dir and masks one event
        type to windows defined by paired start/end marker events.
    """

    if not 0 <= min_overlap <= 1:
        raise ValueError("'min_overlap' must be between 0 and 1.")

    start_names = _as_event_name_list(segs_start)
    end_names = _as_event_name_list(segs_end)
    if not start_names or not end_names:
        logger.error("Both 'segs_start' and 'segs_end' must be provided.")
        return

    logger.debug('')
    logger.debug(f'Time start: {datetime.now()}')
    logger.debug(f"Masking '{event_name}' in annotations from: {xml_dir}")
    logger.debug(f"Using segment markers start={start_names}, end={end_names}")
    logger.debug(f'Saving masked annotations to: {out_dir}')
    logger.debug(f'Minimum event overlap retained: {min_overlap}')
    logger.debug('')

    if isinstance(subs, list):
        subs_list = list(subs)
    elif subs == 'all':
        subs_list = [x for x in listdir(xml_dir) if '.' not in x]
    elif isinstance(subs, str):
        subs_list = [subs]
    else:
        logger.error("'subs' must be a list, a subject id, or 'all'.")
        return

    subs_list.sort()
    for sub in subs_list:
        sub_dir = f'{xml_dir}/{sub}'
        if not path.exists(sub_dir):
            logger.warning(f'{sub_dir} does not exist - skipping {sub}.')
            continue

        if sessions == 'all':
            sub_sessions = [x for x in listdir(sub_dir) if '.' not in x]
        elif isinstance(sessions, list):
            sub_sessions = list(sessions)
        else:
            sub_sessions = [sessions]

        sub_sessions.sort()
        for ses in sub_sessions:
            xdir = f'{xml_dir}/{sub}/{ses}/'
            if not path.exists(xdir):
                logger.warning(f'{xdir} does not exist - skipping {sub}, {ses}.')
                continue

            xml_files = sorted([
                x for x in listdir(xdir)
                if x.endswith('.xml') and not x.startswith('.')
            ])
            if len(xml_files) == 0:
                logger.warning(f'No xml files found for {sub}, {ses} - skipping.')
                continue

            if not path.exists(out_dir):
                mkdir(out_dir)
            if not path.exists(f'{out_dir}/{sub}'):
                mkdir(f'{out_dir}/{sub}')
            if not path.exists(f'{out_dir}/{sub}/{ses}'):
                mkdir(f'{out_dir}/{sub}/{ses}')

            for xml_file in xml_files:
                backup_file = f'{out_dir}/{sub}/{ses}/{xml_file}'
                if path.exists(backup_file):
                    logger.warning(f'{backup_file} already exists - overwriting copied annotations.')
                if not xdir + xml_file == backup_file:
                    shutil.copy(xdir + xml_file, backup_file)

                annot = Annotations(backup_file, rater_name=rater)

                start_events = []
                for name in start_names:
                    start_events.extend(annot.get_events(name=name))

                end_events = []
                for name in end_names:
                    end_events.extend(annot.get_events(name=name))

                windows = _pair_sensible_event_windows(start_events, end_events)
                ignored_starts = max(0, len(start_events) - len(windows))
                ignored_ends = max(0, len(end_events) - len(windows))
                if ignored_starts or ignored_ends:
                    logger.warning(
                        f'Ignored unmatched segment markers for {sub}, {ses} in {xml_file}: '
                        f'{ignored_starts} start, {ignored_ends} end.'
                    )

                windows = _merge_time_windows(windows)
                if not windows:
                    logger.warning(
                        f'No valid segment windows found for {sub}, {ses} in {xml_file} - '
                        f"leaving '{event_name}' unchanged."
                    )
                    continue

                events = annot.get_events(name=event_name)
                if len(events) == 0:
                    logger.warning(f"No '{event_name}' events found for {sub}, {ses} in {xml_file}.")
                    continue

                kept_events = [
                    evt for evt in events
                    if _event_overlap_ratio(evt, windows) >= min_overlap
                ]

                annot.remove_event_type(name=event_name)
                if kept_events:
                    grapho = graphoelement.Graphoelement()
                    grapho.events = kept_events
                    grapho.to_annot(annot, event_name)

                logger.info(
                    f"Masked '{event_name}' for {sub}, {ses} in {xml_file}: "
                    f'kept {len(kept_events)}/{len(events)} events across {len(windows)} windows.'
                )

    return

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
                for x in newevents:
                    x.update(d1)

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
                

def _laplacian_montage(ch_names, montage, oREF, rename, renames):
    """Assign montage coordinates to recording labels without renaming samples.

    ``renames`` maps montage names to recording names, e.g. {'FpZ': 'EEG'}.
    Montage lookup ignores case; recording labels always remain exact.
    """
    dig = mne.channels.make_standard_montage(montage)
    positions = dig.get_positions()
    source_positions = positions['ch_pos']

    def resolve(label, available, kind):
        if label in available:
            return label
        matches = [name for name in available if name.casefold() == label.casefold()]
        if len(matches) != 1:
            raise ValueError(f'{kind} label {label!r} has {len(matches)} matches; '
                             'use an unambiguous montage-to-recording mapping')
        return matches[0]

    sources = {}
    if rename:
        if not isinstance(renames, dict) or not renames:
            raise ValueError('laplacian_rename=True requires a nonempty renames dictionary '
                             'mapping montage names to recording names, e.g. {"FpZ": "EEG"}')
        for source, target in renames.items():
            if not isinstance(source, str) or not isinstance(target, str):
                raise ValueError('renames must map string montage names to string recording names')
            target = resolve(target, ch_names, 'Recording')
            source = resolve(source, source_positions, 'Montage')
            if target in sources and sources[target] != source:
                raise ValueError(f'Conflicting montage positions for recording channel {target!r}')
            sources[target] = source
    if oREF and '_REF' in ch_names:
        if not isinstance(oREF, str):
            raise ValueError('oREF must be a montage electrode name')
        reference = resolve(oREF, source_positions, 'Online reference')
        if '_REF' in sources and sources['_REF'] != reference:
            raise ValueError('oREF and renames assign different positions to _REF')
        sources['_REF'] = reference
    mapped_positions = {}
    for target in ch_names:
        source = sources.get(target)
        if source is None:
            source = resolve(target, source_positions, 'Montage')
        sources[target] = source
        coord = source_positions[source]
        if not np.isfinite(coord).all() or np.allclose(coord, 0):
            raise ValueError(f'Invalid montage position for {target!r} ({source!r})')
        for previous, previous_coord in mapped_positions.items():
            if np.allclose(coord, previous_coord, rtol=0, atol=1e-8):
                raise ValueError(f'Channels {previous!r} and {target!r} share a montage '
                                 'position; Laplacian electrodes must be distinct')
        mapped_positions[target] = coord
    # Preserve montage point order so automatic sphere fitting is independent
    # of recording-label aliases and the requested channel order.
    source_order = {name: i for i, name in enumerate(source_positions)}
    mapped_positions = {name: mapped_positions[name] for name in
                        sorted(ch_names, key=lambda name: source_order[sources[name]])}
    return mne.channels.make_dig_montage(
        ch_pos=mapped_positions, nasion=positions['nasion'], lpa=positions['lpa'],
        rpa=positions['rpa'], coord_frame=positions['coord_frame'])


def laplacian_mne(data, channel, ref_chan=None, oREF=None, laplacian_rename=False,
                  renames=None, montage='standard_alphabetic',
                  logger=create_logger('Laplacian filter')):
    """Compute CSD using all loaded electrodes and return only requested rows.

    Coordinates match case-insensitively; renames maps montage to recording
    labels (e.g. {'FpZ': 'EEG'}). Unknown, ambiguous or duplicate positions
    raise a logged ValueError. ``oREF`` locates a loaded ``_REF`` channel.
    ``ref_chan`` is retained for compatibility, but no second re-reference is
    applied: read_data already references the signal, and CSD is reference-free.
    Input samples and labels are never modified.
    """
    try:
        raw, picks = _prepare_mne_data(data, channel)
        if len(raw.ch_names) < 4:
            raise ValueError('Laplacian filtering needs at least four positioned electrodes '
                             'for automatic sphere fitting')
        dig = _laplacian_montage(raw.ch_names, montage, oREF, laplacian_rename, renames)
        raw.set_montage(dig, verbose=40)
        csd = mne.preprocessing.compute_current_source_density(raw, verbose=40)
        return csd.get_data(picks=picks), 0
    except (TypeError, ValueError, RuntimeError, np.linalg.LinAlgError) as error:
        message = (f'Laplacian filter failed for channel={channel!r}, montage={montage!r}, '
                   f'oREF={oREF!r}, renames={renames!r}: {error}')
        logger.error(message)
        raise ValueError(message) from error


def _prepare_mne_data(data, channel):
    """Validate one trial and return a private RawArray plus exact row picks."""
    ch_names = list(data.axis['chan'][0])
    if not ch_names or any(not isinstance(ch, str) or not ch for ch in ch_names):
        raise ValueError('Channel labels must be nonempty strings')
    if len(set(ch_names)) != len(ch_names):
        raise ValueError(f'Duplicate channel labels in recording: {ch_names!r}')
    values = np.asarray(data.data[0])
    if values.ndim != 2 or values.shape[0] != len(ch_names) or values.shape[1] == 0:
        raise ValueError('Expected a nonempty channels-by-samples array with '
                         f'{len(ch_names)} rows; got shape {values.shape}')
    if np.iscomplexobj(values) or not np.isfinite(values).all():
        raise ValueError('Signal must contain only finite, real samples')
    sfreq = float(data.s_freq)
    if not np.isfinite(sfreq) or sfreq <= 0:
        raise ValueError(f'Sampling frequency must be finite and positive; got {sfreq}')
    selected = [channel] if isinstance(channel, str) else list(channel)
    if not selected or any(not isinstance(ch, str) for ch in selected):
        raise ValueError('channel must be a channel name or a nonempty list of names')
    if len(set(selected)) != len(selected):
        raise ValueError(f'Duplicate requested channels: {selected!r}')
    missing = [ch for ch in selected if ch not in ch_names]
    if missing:
        raise ValueError(f'Requested channels {missing!r} are absent from recording '
                         f'{ch_names!r}; use recording labels, not montage/output aliases')
    info = mne.create_info(ch_names, sfreq, ch_types='eeg', verbose=40)
    # Copy explicitly: MNE filters in place, while callers own the input data.
    raw = mne.io.RawArray(np.array(values, dtype=float, copy=True), info, verbose=40)
    return raw, [ch_names.index(ch) for ch in selected]


def _temporal_filter_mne(data, channel, operation, logger, **filter_kwargs):
    """Filter one trial without requiring electrode coordinates or changing labels.

    Fail before returning any data on invalid input. Returning the original
    Wonambi object with an error flag is unsafe: callers assign the first return
    value directly to ``data.data[0]`` and may overwrite the flag later.
    """
    try:
        raw, picks = _prepare_mne_data(data, channel)
        sfreq = raw.info['sfreq']
        if operation == 'notch_filter' and filter_kwargs.get('method') != 'spectrum_fit':
            freqs = np.atleast_1d(np.asarray(filter_kwargs['freqs'], dtype=float))
            if (freqs.ndim != 1 or not freqs.size or not np.isfinite(freqs).all()
                    or np.any(freqs <= 0) or np.any(freqs >= sfreq / 2)):
                raise ValueError('Notch frequencies must be finite and strictly between '
                                 f'0 and Nyquist ({sfreq / 2:g} Hz); got '
                                 f"{filter_kwargs['freqs']!r}")
        getattr(raw, operation)(**filter_kwargs, verbose=40)
        return raw.get_data(picks=picks), 0
    except (TypeError, ValueError, RuntimeError) as error:
        message = (f'MNE {operation} failed for channel={channel!r}, '
                   f'options={filter_kwargs!r}: {error}')
        logger.error(message)
        raise ValueError(message) from error


def notch_mne(data, channel, freq, logger=create_logger('Notch filter')):
    """Notch-filter one trial and return ``(selected_samples, 0)``.

    ``channel`` uses recording labels; selected rows retain the requested order.
    Temporal filtering needs no electrode positions or montage configuration.
    In particular, FpZ/Fpz and EEG aliases do not need montage renaming here.
    Invalid inputs raise a logged ValueError rather than returning unfiltered data.
    """
    return _temporal_filter_mne(data, channel, 'notch_filter', logger, freqs=freq)


def notch_mne2(data, channel, logger=create_logger('Notch filter')):
    """Remove fitted sinusoids without a montage; see ``notch_mne`` for inputs.

    Uses MNE's spectrum-fit detection, rather than a fixed list of harmonics.
    """
    return _temporal_filter_mne(data, channel, 'notch_filter', logger,
                               freqs=None, filter_length='auto', notch_widths=None,
                               trans_bandwidth=1, method='spectrum_fit')


def bandpass_mne(data, channel, highpass, lowpass,
                 logger=create_logger('Bandpass filter')):
    """Bandpass-filter without a montage; see ``notch_mne`` for input conventions.

    Cutoff validation is delegated to MNE and failures are logged with context.
    """
    return _temporal_filter_mne(data, channel, 'filter', logger,
                               l_freq=highpass, h_freq=lowpass)

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
                     stage, concat_stage, concat_cycle, cycle_idx, logger):
    
    if adap_bands not in ['Fixed', 'Auto', 'Manual']:
        logger.critical("The argument 'adap_bands' must be either one of"
                        f"'Fixed', 'Auto', 'Manual', not: '{adap_bands}'")
        flag = 'error'
   
    if adap_bands == 'Fixed':
        flag = 'None'
        if not frequency:
            frequency = (11,16)
    elif adap_bands == 'Manual':
        logger.debug(f"Checking for spectral peaks in {self.rootpath}/'tracking.tsv' ")
        flag = check_adap_bands(self.rootpath, subs, sessions, chan, logger)
        if flag == 'error':
            logger.critical('Spindle detection finished with ERRORS. See log for details.')
            
        elif flag == 'review':
            logger.info('')
            logger.warning("Some spectral peak entries in 'tracking.tsv' are "
                           "inconsistent or missing. In these cases, detection will "
                           f"revert to fixed bands: {frequency[0]}-{frequency[1]}Hz")
            logger.info('')
    elif adap_bands == 'Auto': 
        if not frequency:
            frequency = (9,16) 
            flag = 'None'
        self.track(subs, sessions, step = 'fooof', show = False, log = False)
        if not type(chan) == type(DataFrame()):
            logger.critical("For adap_bands = Auto, the argument 'chan' must "
                            "be 'None' and specfied in 'tracking.csv'")
            flag = 'error'
        else:
            cat = (int(concat_cycle), int(concat_stage), 1, 1)
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


def infer_polarity(dset, annot, chan, ref_chan, cat=(1,1,1,1), evt_type=None,
                   stage=None, cycle=None, logger=create_logger('Check polarity'),
                   *, reject_artf=False, window_seconds=30, min_windows=5,
                   asymmetry_ratio=1.2, agreement=0.8):
    """Conservatively infer inversion from persistent positive slow-wave asymmetry.

    Both mean peak magnitude and upper-tail magnitude must exceed their negative
    counterparts by 20% in at least 80% of usable 30-second windows (minimum five).
    These are heuristic defaults, not calibrated confidence probabilities.
    Uncertain/insufficient data leave polarity unchanged and are logged.
    """
    if (window_seconds <= 0 or min_windows < 1 or asymmetry_ratio <= 1
            or not 0.5 < agreement <= 1):
        raise ValueError('Invalid polarity inference thresholds')
    segments = fetch(dset, annot, cat=cat, evt_type=evt_type, stage=stage,
                     cycle=cycle, reject_artf=reject_artf)
    segments.read_data(chan, ref_chan)
    s_freq = dset.header['s_freq']
    size = int(window_seconds * s_freq)
    if size < 2:
        raise ValueError('Polarity window must contain at least two samples')
    votes = []
    for segment in segments.segments:
        signal = np.asarray(segment['data'].data[0][0], dtype=float)
        for start in range(0, len(signal) - size + 1, size):
            raw = signal[start:start + size]
            if not np.all(np.isfinite(raw)):
                continue
            centered = raw - np.median(raw)
            scale = 1.4826 * np.median(np.abs(centered))
            # Reject flat windows and conspicuous transients before filtering.
            if scale <= 0 or np.max(np.abs(centered)) > 8 * scale:
                continue
            filtered = transform_signal(centered / scale, s_freq,
                                        'double_sosbutter',
                                        method_opt={'freq': (0.1, 4.5), 'order': 3})
            # Exclude filter edges from the asymmetry measurements.
            edge = int(2 * s_freq)
            filtered = filtered[edge:-edge] if edge else filtered
            if not len(filtered) or not np.all(np.isfinite(filtered)):
                continue
            measures = []
            for component in (filtered, -filtered):
                positive = component[component > 0]
                peaks, _ = find_peaks(component, height=np.std(component))
                if len(positive) < 10 or len(peaks) < 3:
                    break
                tail = positive[positive >= np.percentile(positive, 99)]
                measures.append((np.mean(component[peaks]), np.mean(tail)))
            if len(measures) != 2:
                continue
            positive, negative = map(np.asarray, measures)
            if np.any(negative <= 0) or np.any(positive <= 0):
                continue
            ratios = positive / negative
            votes.append(1 if np.all(ratios >= asymmetry_ratio)
                         else -1 if np.all(ratios <= 1 / asymmetry_ratio) else 0)
    if len(votes) >= min_windows:
        support = np.mean(np.asarray(votes) == 1)
        if support >= agreement:
            logger.debug(f'Signal polarity appears reversed: {support:.0%} of '
                         f'{len(votes)} usable windows support inversion.')
            return True
        if np.mean(np.asarray(votes) == -1) >= agreement:
            logger.debug('Signal polarity appears correct; leaving unchanged.')
            return False
    logger.warning(f'Polarity uncertain ({len(votes)} usable windows); leaving unchanged.')
    return False


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


def check_data_length(data, duration, logger = create_logger('Data check')):
    
    if any([data(trial=x).shape[1]<= duration*data.s_freq for 
            x in range(data.number_of('trial'))]):
        logger.warning('Length of data too short!!')
