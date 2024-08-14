#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 15 11:33:16 2022

@author: nathancross
"""
import mne
import copy
from operator import itemgetter
from os import listdir, path, mkdir
from datetime import datetime
from wonambi import Dataset, graphoelement
from wonambi.attr.annotations import Annotations, create_empty_annotations
import shutil


def remove_evts(xml_dir, out_dir, rater, evt_name = None, 
                part = 'all', visit = 'all'):
    
    print('')
    print(f'Time start: {datetime.now()}')
    print(f'Removing events from files in directory {xml_dir}')  
    print(f'Events = {evt_name}')
    print(f'Saving new annotations files to directory {out_dir}')
    print('')
    
    # Loop through subjects
    if isinstance(part, list):
        None
    elif part == 'all':
            part = listdir(xml_dir)
            part = [ p for p in part if not '.' in p]
    else:
        print("ERROR: 'part' must either be an array of subject ids or = 'all' **CHECK**")
    
    if evt_name == None:
        evt_name = ['spindle']
    
    part.sort()                
    for i, p in enumerate(part):
        # Loop through visits
        if visit == 'all':
            visit = listdir(xml_dir + '/' + '/' + p)
            visit = [x for x in visit if not '.' in x]
        
        print(f'Removing events for Subject {p}..')
        
        visit.sort()
        for v, vis in enumerate(visit): 
            
            # Define files
            xdir = xml_dir + '/' + p + '/' + vis + '/'
            xml_file = [x for x in listdir(xdir) if x.endswith('.xml') if not x.startswith('.')] 
            
            if len(xml_file) == 0:                
                print(f'WARNING: no xml files found for Subject {p}, skipping..')
                
            else:
            
                ## Copy annotations file before beginning
                if not path.exists(out_dir):
                    mkdir(out_dir)
                if not path.exists(out_dir + p ):
                    mkdir(out_dir + p)
                if not path.exists(out_dir + p + '/' + vis):
                    mkdir(out_dir + p + '/' + vis)
                backup = out_dir + p + '/' + vis + '/'
                backup_file = (f'{backup}{xml_file[0]}')
                shutil.copy(xdir + xml_file[0], backup_file)
                            
                # Import Annotations file
                annot = Annotations(backup_file, rater_name=rater)
                
                ### WHOLE NIGHT ###
                # Select and read data
                print('Reading data for ' + p + ', visit ' + vis )
                
                for e, ev in enumerate(evt_name):
                    annot.remove_event_type(name=ev)
                
                
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


def remove_duplicate_evts(annot, evt_name, chan, stage = None):
    
    '''Workaround function because wonambi.attr.annotations.remove_event()
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
    
    '''This function creates a new annoations file and imports staging from a csv 
        file. The staging should be in the 'Alice' format, for further information 
        see: wonambi.attr.annotations.
    '''

    dataset = Dataset(edf_file)
    create_empty_annotations(xml_file, dataset)
    annot = Annotations(xml_file)
    annot.import_staging(hypno_file, 'alice', rater_name=rater, rec_start=dataset.header['start_time'])
    
    
    
    
    
    