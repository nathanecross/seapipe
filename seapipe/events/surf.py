#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 19 14:43:23 2022

@author: nathancross
"""
from safepickle import dump
from os import listdir, mkdir, path
import math
import mne
from numpy import (argsort, asarray, linspace, tile, vstack, zeros)
import shutil
from wonambi import Dataset
from wonambi.trans import fetch
from wonambi.attr import Annotations

'''
    Sleep Undulation Related waveForm (SURF)
    
'''

def marker_to_annot(rec_dir, xml_dir, out_dir, part, visit, rater, chan, 
                    marker_name = None, file_type = 'edf'):
    
    '''
     Reads marker data from EDF directly and exports it into Annotations file 
    '''
    
    # Make output directory
    if not path.exists(out_dir):
        mkdir(out_dir)
    
    # Loop through records
    if isinstance(part, list):
        None
    elif part == 'all':
            part = listdir(xml_dir)
            part = [ p for p in part if not(p.startswith('.'))]
    else:
        print('')
        print("ERROR: 'part' must either be an array of subject ids or = 'all' ")
        print('')
    
    for i, p in enumerate(part):
        
        if not path.exists(out_dir + '/' + p):
            mkdir(out_dir + '/' + p)
        
        if visit == 'all':
            visit = listdir(xml_dir + '/' + p)
            visit = [x for x in visit if not(x.startswith('.'))]
        for j, vis in enumerate(visit): 
            if not path.exists(xml_dir + '/' + p + '/' + vis + '/'):
                print(f'WARNING: input folder missing for Subject {p}, visit {j}, skipping..')
                continue
            else:
                if not path.exists(out_dir + '/' + p + '/' + vis):
                    mkdir(out_dir + '/' + p + '/' + vis)
                rec_file = [s for s in listdir(rec_dir + '/' + p + '/' + vis) if ("."+file_type) in s if not s.startswith(".")]
                xml_file = [x for x in listdir(xml_dir + '/' + p + '/' + vis) if x.endswith('.xml') if not x.startswith(".")] 
                
                if len(xml_file) == 0:
                    print(f'WARNING: annotations does not exist for Subject {p}, visit {j} - check this. Skipping..')
                elif len(xml_file) >1:
                    print(f'WARNING: multiple annotations files exist for Subject {p}, visit {j} - check this. Skipping..')
                else:
                    print(f'Annotating embedded edf markers for Subject {p}, Visit {vis}')
                    
                    
                    # Copy annotations file before beginning
                    xdir = xml_dir + '/' + p + '/' + vis + '/' 
                    odir = out_dir + p + '/' + vis + '/'
                    backup_file = (f'{odir}{xml_file[0]}')
                    shutil.copy(xdir + xml_file[0], backup_file)
                    
                    # Open recording and annotations files
                    dset = Dataset(rec_dir + '/' + p + '/' + vis + '/' + rec_file[0]) 
                    annot = Annotations(backup_file, rater_name=rater)

                    # Read markers from edf
                    markers = dset.read_markers()
                    
                    # Filter out markers
                    if marker_name:
                        markers = [x for x in markers if x['name'] in marker_name]
                    
                    # Save markers to output annotations file
                    annot.add_events(markers,chan=chan)

def erp_analysis(rec_dir, xml_dir, out_dir, part, visit, cycle_idx, chan, ref_chan, oREF, 
                 rater, stage,polar, grp_name, cat, evt_type, filt, window,  
                 file_type ='edf'):

    '''
    
    Analyses ERPs
    
    '''
    
    
    ## Make output directory
    if not path.exists(out_dir):
        mkdir(out_dir)
    
    ## Get number of participants and visits 
    if isinstance(part, list):
        None
    elif part == 'all':
            part = listdir(xml_dir)
            part = [ p for p in part if not(p.startswith('.'))]
    else:
        print('')
        print("ERROR: 'part' must either be an array of subject ids or = 'all' ")
        print('')
        
    for i, p in enumerate(part):
        if visit == 'all':
            visit = listdir(xml_dir + '/' + p)
            visit = [x for x in visit if not(x.startswith('.'))]
        
        
    ## Prepare output array 
    if cycle_idx is None:
        cyc = [1]
    else:
        cyc = cycle_idx
    
    # Axis: sleep cycles
    if cat[0] == 0:
        cyccat = len(cyc)
    else:
        cyccat = 1
          
    # Axis: sleep stages
    if cat[1] == 0:
        stagecat = len(stage)
    else:
        stagecat = 1
    
    # Axis: event types
    if cat[3] == 0:
        evcat = len(evt_type)
    else:
        evcat = 1
    
    # Create output array    
    out = zeros((len(part),len(visit),cyccat,stagecat*evcat), dtype='object')  
           
    
    ## Loop through records
    for i, p in enumerate(part):
        
        if not path.exists(out_dir + '/' + p):
            mkdir(out_dir + '/' + p)
        
        if visit == 'all':
            visit = listdir(xml_dir + '/' + p)
            visit = [x for x in visit if not(x.startswith('.'))]
        for j, vis in enumerate(visit): 
            if not path.exists(xml_dir + '/' + p + '/' + vis + '/'):
                print(f'WARNING: input folder missing for Subject {p}, visit {j}, skipping..')
                continue
            else:
                if not path.exists(out_dir + '/' + p + '/' + vis):
                    mkdir(out_dir + '/' + p + '/' + vis)
                rec_file = [s for s in listdir(rec_dir + '/' + p + '/' + vis) if ("."+file_type) in s if not s.startswith(".")]
                xml_file = [x for x in listdir(xml_dir + '/' + p + '/' + vis) if x.endswith('.xml') if not x.startswith(".")] 
                if len(xml_file) == 0:
                    print(f'WARNING: annotations does not exist for Subject {p}, visit {j} - check this. Skipping..')
                elif len(xml_file) >1:
                    print(f'WARNING: multiple annotations files exist for Subject {p}, visit {j} - check this. Skipping..')
                else:
                    # Open recording and annotations files
                    dset = Dataset(rec_dir + '/' + p + '/' + vis + '/' + rec_file[0]) 
                    annot = Annotations(xml_dir + '/' + p + '/' + vis + '/' + xml_file[0], rater_name=rater)
                    
                    # Check polarity of recording
                    if isinstance(polar, list):
                        polarity = polar[i]
                    else:
                        polarity = 'normal'
                       
                    # Get sleep cycles
                    if cycle_idx is not None:
                        all_cycles = annot.get_cycles()
                        scycle = [all_cycles[i - 1] for i in cycle_idx if i <= len(all_cycles)]
                        
                    else:
                        scycle = [None]  

                    chan_full = []
                    for k, ch in enumerate(chan):
                        chan_full.append(ch + ' (' + grp_name + ')')
                            
                    for l, cyc in enumerate(scycle):
                        
                        ## Select and read data
                        print(f'Reading data for {p}, cycle {l+1}')
                        segments = fetch(dset, annot, cat=cat, chan_full=chan_full, 
                                         cycle=[cyc], evt_type=evt_type, stage=stage,
                                         )
                        
                        events = zeros((1,3))
                        event_dict = {}
                        
                        ## Format events
                        for m, seg in enumerate(segments):
                            evtname = seg['name']
                            stg = seg['stage']
                            event_dict[evtname+stg]=m
                            ev = asarray(list(set([x[0]*dset.header['s_freq'] for x in seg['times']])))
                            ev2 = tile(m,(len(ev)))
                            evts = vstack((ev,ev2,ev2)).T
                            events = vstack((events,evts))

                        events = events.astype(int)[1:-1,:]
                        events = events[argsort(events[:, 0])]
                        
                        
                        
                        ## Format data for analysis
                        print('Preparing p, vis for ERP analysis')
                        print(r"""
                              
                                  |
                                  | /\ 
                                  |/  \  _
                              uV2 |   | / \
                                  |   \/   ^-___-__
                                  |________________
                                         (sec)
                              
                              """)
                        
                        achan = chan + ref_chan
                        d = dset.read_data(chan=achan)
                        
                        # Adjust polarity
                        if polarity == 'normal':
                            data = d.data[0]
                        else:
                            data = d.data[0]*-1
                        
                        # Convert data to mne format & filter
                        ch_names = list(d.axis['chan'][0])
                        dig = mne.channels.make_standard_montage('standard_alphabetic')
                        dig.rename_channels({oREF:'_REF'}, allow_duplicates=False)
                        info = mne.create_info(ch_names, d.s_freq,verbose=False)
                        mneobj = mne.io.RawArray(data,info,verbose=False)
                        dic = [{x:'eeg' for x in mneobj.ch_names}]
                        mneobj.set_channel_types(dic[0],verbose=False)
                        mneobj.info.set_montage(dig,verbose=False)
                        a = mneobj.pick_types(eeg=True).load_data()
                        a.set_eeg_reference(ref_chan,verbose=False)
                        a.filter(l_freq=filt[0], h_freq=filt[1])    
                        
                        # Save epochs around stimulus event (ERP)
                        epochs = mne.Epochs(a, events, event_id=event_dict, 
                                             tmin=window[0], tmax=window[1], preload=True)

                        # Save individual file
                        with open(out_dir + '/' + p + '/' + vis + '/' + p + '_' + 
                                  '_epochs.p', 'wb') as f:
                             dump(epochs, f)

                        # Save individual epochs to master array
                        for av,avalue in enumerate(event_dict):
                            out[i,j,l,av] = epochs[avalue].average()
        
        # Save master array
        with open(out_dir + f'cat_{cat[0]}{cat[1]}{cat[2]}{cat[3]}_evoked.p', 'wb') as f:
             dump(out, f)
             
    return out
                        
                          
                            
def insert_virtual_markers(rec_dir, xml_dir, out_dir, part, visit, rater, chan, name, reptime,
                           file_type ='edf'):
    
    '''
     Inserts virtual markers into Annotations file 
    '''              
    
    # Make output directory
    if not path.exists(out_dir):
        mkdir(out_dir)
    
    # Loop through records
    if isinstance(part, list):
        None
    elif part == 'all':
            part = listdir(xml_dir)
            part = [ p for p in part if not(p.startswith('.'))]
    else:
        print('')
        print("ERROR: 'part' must either be an array of subject ids or = 'all' ")
        print('')
    
    for i, p in enumerate(part):
        
        if not path.exists(out_dir + '/' + p):
            mkdir(out_dir + '/' + p)
        
        if visit == 'all':
            visit = listdir(xml_dir + '/' + p)
            visit = [x for x in visit if not(x.startswith('.'))]
        for j, vis in enumerate(visit): 
            if not path.exists(xml_dir + '/' + p + '/' + vis + '/'):
                print(f'WARNING: input folder missing for Subject {p}, visit {j}, skipping..')
                continue
            else:
                if not path.exists(out_dir + '/' + p + '/' + vis):
                    mkdir(out_dir + '/' + p + '/' + vis)
                #rec_file = [s for s in listdir(rec_dir + '/' + p + '/' + vis) if ("."+file_type) in s if not s.startswith(".")]
                xml_file = [x for x in listdir(xml_dir + '/' + p + '/' + vis) if x.endswith('.xml') if not x.startswith(".")] 
                
                if len(xml_file) == 0:
                    print(f'WARNING: annotations does not exist for Subject {p}, visit {j} - check this. Skipping..')
                elif len(xml_file) >1:
                    print(f'WARNING: multiple annotations files exist for Subject {p}, visit {j} - check this. Skipping..')
                else:
                    print(f'Creating virtual markers for Subject {p}, Visit {vis}')
                    
                    
                    # Copy annotations file before beginning
                    xdir = xml_dir + '/' + p + '/' + vis + '/' 
                    odir = out_dir + p + '/' + vis + '/'
                    backup_file = (f'{odir}{xml_file[0]}')
                    shutil.copy(xdir + xml_file[0], backup_file)
                    
                    # Open recording and annotations files
                    #dset = Dataset(rec_dir + '/' + p + '/' + vis + '/' + rec_file[0]) 
                    annot = Annotations(backup_file, rater_name=rater)  
                         
                    # Get start time of first sleep stage
                    epochs = annot.get_epochs()   
                    first_epoch = [x for x in epochs if x['stage'] in ['NREM1','NREM2','NREM3','REM']][0]
                    onset = int(first_epoch['start'])
                    
                    # Get end time of last sleep stage
                    last_epoch = [x for x in epochs if x['stage'] in ['NREM1','NREM2','NREM3','REM']][-1]
                    offset = int(last_epoch['end'])
                    
                    
                    # Create number of events (based on repition time)
                    num_evts = (offset-onset)/reptime
                    
                    # Create list of start times for events
                    if num_evts.is_integer():
                        starts = linspace(onset, offset, int(num_evts))
                    else:
                        num_evts = math.ceil(num_evts)
                        offset = int(onset + (num_evts*reptime))
                        starts = linspace(onset, offset, num_evts)
                    
                    # Create list of evts  
                    evts = []
                    for x, st in enumerate(starts):
                        evts.append({'name': name,
                         'start': st,
                         'end': st + 0.5,
                         'chan': chan,
                         'stage': '',
                         'quality': 'Good'})

                    # Add events to annotations file
                    annot.add_events(evts)
                    
                    
                    
def marker_to_annotAP(rec_dir, xml_dir, out_dir, part, visit, rater, chan, marker_name=None, file_type='edf'):
    
    '''
     Reads marker data from EDF directly and exports it into Annotations file 
    '''
    
    # Make output directory
    if not path.exists(out_dir):
        mkdir(out_dir)
    
    # Loop through records
    if isinstance(part, list):
        None
    elif part == 'all':
            part = listdir(xml_dir)
            part = [ p for p in part if not(p.startswith('.'))]
    else:
        print('')
        print("ERROR: 'part' must either be an array of subject ids or = 'all' ")
        print('')
    
    for i, p in enumerate(part):
        
        if not path.exists(out_dir + '/' + p):
            mkdir(out_dir + '/' + p)
        
        if visit == 'all':
            visit = listdir(xml_dir + '/' + p)
            visit = [x for x in visit if not(x.startswith('.'))]
        for j, vis in enumerate(visit): 
            if not path.exists(xml_dir + '/' + p + '/' + vis + '/'):
                print(f'WARNING: input folder missing for Subject {p}, visit {j}, skipping..')
                continue
            else:
                if not path.exists(out_dir + '/' + p + '/' + vis):
                    mkdir(out_dir + '/' + p + '/' + vis)
                rec_file = [s for s in listdir(rec_dir + '/' + p + '/' + vis) if ("."+file_type) in s if not s.startswith(".")]
                xml_file = [x for x in listdir(xml_dir + '/' + p + '/' + vis) if x.endswith('.xml') if not x.startswith(".")] 
                
                if len(xml_file) == 0:
                    print(f'WARNING: annotations does not exist for Subject {p}, visit {j} - check this. Skipping..')
                elif len(xml_file) >1:
                    print(f'WARNING: multiple annotations files exist for Subject {p}, visit {j} - check this. Skipping..')
                else:
                    print(f'Annotating embedded edf markers for Subject {p}, Visit {vis}')
                    
                    
                    # Copy annotations file before beginning
                    xdir = xml_dir + '/' + p + '/' + vis + '/' 
                    odir = out_dir + p + '/' + vis + '/'
                    backup_file = (f'{odir}{xml_file[0]}')
                    shutil.copy(xdir + xml_file[0], backup_file)
                    
                    # Open recording and annotations files
                    dset = Dataset(rec_dir + '/' + p + '/' + vis + '/' + rec_file[0]) 
                    annot = Annotations(backup_file, rater_name=rater)

                    # Read markers from edf
                    markers = dset.read_markers()
                    
                    # Filter out markers
                    if marker_name:
                        markers = [x for x in markers if x['name'] in marker_name]
                    
                    # Select every 2nd element    
                    markers = markers[1::2]
                    
                    # Save markers to output annotations file
                    annot.add_events(markers,chan=chan)                    
                                  
                                    
                                    
                                    
                                    
                                    