# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 10:34:53 2021

@author: labottdv
"""
from datetime import datetime
from itertools import product
from numpy import (asarray, concatenate, empty, nan, ndarray, sum, zeros)
from os import listdir, mkdir, path
from pandas import concat, DataFrame, read_csv
from safepickle import dump
import shutil
import sys
import traceback
from wonambi import Dataset
from wonambi.attr import Annotations
from wonambi.detect import consensus, DetectSpindle
from wonambi.trans import (fetch, get_times)
from wonambi.trans.analyze import event_params, export_event_params


def whale_it(rec_dir, xml_dir, out_dir, method, chan, ref_chan, rater, cat, stage, 
             grp_name, cycle_idx=None, frequency=(11, 16), adap_bands=False, duration=(0.5, 3), 
             part='all', visit='all'):
    
    '''
    Runs one (or multiple) automatic spindle detection algorithms and saves the 
    detected events to a new annotations file.
    
    INPUTS:
        rec_dir ->   (string) Path to directory containing edf files.
        xml_dir ->   (string) Path to directory containing xml files
        out_dir ->   (string) Path to output directory, where annotations files with events
                     will be saved.
        method ->    List of names of automated detection algorithms to detect 
                     events with. e.g. ['Lacourse2018','Moelle2011']
        chan ->      List of channel names in data (edf) file to detect on e.g. ['Cz']
        ref_chan ->  List of reference channels e.g. ['M1','M2'] (can be set to None)
        rater ->     (string) Rater name as listed in the annotations file e.g. 'Nathan' 
                     (can be set to None)
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
        stage ->     List of sleep stages to detect events in e.g. ['NREM2','NREM3']
        grp_name ->  (string) the channel group name as to be specified in the montage
                     (for later viewing)
        frequency -> Tuple of 2 digits, for the low and high frequency limits of any 
                     events to be detected.  
        duration ->  Tuple of 2 digits, for the minimum and maximum duration of any 
                     detected events.
        part ->      List of participants to perform detection on e.g. ['SP001']
                     Can also be set to 'all' (string) at which it will run on all
                     the participants it can find in rec_dir
        visit ->     List of visits to perform detection on e.g. ['V1']
                     Can also be set to 'all' (string) at which it will run on all
                     the visits it can find in rec_dir/<part>/
        
    
    '''
    
    
    # First we set up the output directory
    # a. Check for output folder, if doesn't exist, create
    if path.exists(out_dir):
            print(out_dir + " already exists")
    else:
        mkdir(out_dir)
    
    # loop through records
    if isinstance(part, list):
        None
    elif part == 'all':
            part = listdir(rec_dir)
            part = [ p for p in part if not '.' in p]
    else:
        print("ERROR: 'part' must either be an array of subject ids or = 'all' ")       
    
    print(r"""Whaling it... 
                          .
                       ":"
                     ___:____     |"\/"|
                   ,'        `.    \  /
                   |  O        \___/  |
                 ~^~^~^~^~^~^~^~^~^~^~^~^~
                 """)    

    for i, p in enumerate(part):
        print(p)
        # loop through visits
        if visit == 'all':
            visit = listdir(rec_dir + '/' + p)
            visit = [x for x in visit if not '.' in x]
            
        if i<1:
            detable=zeros((len(part),len(visit)*len(chan)*len(method)))
            header=[]
        
        for v, vis in enumerate(visit):
            ## Define files
            rdir = rec_dir + p + '/' + vis + '/'
            xdir = xml_dir + p + '/' + vis + '/'
            edf_file = [x for x in listdir(rdir) if x.endswith('.edf') or x.endswith('.rec') or x.endswith('.eeg')]
            xml_file = [x for x in listdir(xdir) if x.endswith('.xml')]            
            
            ## Copy annotations file before beginning
            if not path.exists(out_dir):
                mkdir(out_dir)
            if not path.exists(out_dir + p ):
                mkdir(out_dir + p)
            if not path.exists(out_dir + p + '/' + vis):
                mkdir(out_dir + p + '/' + vis)
            backup = out_dir + p + '/' + vis + '/'
            backup_file = (f'{backup}{p}_{vis}_spindles.xml')
            shutil.copy(xdir + xml_file[0], backup_file)
            
            ## Now import data
            dset = Dataset(rdir + edf_file[0])
            annot = Annotations(backup_file, rater_name=rater)
            
            # Get sleep cycles (if any)
            if cycle_idx is not None:
                all_cycles = annot.get_cycles()
                cycle = [all_cycles[y - 1] for y in cycle_idx if y <= len(all_cycles)]
            else:
                cycle = None
     
            ## Select and read data
            for c, ch in enumerate(chan):
                print(f'Reading data for {p}, visit {vis}, ' + str(ch))
                segments = fetch(dset, annot, cat=cat, stage=stage, cycle=cycle, reject_epoch=True, 
                                  reject_artf=['Artefact', 'Arou', 'Arousal'])
                segments.read_data([ch], ref_chan, grp_name=grp_name)
                
                ## Loop through methods (i.e. WHALE IT!)
    
                for m, meth in enumerate(method):
                    print(meth)
                    
                    ### define detection
                    
                    # Check for adapted bands
                    if adap_bands is True:
                        freq = frequency[ch][p + '_' + vis]
                        print(f'Using adapted bands for {p}: {round(freq[0],2)}-{round(freq[1],2)} Hz')
                    else:
                        freq = frequency
                    
                    
                    detection = DetectSpindle(meth, frequency=freq, duration=duration)
                        
                    ### run detection and save to Annotations file
                    all_spin = []
                    
                    if cat[0] == 1 and cat[1] == 0:
                        for s, seg in enumerate(segments):
                            print(f'Detecting events, stage {stage[s]}')
                            spindles = detection(seg['data'])
                            spindles.to_annot(annot, meth)
                            all_spin.append(spindles)
                            if len(spindles.events) == 0:
                                print(f'WARNING: No events detected by {meth} for {p}, {vis}')
                            detable[i,v*(len(chan)+len(method))+c+m-1] = len(spindles.events)
                        if i<1:
                            header.extend([f"{vis}_{ch}_{meth}_{stage[s]}"])
                    else:
                        for s, seg in enumerate(segments):
                            print('Detecting events, segment {} of {}'.format(s + 1, 
                                  len(segments)))
                            spindles = detection(seg['data'])
                            spindles.to_annot(annot, meth)
                            all_spin.append(spindles)
                            if len(spindles.events) == 0:
                                print(f'WARNING: No events detected by {meth} for {p}, {vis}')
                            detable[i,v*(len(chan)+len(method))+c+m-1] = len(spindles.events)
                        if i<1:
                            header.extend([f"{vis}_{ch}_{meth}_allstage"])
       
    ## all_spin contains some basic spindle characteistics
    print('Detection complete and saved.')   
    
    outable = DataFrame(detable)
    outable.columns=header  
    outable.index=part
    DataFrame.to_csv(outable, f'{out_dir}/detection_audit.csv', sep=',') 
    
    
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
                backup_file = (f'{backup_dir}/{pre}{keyword}.{ext}')
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


def whale_farm(rec_dir, xml_dir, out_dir, chan, grp_name, rater, stage=None, keyword = None,
               ref_chan=[], evt_name='spindle', segs=None, cycle_idx=None, adap_bands=False,
               frequency=(11,16),
               part='all', visit='all', param_keys=None, exclude_poor=False, 
               reject_artf=['Artefact', 'Arou', 'Arousal'], epoch_dur=30, n_fft_sec=4, 
               Ngo={'run':False}):
    
    
    if not path.exists(out_dir):
            mkdir(out_dir)

    logfile = out_dir + '/' + 'whales_log.txt'
    log = open(logfile, "a")
    
    try:
        print('')
        print(f'{datetime.now()}')
        print('Running whale_farm')  
        print('')
        with open(logfile, 'a') as f:
                        print('Running whale_farm', file=f) 
                        
        # Loop through subjects
        if isinstance(part, list):
            None
        elif part == 'all':
                part = listdir(rec_dir)
                part = [ p.split('.')[0] for p in part if '.set' in p if 'PSG' not in p]
        else:
            print("ERROR: 'part' must either be an array of subject ids or = 'all' **CHECK**")
            with open(logfile, 'a') as f:
                print('', file=f) 
                print("ERROR: 'part' must either be an array of subject ids or = 'all' ", file=f)  
        
        # Set base parameters
        params = {}
        sublist = []
        header = []
        
        if param_keys is None:
            param_keys = ['count', 'density','dur', 'ptp', 'energy', 'peakef'] # Default spindle metrics to extract
        if cycle_idx is not None:
                for m, param in enumerate(param_keys):
                    params[param] = zeros((len(part),len(chan)*len(cycle_idx)*len(stage)+len(chan)+len(stage)))
                    params[param].fill(nan)
        else:
                for m, param in enumerate(param_keys):
                    params[param] = zeros((len(part),len(chan)*len(stage)+len(chan)))
                    params[param].fill(nan)
        
        part.sort()                
        for i, p in enumerate(part):
            dat = []
            sublist.append(p)
            
            # Define files
            edf_file= [x for x in listdir(rec_dir) if p + '.set' in x]
            if len(edf_file)>1:
                edf_file=[x for x in edf_file if len(x)<6]
            
            # Loop through visits
            if visit == 'all':
                visiti = listdir(xml_dir + '/' + p)
                visiti = [x for x in visiti if not '.' in x]
            else:
                visiti = visit
            # Update size of file based on number of visits(sessions)
            if i == 0:
                for m, param in enumerate(param_keys):
                        params[param] = zeros((len(part), len(visiti)*len(params[param][0])))
                        params[param].fill(nan)
            
            visiti.sort()
            for v, vis in enumerate(visiti):  
                print(f'Extracting {evt_name} parameters for Subject {p}, visit {vis}..')
                with open(logfile, 'a') as f:
                    print('', file=f) 
                    print(f'Extracting {evt_name} parameters for Subject {p}, visit {vis}..', file=f)  
                

                
                if keyword is not None:
                    xml_file = [x for x in listdir(xml_dir + '/' + p + '/' + vis) 
                            if x.endswith('.xml') if not x.startswith('.') if keyword in x] 
                    
                else:
                    xml_file = [x for x in listdir(xml_dir +  '/' + p + '/' + vis ) 
                            if x.endswith('.xml') if not x.startswith('.')]

                    
                if len(xml_file) == 0:                
                    print(f'WARNING:{evt_name} has not been detected for Subject {p}, visit {vis} - skipping..')
                    with open(logfile, 'a') as f:
                        print('', file=f) 
                        print(f'WARNING: {evt_name} has not been detected for Subject {p}, visit {vis} - skipping..', file=f)  
                    
                else:
                    xml_file_path = xml_dir + '/' + p + '/' + vis + '/' + xml_file[0]
                    
                    # Open dataset
                    dataset = Dataset(rec_dir + edf_file[0])
            
                    # Import Annotations file
                    annot = Annotations(xml_file_path, 
                                        rater_name=rater)
                
                    # Get sleep cycles (if any)
                    if cycle_idx is not None:
                        all_cycles = annot.get_cycles()
                        cycle = [all_cycles[y - 1] for y in cycle_idx if y <= len(all_cycles)]
                    else:
                        cycle = None
                    
                    
                    # Run through channels
                    for ch, channel in enumerate(chan):
                        if Ngo['run'] == False:
                            if grp_name is not None:
                                chan_ful = [channel + ' (' + grp_name + ')']
                            else:
                                chan_ful = [channel]
                                print(f'channel is {chan_ful}')
                        else:
                            chan_ful = Ngo['chan']
                        
                        if adap_bands is True:
                            freq = frequency[channel][p + '_' + vis]
                            print(f'Using adapted bands for {p}: {round(freq[0],2)}-{round(freq[1],2)} Hz')
                        else:
                            freq = frequency
                        print(f'freq={freq}')
                        
                        # Create header for output file
                        for m, param in enumerate(param_keys):
                            if i == 0:
                                header.append(param + '_' + channel + '_wholenight' + '_visit_' + vis)
                                
                        if segs is not None:
                            for s, seg in enumerate(segs):
                                for m, param in enumerate(param_keys):
                                                if i == 0:
                                                    header.append(param + '_' + channel + '_' + seg[0] + '_visit_' + vis)  
                        else:
                            for s, st in enumerate(stage):
                                for m, param in enumerate(param_keys):
                                                if i == 0:
                                                    header.append(param + '_' + channel + '_' + st + '_visit_' + vis)
                                if cycle_idx is not None:
                                    for cy, cyc in enumerate(cycle_idx):
                                        for m, param in enumerate(param_keys):
                                                if i == 0:
                                                    header.append(param + '_' + channel + '_' + st + '_cycle' + str(cy+1) + '_visit_' + vis)
                    
                        
                        ### WHOLE NIGHT ###
                        # Select and read data
                        print('Reading data for ' + p + ', visit ' + vis + ' '+ str(channel))
                        with open(logfile, 'a') as f:
                            print('', file=f) 
                            print('Reading data for ' + p + ', visit ' + vis + ' ' + str(channel), file=f)
                        segments = fetch(dataset, annot, cat=(0,0,0,0), evt_type=[evt_name], cycle=cycle, 
                                         chan_full=chan_ful, reject_epoch=True, 
                                         reject_artf = reject_artf)
                        segments.read_data([channel], ref_chan, grp_name=grp_name)
                        if type(chan_ful)!=str and len(chan_ful)>1:
                            chan_ful = chan_ful[0]
                        
                        if len(segments) == 0:
                            print('')
                            print(f"WARNING: Events haven't been detected for {p}, {vis} on channel {channel}, skipping...")
                            print('')
                            with open(logfile, 'a') as f:
                                print('', file=f) 
                                print("WARNING: Events haven't been detected for {p}, {vis}  on channel {channel}, skipping...", file=f)
                        
                        # Calculate event density (whole night)
                        poi = get_times(annot, stage=stage, cycle=cycle, chan=[channel], exclude=exclude_poor)
                        total_dur = sum([x[1] - x[0] for y in poi for x in y['times']])
                        evts = annot.get_events(name=evt_name, time=None, chan = chan_ful, stage = stage)
                        count = len(evts)
                        density = len(evts) / (total_dur / epoch_dur)
                        print('')
                        print('----- WHOLE NIGHT -----')
                        print(f'No. Segments = {len(segments)}, Total duration (s) = {total_dur}')
                        print(f'Density = {density} per epoch')
                        print('')
                        with open(logfile, 'a') as f:
                            print('', file=f) 
                            print('----- WHOLE NIGHT -----', file=f) 
                            print(f'No. Segments = {len(segments)}, Total duration (s) = {total_dur}', file=f) 
                            print(f'Density = {density} per epoch', file=f) 
                        dat.append(count)
                        dat.append(density)
                        

                        
                        # Set n_fft
                        n_fft = None
                        if segments and n_fft_sec is not None:
                            s_freq = segments[0]['data'].s_freq
                            n_fft = int(n_fft_sec * s_freq)
                        
                        # Export event parameters (whole night)
                        data = event_params(segments, params='all', band=freq, n_fft=n_fft)
                        if data:
                                if not path.exists(out_dir + '/' + p):
                                    mkdir(out_dir + '/' + p)
                                if not path.exists(out_dir + '/' + p + '/' + vis):
                                    mkdir(out_dir + '/' + p + '/' + vis)
                                out_part = out_dir + '/' + p + '/' + vis
                                data = sorted(data, key=lambda x: x['start'])
                                outputfile = out_part + '/' + p + '_' + vis + '_' + channel + '_' + evt_name + '.csv'
                                print('Writing to ' + outputfile)
                                with open(logfile, 'a') as f:
                                    print('', file=f) 
                                    print('Writing to ' + outputfile, file=f) 
                                export_event_params(outputfile, data, count=len(evts), 
                                                    density=density)
                        else:
                            print('No valid data found.')
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
                                print('')
                                print(f'----- Segment {seg} -----')
                                print(f'No. Segments = {len(segments)}, Total duration (s) = {total_dur}')
                                print(f'Density = {density} per epoch')
                                print('')
                                with open(logfile, 'a') as f:
                                    print('', file=f) 
                                    print(f'----- Segment {seg} -----', file=f) 
                                    print(f'No. Segments = {len(segments)}, Total duration (s) = {total_dur}', file=f) 
                                    print(f'Density = {density} per epoch', file=f) 
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
                                        outputfile = out_dir + '/' + p + '_' + vis + '_' + channel + '_' + evt_name + '.csv'
                                        print('Writing to ' + outputfile)
                                        with open(logfile, 'a') as f:
                                            print('', file=f) 
                                            print('Writing to ' + outputfile, file=f) 
                                        export_event_params(outputfile, data, count=len(evts), 
                                                            density=density)
                                else:
                                    print('No valid data found.')
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
                            for s, st in enumerate(stage):
                                
                                if not isinstance(chan_ful, ndarray):
                                    chan_ful = [chan_ful]
                                
                                segments = fetch(dataset, annot, cat=(0,0,0,0), evt_type=[evt_name], 
                                                 stage = [st], cycle=cycle, 
                                                 chan_full=chan_ful, reject_epoch=True, 
                                                 reject_artf = reject_artf, min_dur=0.5)
                                segments.read_data([channel], ref_chan, grp_name=grp_name)
                                if len(chan_ful) > 1:
                                    chan_ful = chan_ful[0]
            
                                # Calculate event density (per stage)
                                poi = get_times(annot, stage=[st], cycle=cycle, chan=[channel], exclude=exclude_poor)
                                total_dur = sum([x[1] - x[0] for y in poi for x in y['times']])
                                evts = annot.get_events(name=evt_name, time=None, chan = chan_ful, stage = st)
                                count = len(evts)
                                density = len(evts) / (total_dur / epoch_dur)
                                print('')
                                print(f'---- STAGE {st} ----')
                                print(f'No. Segments = {len(segments)}, Total duration (s) = {total_dur}')
                                print(f'Density = {density} per epoch')
                                print('')
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
            
                                            
                                            segments.read_data([channel], ref_chan, grp_name=grp_name)
                                            
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
                                            print('')
                                            print(f'---- STAGE {st}, CYCLE {cy+1} ----')
                                            print(f'No. Segments = {len(segments)}, Total duration (s) = {total_dur}')
                                            print(f'Density = {density} per epoch')
                                            print('')
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
                                        print(f'No STAGE {st} in CYCLE {cy+1}')
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
        DataFrame.to_csv(output, f'{out_dir}/{evt_name}_dataset.csv', sep=',')            
        with open(f'{out_dir}/{evt_name}.p', 'wb') as f:
            dump(params, f)
            
        return output
        
    except Exception as e:
            print(e)
            tb = sys.exc_info()[2]
            tbinfo = traceback.format_tb(tb)[0]
            print(tbinfo)
            with open(logfile, 'a') as log:
                print(tbinfo, file=log)
                print(e, file=log)
                
def whale_farm_from_csvs(rec_dir, out_dir, chan, evt_name, part='all', visit='visit'):

    if not path.exists(out_dir):
            mkdir(out_dir)
    
    # Loop through subjects
    if isinstance(part, list):
        None
    elif part == 'all':
            part = listdir(rec_dir)
            part = [p for p in part if not '.' in p]
    else:
        print("ERROR: 'part' must either be an array of subject ids or = 'all' **CHECK**")
        
    # Loop through visits
    part.sort()                
    for i, p in enumerate(part):
        if visit == 'all':
            visit = listdir(rec_dir + '/' + p)
            visit = [x for x in visit if not '.' in x]    
    
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
        visit.sort()
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
    dataset.to_csv(out_dir + f'{evt_name}_{chan}_fromcsvs.csv')
                