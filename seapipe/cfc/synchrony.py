# -*- coding: utf-8 -*-
"""


Author: Jordan O'Byrne & formatted by Nathan Cross
"""

from numpy import nan, ones
from os import listdir, mkdir, path
from pandas import DataFrame
from shutil import copy
from wonambi import Dataset
from wonambi.attr import Annotations
from wonambi.detect import match_events
from wonambi.trans import fetch



def event_sync(rec_dir, xml_dir, out_dir, part, visit, cat, 
               evttype_target, evttype_probe, iu_thresh, evttype_tp_target, 
               evttype_fn, chan=None, stage=None, grp=None, rater=None):
    
    
    '''
    This function divides a target event type from XML into events with and without a probe
    '''
    
    
    if path.exists(out_dir):
            print(out_dir + " already exists")
    else:
        mkdir(out_dir)
    
    # loop through records
    if isinstance(part, list):
        None
    elif part == 'all':
            part = listdir(xml_dir)
            part = [ p for p in part if not '.' in p]
    else:
        print("ERROR: 'part' must either be an array of subject ids or = 'all' ")       
    
    if visit == 'all': 
        visit = [x for x in listdir(xml_dir + '/' + part[0]) if not '.' in x]
    
    part.sort()
    visit.sort()
    
    print(f"""Merging events... "{evttype_target}" - with - "{evttype_probe}" - into "{evttype_tp_target}"
                    ___
                 , | l | 
                (( | l | ))
                   | l | '
                    \_/
                   /...\--.      _  
                   =====  `-----(_=
                 """)    
    
    # Create output 
    stats_header = ['Recall', 'Precision', 'F1 score', evttype_tp_target, f'{evttype_probe}-', evttype_fn]
    stats = ones((len(part)*len(visit), len(stats_header))) * nan
    
    # Check for channel specification
    if not chan:
        chan = [None]
        
    if not stage:
        stage = [None]
    
    # Concatenation
    cat = list(cat)
    cat[2] = 0
    cat[3] = 0 # force event type non-concatenation (required for analysis)
    cat = tuple(cat)
    
    if cat[1] == 1 and len(stage)>1:
        stage = [stage]
    

    for s, stg in enumerate(stage):
        for ch, channel in enumerate(chan):
            ids = []
            print(f'Channel {channel}; Stage {stg}')
  
            if channel:
                channel = channel + ' (' + grp + ')'

            for i, p in enumerate(part):
                
                
                for v, vis in enumerate(visit):
                    
                    print(f'Subject: {p}, Visit: {vis}')
                    print(f'{channel}')
                    ids.append(f'{p}_{vis}')
                    
                    ## Define files
                    rdir = rec_dir + p + '/' + vis + '/'
                    xdir = xml_dir + p + '/' + vis + '/'
                    edf_file = [x for x in listdir(rdir) if x.endswith('.edf') 
                                or x.endswith('.rec') or x.endswith('.eeg')]
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
                        copy(xdir + xml_file[0], backup_file)
                    else:
                        backup = out_dir + p + '/' + vis + '/'
                        backup_file = (f'{backup}{p}_{vis}_spindles.xml')
                    ## Import data
                    dset = Dataset(rdir + edf_file[0])
                    annot = Annotations(backup_file, rater_name=rater)
                    
                    
                    # Get events
                    # target
                    segments = fetch(dset, annot, cat=cat, evt_type=[evttype_target], 
                                     stage=stg, chan_full=[channel], reject_epoch=True, 
                                     reject_artf=['Artefact', 'Arou', 'Arousal'])
                    evt_target_seg = segments.segments
                    print(f'#targets = {len(evt_target_seg)}')
                    for s,seg in enumerate(evt_target_seg):
                        evt_target_seg[s]['start'] = seg['times'][0][0]
                        evt_target_seg[s]['end'] = seg['times'][0][1]
                        evt_target_seg[s]['chan'] = [seg['chan']]
                        evt_target_seg[s]['quality'] = 'Good'
                    
                    # target
                    segments = fetch(dset, annot, cat=cat, evt_type=[evttype_probe], 
                                     stage=stg, chan_full=[channel], reject_epoch=True, 
                                     reject_artf=['Artefact', 'Arou', 'Arousal'])
                    evt_probe_seg = segments.segments
                    print(f'#probes = {len(evt_probe_seg)}')
                    for s,seg in enumerate(evt_probe_seg):
                        evt_probe_seg[s]['start'] = seg['times'][0][0]
                        evt_probe_seg[s]['end'] = seg['times'][0][1]
                        evt_probe_seg[s]['chan'] = [seg['chan']]
                        evt_probe_seg[s]['quality'] = 'Good'
                        
                    
                    # Assess co-occurrence
                    matched = match_events(evt_probe_seg, evt_target_seg, iu_thresh)
         
                    # Write true positive target events to xml
                    matched.to_annot(annot, 'tp_std', evttype_tp_target)
                    
                    if evttype_fn is not None:
                        matched.to_annot(annot, 'fn', evttype_fn)
        


def event_sync_dataset(rec_dir, xml_dir, out_dir, part, visit, cat, outfile_suffix, 
               evttype_target, evttype_probe, iu_thresh, evttype_tp_target, 
               evttype_fn, chan=None, stage=None, grp='eeg', rater=None):
    
    
    '''
        This function accepts 2 event types: {target} and {probe}
        The {target} is split into 2 depending on the presence  of the probe at
        the same time. 
    '''
    
    
    if path.exists(out_dir):
            print(out_dir + " already exists")
    else:
        mkdir(out_dir)
    
    # loop through records
    if isinstance(part, list):
        None
    elif part == 'all':
            part = listdir(xml_dir)
            part = [ p for p in part if not '.' in p]
    else:
        print("ERROR: 'part' must either be an array of subject ids or = 'all' ")       
    
    if visit == 'all': 
        visit = [x for x in listdir(xml_dir + '/' + part[0]) if not '.' in x]
    
    part.sort()
    visit.sort()
    
    print(f" Extracting events and creating dataset...")    
    
    # Create output 
    stats_header = ['Recall', 'Precision', 'F1 score', evttype_tp_target, f'{evttype_probe}-', evttype_fn]
    stats = ones((len(part)*len(visit), len(stats_header))) * nan
    
    # Check for channel specification
    if not chan:
        chan = [None]
        
    if not stage:
        stage = [None]
    
    # Concatenation
    cat = list(cat)
    cat[2] = 0
    cat[3] = 0 # force event type non-concatenation (required for analysis)
    cat = tuple(cat)
    
    if cat[1] == 1 and len(stage)>1:
        stage = [stage]
    

    for s, stg in enumerate(stage):
        for ch, channel in enumerate(chan):
            ids = []
            print(f'Channel {channel}; Stage {stg}')
            
            
                
            if channel:
                if not stg:
                    stats_file = out_dir + channel + '_' + outfile_suffix
                elif len(stg)<2:
                    stats_file = out_dir + channel + '_' + stg[0] + '_' + outfile_suffix
                else:
                    if type(stg) is not list:
                        stg = [stg]
                    stages = '_'.join(stg)
                    stats_file = out_dir + channel + '_' + stages + '_' + outfile_suffix
                channel = channel + ' (' + grp + ')'
            else:
                stats_file = out_dir + outfile_suffix 
    
            
            for i, p in enumerate(part):
                
                
                for v, vis in enumerate(visit):
                    
                    print(f'Subject: {p}, Visit: {vis}')
                    print(f'{channel}')
                    ids.append(f'{p}_{vis}')
                    
                    ## Define files
                    rdir = rec_dir + p + '/' + vis + '/'
                    xdir = xml_dir + p + '/' + vis + '/'
                    edf_file = [x for x in listdir(rdir) if x.endswith('.edf') 
                                or x.endswith('.rec') or x.endswith('.eeg')]
                    xml_file = [x for x in listdir(xdir) if x.endswith('.xml')]            
                    
                    ## Import data
                    dset = Dataset(rdir + edf_file[0])
                    annot = Annotations(xdir + xml_file[0], rater_name=rater)
                    
                    # Get events
                    # target
                    segments = fetch(dset, annot, cat=cat, evt_type=[evttype_target], 
                                     stage=stg, chan_full=[channel], reject_epoch=True, 
                                     reject_artf=['Artefact', 'Arou', 'Arousal'])
                    evt_target_seg = segments.segments
                    print(f'#targets = {len(evt_target_seg)}')
                    for s,seg in enumerate(evt_target_seg):
                        evt_target_seg[s]['start'] = seg['times'][0][0]
                        evt_target_seg[s]['end'] = seg['times'][0][1]
                        evt_target_seg[s]['chan'] = [seg['chan']]
                        evt_target_seg[s]['quality'] = 'Good'
                    
                    # target
                    segments = fetch(dset, annot, cat=cat, evt_type=[evttype_probe], 
                                     stage=stg, chan_full=[channel], reject_epoch=True, 
                                     reject_artf=['Artefact', 'Arou', 'Arousal'])
                    evt_probe_seg = segments.segments
                    print(f'#probes = {len(evt_probe_seg)}')
                    for s,seg in enumerate(evt_probe_seg):
                        evt_probe_seg[s]['start'] = seg['times'][0][0]
                        evt_probe_seg[s]['end'] = seg['times'][0][1]
                        evt_probe_seg[s]['chan'] = [seg['chan']]
                        evt_probe_seg[s]['quality'] = 'Good'
                        
                    # Assess co-occurrence
                    matched = match_events(evt_probe_seg, evt_target_seg, iu_thresh)
                    
                    # store stats in master table
                    stats[(i*len(visit))+v, 0] = matched.recall
                    stats[(i*len(visit))+v, 1] = matched.precision
                    stats[(i*len(visit))+v, 2] = matched.f1score
                    stats[(i*len(visit))+v, 3] = matched.n_tp
                    stats[(i*len(visit))+v, 4] = matched.n_fp
                    stats[(i*len(visit))+v, 5] = matched.n_fn
                
            # export stats table
            print(f'Saving {stats_file}')
            df = DataFrame(data=stats, index=ids, columns=stats_header)
            df.to_csv(stats_file)
