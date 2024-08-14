#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 11:33:57 2022

@author: nathancross hello
"""

from math import floor
from numpy import where, zeros
from os import listdir, mkdir, path
from pandas import DataFrame
from wonambi import Dataset
from wonambi.attr import Annotations
from wonambi.trans import (fetch, math)


def perievent_histogram(rec_dir, xml_dir, out_dir, part, visit, stage, chan, ref_chan, cat,
                        grp_name, rater, standard, detection, chan_specific_markers,
                        window_size):

    # Make output directory
    if path.exists(out_dir):
        print(out_dir + " already exists")
    else:
        mkdir(out_dir)
    
    # Loop through channels
    for ch, channel in enumerate(chan):
        chan_ful = [channel + ' (' + grp_name + ')']
        
        # Loop through subjects and visits
        if isinstance(part, list):
            None
        elif part == 'all':
                part = listdir(xml_dir)
                part = [ p for p in part if not '.' in p]
        else:
            print("ERROR: 'part' must either be an array of subject ids or = 'all' **CHECK**")
            
        # Loop through visits
        visname = []
        if visit == 'all':
            visname = 'all'
            visit = listdir(xml_dir + '/' + part[0])
            visit = [x for x in visit if not '.' in x]
        else:
            visname.append([x for x in visit])
            
        # Make output array
        counts = zeros((len(part)*len(visit),int(window_size/0.1)))
        subvis = []
        index = -1
        
        part.sort()              
        for pp, p in enumerate(part):
            
            print(f'Running Subject {p}..')
            
            visit.sort()
            for v, vis in enumerate(visit): 
                
                # Define files
                rdir = rec_dir + '/' + p + '/' + vis + '/'
                #edf_file = [x for x in listdir(rdir) if x.endswith('.edf') if not x.startswith('.')]
                edf_file = [x for x in listdir(rdir) if x.endswith('.edf') or x.endswith('.rec') or x.endswith('.eeg')]
                xdir = xml_dir + '/' + p + '/' + vis + '/'
                xml_file = [x for x in listdir(xdir) if x.endswith('.xml') if not x.startswith('.')] 
                
                if len(xml_file) == 0:                
                    print(f'WARNING: no Annotations file found for Subject {p}, skipping..')
                else:
                    
                    # Import Annotations file
                    annot = Annotations(xdir + xml_file[0], rater_name=rater)
                    dataset = Dataset(rdir + edf_file[0])
                    s_freq = dataset.header['s_freq']
                                    
                    # Select and read data
                    print('Reading data for ' + p + ', visit ' + vis + ', channel '+ channel)
                    
                    # Retrieve standard events (to be compared against)
                    if chan_specific_markers:
                        stand_evts = fetch(dataset, annot, cat=cat, evt_type=standard, 
                                          stage=stage, cycle=None, chan_full=chan_ful, 
                                          reject_epoch=True)
                    else:
                        stand_evts = fetch(dataset, annot, cat=cat, evt_type=standard, 
                                          stage=stage, cycle=None,  reject_epoch=True)
                    
                    
                    starts = [x['times'][0][0] for x in stand_evts.segments]
                    
                    # Retrieve detection events (events to compare against standard)
                    if chan_specific_markers:
                        segments = fetch(dataset, annot, cat=cat, evt_type=detection, 
                                         stage=stage, cycle=None, chan_full=chan_ful, 
                                         reject_epoch=True)
                    else:
                        segments = fetch(dataset, annot, cat=cat, evt_type=detection, 
                                         stage=stage, cycle=None,  reject_epoch=True)
                    
                    segments.read_data([channel], ref_chan, grp_name=grp_name)
                    
                    
                    # Find the peaks of spindle events (relative to the start of the segment)
                    pktime_rel = zeros((len(segments)))
                    for i, seg in enumerate(segments):
                        dat = seg['data'] 
                        pk = math(dat, operator=_amax, axis='time').data[0][0]
                        pktime = where(seg['data'].data[0][0] == pk)[0][0]/s_freq
                        pktime_rel[i] = seg['data'].axis['time'][0][0] + pktime
                        
                    
                    histmat = zeros((len(pktime_rel),int(window_size/0.1)))
                    for y,pk in enumerate(pktime_rel):
                        for x,st in enumerate(starts):
                            
                            halfwin = window_size/2
                            # X seconds before
                            if -halfwin < pk-st < 0:
                                ind = floor((pk - st + halfwin)/0.1)
                                histmat[y,ind] = 1
                            # X seconds after
                            if 0 <= pk - st < halfwin:
                                ind = 40 + floor((pk - st)/0.1)
                                histmat[y,ind] = 1

                    index +=1
                    counts[index,:] = sum(histmat,axis=0)
                    subvis.append(p+'_'+vis)
        
        print(f"Saving histogram data to '{out_dir}peth_{channel}.csv'")
        output = DataFrame(counts, index = subvis)
        output.to_csv(f'{out_dir}peth_visit_{visname}_{channel}.csv')
                    
                        
def _amax(x, axis, keepdims=None):
    return amax(x, axis=axis)
