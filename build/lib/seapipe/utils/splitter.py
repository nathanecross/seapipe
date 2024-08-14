#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 17:38:31 2022

@author: nathancross
"""

from os import listdir, mkdir, path
import shutil
from wonambi.attr import Annotations


def extract_grouped_markers(xml_dir, out_dir, part, visit, rater, chan, 
                            target, splitter, out_event, file_type ='edf'):
    
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
                xml_file = [x for x in listdir(xml_dir + '/' + p + '/' + vis) if x.endswith('.xml') if not x.startswith(".")] 
                
                if len(xml_file) == 0:
                    print(f'WARNING: annotations does not exist for Subject {p}, visit {j} - check this. Skipping..')
                elif len(xml_file) >1:
                    print(f'WARNING: multiple annotations files exist for Subject {p}, visit {j} - check this. Skipping..')
                else:
                    print(f'Creating virtual markers for Subject {p}, Visit {vis}')
                    
                    print(f'Grouping markers for subject {p}, visit {vis}')
                    
                    # Copy annotations file before beginning
                    xdir = xml_dir + '/' + p + '/' + vis + '/' 
                    odir = out_dir + p + '/' + vis + '/'
                    backup_file = (f'{odir}{xml_file[0]}')
                    shutil.copy(xdir + xml_file[0], backup_file)
                    
                    # Open recording and annotations files
                    annot = Annotations(backup_file, rater_name=rater)  
                    
                    # Get list of events for target and splitters
                    a = annot.get_events(target)
                    b = annot.get_events(splitter)
                    
                    # Generate list of splitter start times
                    splittimes = [x['start'] for x in b]
                    
                    # Check if splitter begins before target
                    if a[0]['start'] < splittimes[0]:
                        print('WARNING: No splitter found before first target')
                    
                    #
                    out = []
                    for e1,evt1 in enumerate(a):
                       if e1 >0:
                           splitter_trgt = [x for x in b if a[e1-1]['end'] < x['start'] 
                                            < evt1['start']]
                            
                           if len(splitter_trgt) > 0:
                               out.append(evt1)
                    
     
                    for item in out:
                        item['name'] = out_event
                    
                    annot.add_events(out)
                    
                    
                    
                    
                    