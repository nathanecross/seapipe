# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 11:00:11 2021

@author: labottdv
"""

from os import listdir, mkdir, path, stat
import shutil
from wonambi import Dataset 
from wonambi.attr import Annotations
from wonambi.detect import DetectSlowWave
from wonambi.trans import fetch


def slow_it(rec_dir, xml_dir, out_dir, method, chan, rater, cat, stage, ref_chan, grp_name, 
            event_name, cycle_idx=None, duration=(0.2, 2), part='all', visit='all', 
            average_channels = False, invert = False):
    
    # First we set up the output directory
    # a. check for derivatives folder, if doesn't exist, create
    if not path.exists(out_dir ):
        mkdir(out_dir)
    
    # loop through records
    if isinstance(part, list):
        None
    elif part == 'all':
            part = listdir(rec_dir)
            part = [ p for p in part if not '.' in p]
    else:
        print("ERROR: 'part' must either be an array of subject ids or = 'all' ")       
    
    print(r"""Detecting slow waves... 
                          .
``'-.,_,.-'``'-.,_,.='``'-.,_,.-'``'-.,_,.='``
                 """)    
    
    for i, p in enumerate(part):
        print(p)
        
        # loop through visits
        if visit == 'all':
            visit = listdir(rec_dir + '/' + p)
            visit = [x for x in visit if not '.' in x]
        
        for v, vis in enumerate(visit):
            ## Define files
            rdir = rec_dir + p + '/' + vis + '/'
            edf_file = [x for x in listdir(rdir) if x.endswith('.edf') or x.endswith('.rec') or x.endswith('.eeg')]
            xdir = xml_dir + p + '/' + vis + '/'
            xml_file = [x for x in listdir(xdir) if x.endswith('.xml')]            
            
            ## Copy annotations file before beginning

            if not path.exists(out_dir + p):
                mkdir(out_dir + p)
            if not path.exists(out_dir + p + '/' + vis ):
                mkdir(out_dir + p + '/' + vis)
            backup = out_dir + p + '/' + vis  + '/'
            backup_file = (f'{backup}{p}_{vis}_{event_name}.xml' )
            shutil.copy2(xdir + xml_file[0], backup_file)
            
            ## Now import data
            dset = Dataset(rdir + edf_file[0])
            annot = Annotations(backup_file, rater_name=rater)
            
            # Get sleep cycles (if any)
            if cycle_idx is not None:
                all_cycles = annot.get_cycles()
                cycle = [all_cycles[y - 1] for y in cycle_idx if y <= len(all_cycles)]
            else:
                cycle = None
            
            print(f'Reading data for {p}, visit {vis}')
            ## Select and read data
            segments = fetch(dset, annot, cat, stage=stage, cycle=cycle, 
                                 reject_epoch=True, reject_artf=['Artefact', 'Arou', 'Arousal'])
            segments.read_data(chan, ref_chan, grp_name=grp_name, 
                                       average_channels=average_channels)

            ## Loop through methods 
            for m, meth in enumerate(method):
                print(meth)
                ### define detection
                detection = DetectSlowWave(meth, duration=duration)
                detection.invert = invert

                ### run detection and save to Annotations file
                all_slow = []
                
                for i, seg in enumerate(segments):
                    print('Detecting events, segment {} of {}'.format(i + 1, 
                          len(segments)))
                    swaves = detection(seg['data'])
                    swaves.to_annot(annot, event_name)
                    all_slow.append(swaves)
                    
    ## all_spin contains some basic spindle characteistics
    print('Detection complete and saved.')        
    return