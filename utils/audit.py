#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 13:36:12 2023

@author: nathancross
"""
import logging
import sys
from os import listdir, mkdir, rename
from numpy import array, delete, zeros
from pandas import DataFrame
from wonambi import Dataset

def check_dataset(indir, filetype = '.edf', outfile = True):
    
    """Audits the directory specified by <in_dir> to check if the dataset is
    BIDS compatible, how many sessions, recordings (e.g. edfs) and annotations
    files there are per participant.
    You can specify  an optional output filename that will contain the printout.
    """
    
    # Set logging 
    logging.basicConfig(level=logging.DEBUG,
                        format="%(message)s",
                        handlers=[logging.StreamHandler(sys.stdout)])
    logger = logging.getLogger() 
    
    # Begin directory check
    logger.info('')
    logger.info(f'Checking dataset in directory: {indir}')
    logger.info('')
    
    # Extract participants inside <indir>
    part = [x for x in listdir(indir) if not '.' in x]
    part.sort()
    
    # Initiate arrays
    nsd = zeros(len(part), dtype=int)
    nedf = zeros(len(part), dtype=int)
    nxml = zeros(len(part), dtype=int)
    bids = zeros((len(part),1), dtype=bool)
    
    # Check number of subdirectories per participant
    for p, pt in enumerate(part):
        dirs = [x for x in listdir(f'{indir}/{pt}') if not '.' in x]
        files = [x for x in listdir(f'{indir}/{pt}') if '.' in x if not '.DS_Store' in x]
        nsd[p] = len(dirs)
        annots = 0
        edfs = 0
        
        if len(dirs) < 1:
            if len(files) > 0: 
                nxml[p] = len([x for x in files if '.xml' in x])
                nedf[p] = len([x for x in files if filetype in x])
                logger.info('')
                logger.info(f'WARNING: {pt} has 0 sessions directories')
                logger.info('')
            else:
                logger.info('')
                logger.info(f'WARNING: {pt} has no files')
                logger.info('')
                
        else:    
            for d, dr in enumerate(dirs):
                dirs2 = [x for x in listdir(f'{indir}/{pt}/{dr}/eeg/') if not '.' in x]
                files2 = [x for x in listdir(f'{indir}/{pt}/{dr}/eeg/') if '.' in x]
                dirscheck = []
                if len(dirs2) < 1:
                    if len(files2) > 0: 
                        annots += len([x for x in files2 if '.xml' in x])
                        edfs += len([x for x in files2 if filetype in x])
                        dirscheck.append(0)
                    else:
                        logger.info('')
                        logger.info(f'WARNING: {pt} has no files')
                        logger.info('')
                        dirscheck.append(1)
            if set(dirscheck) == {0}:
                bids[p] = 1
                nedf[p] = edfs
                nxml[p] = annots
                
    if len(set(nsd)) > 1:  
        logger.info('')
        logger.info('WARNING: Not all participants have the same number of sessions')
        logger.info('')
    
    # Create audit dataframe
    subdirs = DataFrame(zeros((len(part),4)), 
                        index=part,
                        columns=['BIDS?','#sessions','#recordings','#annotations'])
    subdirs.index.name = 'Participants'
    
    subdirs[subdirs.columns[0]] = bids
    subdirs[subdirs.columns[1]] = nsd
    subdirs[subdirs.columns[2]] = nedf
    subdirs[subdirs.columns[3]] = nxml
    
    if outfile:
        subdirs.to_csv(f'{indir}/{outfile}')

    return subdirs

def make_bids(in_dir, origin = 'SCN'):
    
    """Converts the directory specified by <in_dir> to be BIDS compatible.
    You can specify the origin format of the data. For now, this only converts
    from the Sleep Cognition Neuroimaging laboratory format, but please contact 
    me (nathan.cross.90@gmail.com)if you would like more formats.
    """
    
    if origin=='SCN':
        parts = [x for x in listdir(in_dir) if '.' not in x]
        
        for p, part in enumerate(parts):
            
            src = f'{in_dir}/{part}'
            dst = f'{in_dir}/sub-{part}'
            rename(src, dst)
            
            sess = [x for x in listdir(dst) if '.' not in x]
            
            for s, ses in enumerate(sess):
                src = f'{in_dir}/sub-{part}/{ses}'
                dst = f'{in_dir}/sub-{part}/ses-{ses}/'
                rename(src, dst)
                
                mkdir(f'{in_dir}/sub-{part}/ses-{ses}/eeg/')
                
                files = [x for x in listdir(dst) if '.edf' in x] 
                
                for f, file in enumerate(files):
                    src = f'{in_dir}/sub-{part}/ses-{ses}/{file}'
                    
                    newfile = file.split('_')[0]
                    dst = f'{in_dir}/sub-{part}/ses-{ses}/eeg/sub-{newfile}_ses-{ses}_eeg.edf'
                    rename(src, dst)
                    
def extract_channels(in_dir, exclude=['A1','A2','M1','M2'], quality=False):
    
    """Reads channel information from the files in the directory specified by 
    <in_dir> and writes them to the BIDS compatible channels.tsv file per participant
    and session.
    You can specify whether to exclude any channels, if preferrable.
    """
    
    parts = [x for x in listdir(in_dir) if '.' not in x]
    
    for p, part in enumerate(parts):
        ppath = f'{in_dir}/{part}'
        sess = [x for x in listdir(ppath) if '.' not in x]
        
        for s, ses in enumerate(sess):
            spath = f'{ppath}/{ses}/eeg/'
            files = [x for x in listdir(spath) if '.edf' in x] 
            
            for f, file in enumerate(files):
                src = f'{spath}/{file}'
                
                data = Dataset(src)
                chans = data.header['orig']['label'] #data.header['chan_name']
                types = array([x.split('-')[0] for x in data.header['orig']['transducer']])
                units = array(data.header['orig']['physical_dim'])
                
                if exclude:
                    ex = [chans.index(x) for x in exclude if x in chans]
                    chans = delete(array(chans), ex)
                    types = delete(types, ex)
                    units = delete(units, ex)
                else:
                    chans = array(chans)
                
                # Save dataframe
                df = DataFrame(chans)
                df.columns = ['name']
                df['type'] = types
                df['units'] = units
                df['status'] = 'N/A'
                df['status_description'] = 'N/A'
                df.to_csv(f"{spath}{part}_{ses}_channels.tsv", sep = "\t", 
                          header=True, index=False)
                
                
                
                
                
                
                
                
                
                
                
                