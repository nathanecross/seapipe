#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 13:36:12 2023

@author: nathancross
"""
from copy import deepcopy
from datetime import datetime
from os import listdir, mkdir, path, rename, walk
from numpy import array, delete, zeros
from pandas import DataFrame, read_csv
from wonambi import Dataset
from wonambi.attr import Annotations
from .logs import create_logger
from .load import load_channels, rename_channels


def check_dataset(rootpath, outfile = False, filetype = '.edf', tracking = False, 
                  logger = None):
    
    """Audits the directory specified by <in_dir> to check if the dataset is
    BIDS compatible, how many sessions, recordings (e.g. edfs) and annotations
    files there are per participant.
    You can specify  an optional output filename that will contain the printout.
    """
    
    
    if not logger:
        logger = create_logger("Audit")
    
    datapath = path.join(rootpath, 'DATA')
    if not path.exists(datapath):
        logger.error(f'PATH: {datapath} does not exist. Check documentation for how to arrange data:\nhttps://seapipe.readthedocs.io/en/latest/index.html\n')
        return DataFrame()
    
    logger.debug(f'Checking dataset in directory: {datapath}')
    
    # Extract participants inside <indir>
    if not tracking:
        part = [x for x in listdir(datapath) if path.isdir(path.join(datapath, x))]
        part.sort()
    else:
        logger.debug(f'Reading participant names from {tracking}')
        part = read_csv(f'{rootpath}/{tracking}' , sep='\t')
        part = part['sub'].drop_duplicates()
    part.sort()
    
    nsd = []
    nedf = []
    bids = []
    finalbids = 0
    filesize = 0
    
    for pt in part:
        
        real_files = [x for x in listdir(path.join(datapath, pt))if not x.startswith('.')]
        dirs = [x for x in real_files if path.isdir(path.join(datapath, pt, x))]
        files = [x for x in real_files if path.isfile(path.join(datapath, pt, x))]
        
        nsd.append(len(dirs))
        annots = 0
        edfs = 0
        
        if len(dirs) < 1:
            finalbids += 1
            if len(files) > 0:
                nedf.append(len([x for x in files if filetype in x]))
                logger.critical(f'{pt} has 0 sessions directories.\n')
            else:
                logger.critical(f'{pt} has no files\n')
        else:
            for dr in dirs:
                try:
                    eeg_dir = path.join(datapath, pt, dr, 'eeg')
                    real_files = [x for x in listdir(eeg_dir) if not x.startswith('.')]
                    dirs2 = [x for x in real_files if path.isdir(path.join(eeg_dir, x))]
                    files2 = [x for x in real_files if path.isfile(path.join(eeg_dir, x))]
                    if len(dirs2) < 1:
                        if len(files2) == 1:
                            annots += len([x for x in files2 if '.xml' in x])
                            edfs += len([x for x in files2 if filetype in x])
                            filesize += sum([path.getsize(path.join(eeg_dir, x)) for x in files2 if filetype in x])
                        elif len(files2) > 1:
                            finalbids += 1
                            logger.critical(f'BIDS incompatibility. >1 file found for {pt}, {dr}. There should only be 1 recording per session directory\n')
                        else:
                            logger.warning(f'{pt}, {dr} has no files\n')
                except:
                    finalbids += 1
                    logger.critical(f"BIDS incompatibility. No 'eeg' directory found for {pt}, {dr}\n")
            
            bids.append(all([len(dirs2) < 1 for dirs2 in dirs]))
            nedf.append(edfs)
    
    if len(set(nsd)) > 1:
        logger.warning('Not all participants have the same number of sessions\n')
    
    subdirs = DataFrame({'BIDS?': bids, '#sessions': nsd, '#recordings': nedf},
                        index=part)
    subdirs[''] = ['!!' if c1 != c2 
                   else '!!' if c1 == 0 
                   else '!!' if c2 == 0 
                   else '' 
                   for c1, c2 in zip(subdirs['#sessions'], subdirs['#recordings'])]
    
    if outfile:
        subdirs.to_csv(outfile)
    
    if finalbids == 0:
        logger.info('\n                      Summary:')
        logger.info(f"                      {sum(subdirs['#recordings'])} files, {filesize / (10**9):,.2f} GB")
        logger.info(f"                      Subjects: {subdirs.shape[0]}")
        logger.info(f"                      Sessions: {max(subdirs['#sessions'])}\n")
        logger.debug('The dataset appears compatible for SEAPIPE analysis.\n')
    else:
        logger.critical('The dataset DOES NOT appear compatible for SEAPIPE analysis.\n')
    
    return subdirs



def make_bids(in_dir, origin = 'SCN', filetype = '.edf'):
    
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
                
                # EDFs
                files = [x for x in listdir(dst) if filetype in x] 
                for f, file in enumerate(files):
                    src = f'{in_dir}/sub-{part}/ses-{ses}/{file}'
                    dst = f'{in_dir}/sub-{part}/ses-{ses}/eeg/sub-{part}_ses-{ses}_eeg{filetype}'
                    rename(src, dst)
                
                # XMLs
                odir = '/'.join(in_dir.split('/')[0:-1]) + '/OUT/'
                if not path.exists(odir):
                    mkdir(odir)
                odir = f'{odir}/staging/'
                if not path.exists(odir):
                    mkdir(odir)
                odir = f'{odir}/sub-{part}/'
                if not path.exists(odir):
                    mkdir(odir)
                odir = f'{odir}/ses-{ses}/'
                if not path.exists(odir):
                    mkdir(odir)
                
                dst = f'{in_dir}/sub-{part}/ses-{ses}/'
                files = [x for x in listdir(dst) if '.xml' in x]
                for f, file in enumerate(files):
                    src = f'{in_dir}/sub-{part}/ses-{ses}/{file}'
                    
                    if len(file.split('_')) > 1:
                        newfile = file.split('_')[0]
                    else:
                        newfile = file.split('.')[0]
                    
                    dst = f'{odir}/sub-{newfile}_ses-{ses}_eeg.xml'
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
                
                
def track_processing(self, step, subs, tracking, df, chan, stage, show=False, 
                     log=True):

    ## Set up logging
    lg = create_logger('Tracking')
    
    ## Ensure correct format of chan and stage
    if isinstance(chan, str):
        chan = [chan]
    if isinstance(stage, str):
        stage = [stage]
    
    ## Track sleep staging
    if 'staging' in step or 'stage' in step:
         stage_df = []
         stage_dict = {}
         spath = self.outpath + '/staging/'
         for sub in subs:
            try:
                stage_ses = next(walk(f'{spath}/{sub}'))[1]
                stage_dict[sub] = dict([(x,[]) if x in stage_ses else (x,'-') 
                                     for x in tracking['ses'][sub]])
                stage_df.append([x if x in stage_ses else '-' 
                              for x in tracking['ses'][sub]])
            except:
                stage_df.append(['-'])
        
         # Update tracking
         tracking['staging'] = stage_dict 
         df['staging'] = stage_df
         
         # Check for Artefact or Arousal events
         if list(map(list, list(set(map(tuple, stage_dict.values()))))) == [['-']]:
            lg.warning('Staging has NOT been run.')
         else:
            for sub in stage_dict.keys():
                stage_ses = [x for x in stage_dict[sub].keys()]
                for ses in stage_ses:
                    tracking['staging'][sub][ses] = ['stage']
                    try:
                        xml = [x for x in listdir(f'{spath}/{sub}/{ses}') if '.xml' in x]
                        if len(xml) == 0:
                            if log:
                                lg.warning(f'No staging found for {sub}, {ses}')
                        elif len(xml) > 2: 
                            if log:
                                lg.warning(f'>1 staging files found for {sub}, {ses} - only 1 staging file is allowed.')
                        else:
                            xml = xml[0]
                            annot = Annotations(f'{spath}/{sub}/{ses}/{xml}')
                            events = sorted(set([x['name'] for x in annot.get_events()]))
                            for event in events:
                                if event in ['Arou', 'Arousal', 'Artefact']:
                                    tracking['staging'][sub][ses].append(event)  
                    except:
                        if log:
                            lg.warning(f'No staging found for {sub}, {ses}')    
                            
         
    ## Track spindle detection                 
    if 'spindles' in step or 'spindle' in step:
        spin_dict = {}
        spath = self.outpath + '/spindle/'
        df['spindle'] = [['-']] * len(df)
        spin_df = deepcopy(df['spindle'])
        
        for sub in subs:
           try:
               stage_ses = next(walk(f'{spath}/{sub}'))[1]
               spin_dict[sub] = dict([(x,{}) if x in stage_ses else (x,'-') 
                                    for x in tracking['ses'][sub]])
               spin_df.loc[sub] = [x if x in stage_ses else '-' 
                             for x in tracking['ses'][sub]]
           except:
               spin_dict[sub] = {'-':'-'}


        # Update tracking
        tracking['spindle'] = spin_dict 
        
        # Check for events
        if list(map(list, list(set(map(tuple, spin_dict.values()))))) == [['-']]:
            if log:
                lg.debug('Spindle detection has NOT been run.')
        else:
            methods = ['Lacourse2018','Moelle2011','Ferrarelli2007','Nir2011',
                       'Wamsley2012','Martin2013','Ray2015','FASST','FASST2',
                       'UCSD','Concordia','Lacourse2018_adap','Moelle2011_adap',
                       'Ferrarelli2007_adap','Nir2011_adap','Wamsley2012_adap',
                       'Martin2013_adap','Ray2015_adap','FASST_adap','FASST2_adap',
                       'UCSD_adap','Concordia_adap']
            for sub in spin_dict.keys():
                for ses in spin_dict[sub]:
                    if not spin_dict[sub][ses] == '-': 
                        xml = [x for x in listdir(f'{spath}/{sub}/{ses}') if '.xml' in x]
                        if len(xml) == 0:
                            if log:
                                lg.warning(f'No spindle annotations found for {sub}, {ses}')
                        elif len(xml) > 2: 
                            if log:
                                lg.warning(f'>1 spindle annotation files found for {sub}, {ses}..')
                        else:
                            xml = xml[0]
                            annot = Annotations(f'{spath}/{sub}/{ses}/{xml}')
                            events = [x for x in annot.get_events() if x['name'] in methods]
                            chans = sorted(set([x['chan'][0] for x in events]))
                            if chan:
                                chans = [x for x in chans for y in chan if y in x]
                            if len(chans) == 0:
                                if log:
                                    lg.warning(f'Spindles have NOT been detected for {sub}, {ses}.')
                                spin_dict[sub][ses] = ('-')
                                spin_df.loc[sub] = list(map(lambda x: x.replace(ses,'-'),spin_df.loc[sub]))
                                break
                            else:
                                for chan in chans:
                                    tracking['spindle'][sub][ses][chan] = []
                                    methlist = sorted(set([x['name'] for x in events]))
                                    if len(methlist) > 0:
                                        for method in methlist:
                                            update = datetime.fromtimestamp(path.getmtime(f'{spath}/{sub}/{ses}/{xml}')).strftime("%m-%d-%Y, %H:%M:%S")
                                            tracking['spindle'][sub][ses][chan].append({'Method':method,
                                                                                 'Stage':'',      # FLAG FOR UPDATE
                                                                                 'Cycle':'',      # FLAG FOR UPDATE
                                                                                 'File':f'{spath}/{sub}/{ses}/{xml}',
                                                                                 'Updated':update}) 

        df['spindle'] = spin_df
    
    
    ## Track slow oscillation detection                 
    if 'slow wave' in step or 'slow oscillation' in step or 'so' in step:
        so_dict = {}
        spath = self.outpath + '/slow_oscillation/'
        df['slow_osc'] = [['-']] * len(df)
        so_df = deepcopy(df['slow_osc'])
        
        for sub in subs:
           try:
               stage_ses = next(walk(f'{spath}/{sub}'))[1]
               so_dict[sub] = dict([(x,{}) if x in stage_ses else (x,'-') 
                                    for x in tracking['ses'][sub]])
               so_df.loc[sub] = [x if x in stage_ses else '-' 
                             for x in tracking['ses'][sub]]
           except:
               so_dict[sub] = {'-':'-'}

        # Update tracking
        tracking['slow_osc'] = so_dict 
        
        # Check for events
        if list(map(list, list(set(map(tuple, so_dict.values()))))) == [['-']]:
            if log:
                lg.debug('Slow oscillation detection has NOT been run.')
        else:
            methods = ['Massimini2004','AASM/Massimini2004','Ngo2015','Staresina2015',]
            for sub in so_dict.keys():
                for ses in so_dict[sub]:
                    if not so_dict[sub][ses] == '-': 
                        xml = [x for x in listdir(f'{spath}/{sub}/{ses}') if '.xml' in x]
                        if len(xml) == 0:
                            if log:
                                lg.warning(f'No slow oscillation annotations found for {sub}, {ses}')
                        elif len(xml) > 2: 
                            if log:
                                lg.warning(f'>1 slow oscillation annotation files found for {sub}, {ses}..')
                        else:
                            xml = xml[0]
                            annot = Annotations(f'{spath}/{sub}/{ses}/{xml}')
                            events = [x for x in annot.get_events() if x['name'] in methods]
                            chans = sorted(set([x['chan'][0] for x in events]))
                            if chan:
                                chans = [x for x in chans for y in chan if y in x]
                            if len(chans) == 0:
                                if log:
                                    lg.warning(f'Slow oscillations have NOT been detected for {sub}, {ses}.')
                                so_dict[sub][ses] = ('-')
                                so_df.loc[sub] = list(map(lambda x: x.replace(ses,'-'),so_df.loc[sub]))
                                break
                            else:
                                for chan in chans:
                                    tracking['slow_osc'][sub][ses][chan] = []
                                    methlist = sorted(set([x['name'] for x in events]))
                                    if len(methlist) > 0:
                                        for method in methlist:
                                            update = datetime.fromtimestamp(path.getmtime(f'{spath}/{sub}/{ses}/{xml}')).strftime("%m-%d-%Y, %H:%M:%S")
                                            tracking['slow_osc'][sub][ses][chan].append({'Method':method,
                                                                                 'Stage':'',      # FLAG FOR UPDATE
                                                                                 'Cycle':'',      # FLAG FOR UPDATE
                                                                                 'File':f'{spath}/{sub}/{ses}/{xml}',
                                                                                 'Updated':update}) 

        df['slow_osc'] = so_df
    
    
    ## Track fooof detection                 
    if 'fooof' in step or 'specparams' in step:
        fooof_dict = {}
        spath = self.outpath + '/fooof/'
        df['fooof'] = [['-']] * len(df)
        fooof_df = deepcopy(df['fooof'])
        
        for sub in subs:
           try:
               stage_ses = next(walk(f'{spath}/{sub}'))[1]
               fooof_dict[sub] = dict([(x,{}) if x in stage_ses else (x,'-') 
                                    for x in tracking['ses'][sub]])
               fooof_df.loc[sub] = [x if x in stage_ses else '-' 
                             for x in tracking['ses'][sub]]
           except:
               fooof_dict[sub] = dict([(ses,'-') for ses in self.tracking['ses'][sub]])


        # Update tracking
        tracking['fooof'] = fooof_dict 
        
        # Check for events
        if list(map(list, list(set(map(tuple, fooof_dict.values()))))) == [['-']]:
            if log:
                lg.debug('FOOOF detection has NOT been run.')
        else:
            for sub in fooof_dict.keys():
                for ses in fooof_dict[sub]:
                    if not fooof_dict[sub][ses] == '-': 
                        files = [x for x in listdir(f'{spath}/{sub}/{ses}') if '.csv' in x]
                        
                        chans = sorted(set([file.split('_')[2] for file in files]))
                        if chan:
                            chans = [x for x in chans for y in chan if y in x]
                        if len(chans) == 0:
                            if log:
                                lg.warning(f'FOOOF has NOT been run for {sub}, {ses}.')
                            fooof_dict[sub][ses] = ('-')
                            fooof_df.loc[sub] = list(map(lambda x: x.replace(ses,'-'),fooof_df.loc[sub]))
                            break
                        else:
                            for chan in chans:
                                tracking['fooof'][sub][ses][chan] = []
                                chan_files = [file for file in files if chan in file]
                                for chanfile in chan_files:
                                    update = datetime.fromtimestamp(path.getmtime(f'{spath}/{sub}/{ses}/{chanfile}')).strftime("%m-%d-%Y, %H:%M:%S")
                                    tracking['fooof'][sub][ses][chan].append({'Stage':chanfile.split('_')[3],      
                                                                              'Cycle':'',      # FLAG FOR UPDATE
                                                                              'Bandwidth':chanfile.split('_')[-1].split('.csv')[0],
                                                                              'File':f'{spath}/{sub}/{ses}/{chanfile}',
                                                                              'Updated':update})
        df['fooof'] = fooof_df

    return df, tracking               
                
                
def check_fooof(self, frequency, chan, ref_chan, stage, cat, cycle_idx, logger):

    bandwidth = f'{frequency[0]}-{frequency[1]}Hz'
    review = []
    for sub in self.tracking['fooof']:
        sessions = list(self.tracking['fooof'][sub].keys())
        for ses in sessions:
            if not self.tracking['fooof'][sub][ses] == '-':
                flag, chanset = load_channels(sub, ses, chan, ref_chan, 0, 
                                              logger, verbose=0)
                if flag>0:
                    return 'error', None, None, None
                newchans = rename_channels(sub, ses, chan, logger)
                
                for c, ch in enumerate(chanset):
                    if newchans:
                        fnamechan = newchans[ch]
                    else:
                        fnamechan = ch
                    try:
                        fooof = self.tracking['fooof'][sub][ses][fnamechan]
                    except:
                        logger.warning(f'No fooof for {sub}, {ses}, {ch}')
                        break
                    if cat[0] + cat[1] == 2: # whole night
                        num_files = 1
                        stagename = '-'.join(stage)
                        files = [x['File'] for x in fooof if stagename in x['Stage'] 
                                 if bandwidth in x['Bandwidth']]
                    elif cat[0] + cat[1] == 0: # stage*cycle
                        # num_files = len(stage)*len(cycle_idx)
                        # files = []
                        # for stg in stage:
                        #     for cyc in cycle_idx:
                        #         files.append([x['File'] for x in fooof 
                        #                       if stage in x['Stage'] 
                        #                       if cyc in x['Cycle']
                        #                       if bandwidth in x['Bandwidth']])
                        logger.error('Adapted bands for stage*cycle has not yet been implemented')
                        return 'error', None, None, None
                    elif cat[0] == 0:
                        # num_files = len(cycle_idx)
                        # files = []
                        # for cyc in cycle_idx:
                        #     files.append([x['File'] for x in fooof 
                        #                   if stage in x['Stage'] 
                        #                   if cyc in x['Cycle']
                        #                   if bandwidth in x['Bandwidth']])
                        logger.error('Adapted bands for per_cycle has not yet been implemented')
                        return 'error', None, None, None
                    elif cat[1] == 0:
                        num_files = len(stage)
                        files = []
                        for stg in stage:
                            [files.append(x['File']) for x in fooof if stg in x['Stage']
                                          if bandwidth in x['Bandwidth']]
                    
                    if num_files != len(files):
                        flag +=1
            else:
                flag = 1
        
            if flag>0:
                review.append([sub, ses])
    
    for row in chan.index:
        subses = [chan['sub'].loc[row], chan['ses'].loc[row]]
        if not subses in review:
            chan = chan.drop([row])
    
    sub = list(chan['sub'])
    ses = list(chan['ses'])
    
    return 'review', chan, sub, ses
        
                
                
                
                
                