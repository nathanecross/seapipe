#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 29 18:02:41 2024

@author: ncro8394
"""
from itertools import chain
from os import listdir, mkdir, path
from numpy import char, reshape
from pandas import DataFrame, read_csv, read_excel
from wonambi import Dataset
from wonambi.attr import Annotations, create_empty_annotations
from .logs import create_logger

def read_tracking_sheet(filepath, logger):
        
        track_file = [x for x in listdir(filepath) if 'tracking' in x]
        
        if len(track_file) > 1:
            logger.error('>1 tracking file found.')
            logger.warning("Is the tracking file currently open?")
            return 'error'
        elif len(track_file) == 0:
            logger.error('No tracking file found.')
            return 'error'
        else:
            track_file = track_file[0]
            if '.tsv' in track_file:
                track = read_csv(f'{filepath}/{track_file}' , sep='\t')
            elif '.xls' in track_file:
                track = read_excel(f'{filepath}/{track_file}')
    
        return track

def select_input_dirs(outpath, xml_dir, evt_name=None, 
                      logger = create_logger('Select inputs')):
    if not xml_dir:
        if evt_name in ['spindle', 'Ferrarelli2007', 'Nir2011', 'Martin2013', 
                        'Moelle2011', 'Wamsley2012', 'Ray2015', 'Lacourse2018', 
                        'FASST', 'FASST2', 'Concordia','UCSD','spindle_adap', 
                        'Ferrarelli2007_adap', 'Nir2011_adap', 'Martin2013_adap', 
                        'Moelle2011_adap', 'Wamsley2012_adap', 'Ray2015_adap', 
                        'Lacourse2018_adap', 'FASST_adap', 'FASST2_adap', 
                        'Concordia_adap','UCSD_adap']:
            xml_dir = f'{outpath}/spindle'
            
        elif evt_name in ['Ngo2015','Staresina2015','Massimini2004','slowwave',
                          'slowosc','SO']:
            xml_dir = f'{outpath}/slowwave'
            
        elif evt_name in ['event_pac']:
            xml_dir = f'{outpath}/event_pac'
            
        elif evt_name in ['pac']:
            xml_dir = f'{outpath}/pac' 
            
        elif evt_name in ['staging', 'macro', None]:
            xml_dir = [x for x in listdir(outpath) if 'staging' in x]
            
            if len(xml_dir) > 1:
                xml_dir = [x for x in xml_dir if 'manual' in x]
                logger.warning(f">1 derivatives directories have 'staging' in the "
                               f"name. Will default to {xml_dir}")
                
            if len(xml_dir) == 1:
                xml_dir = f'{outpath}/{xml_dir[0]}'

            else:
                logger.error(" 'xml_dir' wasn't specified and it cannot be "
                                "determined from inside the derivatives directory. "
                                "Please specify manually.") 
                xml_dir = ''
        else:
            xml_dir = f'{outpath}/{evt_name}'
        
    return xml_dir

def select_output_dirs(outpath, out_dir, evt_name=None, 
                      logger = create_logger('Select outputs')):
            
    if not out_dir:
        if evt_name in ['spindle', 'Ferrarelli2007', 'Nir2011', 'Martin2013', 
                        'Moelle2011', 'Wamsley2012', 'Ray2015', 'Lacourse2018', 
                        'FASST', 'FASST2', 'Concordia','UCSD','spindle_adap', 
                        'Ferrarelli2007_adap', 'Nir2011_adap', 'Martin2013_adap', 
                        'Moelle2011_adap', 'Wamsley2012_adap', 'Ray2015_adap', 
                        'Lacourse2018_adap', 'FASST_adap', 'FASST2_adap', 
                        'Concordia_adap','UCSD_adap']:
            out_dir = f'{outpath}/spindle'
        elif evt_name in ['Ngo2015','Staresina2015','Massimini2004','slowwave',
                          'slowosc','SO']:
            out_dir = f'{outpath}/slowwave'
        elif evt_name in ['pac']:
            out_dir = f'{outpath}/pac'
            
        elif evt_name in ['staging', 'macro', None]:
            out_dir = [x for x in listdir(outpath) if 'staging' in x]
            
            if len(out_dir) > 1:
                out_dir = [x for x in out_dir if 'manual' in x]
                
            if len(out_dir) == 1:
                out_dir = f'{outpath}/{out_dir[0]}'
                logger.warning(f">1 derivatives directories have 'staging' in "
                               f"the name. Will default to {out_dir}")
                
            else:
                logger.error(" 'xml_dir' wasn't specified and it cannot be "
                                "determined from inside the derivatives directory. "
                                "Please specify manually.") 
                xml_dir = ''    
        else:
            out_dir = f'{outpath}/{evt_name}'
    
    if not path.exists(out_dir):
        mkdir(out_dir)
        
    return out_dir

def load_stages(in_dir, xml_dir, subs = 'all', sessions = 'all', filetype = '.edf',
                stage_key = None, logger = create_logger('Load stages')):
    
    '''
        Extracts stages from the BIDS formatted dataset, in which
        staging has been listed in a file *acq-PSGScoring_events.tsv, and
        saves the information in an annotations file
        
        Parameters
        ----------
                    in_dir   :  str
                                The path to the BIDS dataset containing EEG recordings
                                and staging *acq-PSGScoring_events.tsv files
                    xml_dir  :  str
                                The derivatives path to store the annotations (.xml) 
                                file under the <sub>/<ses> structure
                    subs     :  str or ist of str
                                The participant ids to run this function on. Can 
                                be set to 'all', and all participants will be 
                                formatted.
                    sessions :  str or list of str
                                The participant ids to run this function on. Can 
                                be set to 'all', and all participants will be 
                                formatted.   
                    filetype :  str
                                The extension of EEG recording files 
                                (default = '.edf')
                    stage_key : dict or NoneType
                                Key for staging names to be saved into annotations
                                file (default is set to be compatible with Wonambi)
                    logger   :  logger for logging
                    

        Returns
        -------
                    flag : 0 (no errors)
                           1+ (errors or warnings)
        
        
    '''
    
    flag = 0
    if not stage_key:
        stage_key = {0: 'Wake',
                     1: 'NREM1',
                     2: 'NREM2',
                     3: 'NREM3',
                     4: 'REM'}
    
    # Define subs    
    if subs == 'all':
        subs = [x for x in listdir(in_dir) if '.' not in x]
    elif not isinstance(subs, list):  
        logger.error(f"'subs' must be a list of subject ids or 'all', not {subs}.")
        return 
    
    # Loop through
    subs.sort()
    for s, sub in enumerate(subs):
        
        # Make subject output directory
        if not path.exists(f'{xml_dir}/{sub}'):
            mkdir(f'{xml_dir}/{sub}')
            
        # Define sessions
        subdir = f'{in_dir}/{sub}'
        if sessions == 'all':
            sessions = [x for x in listdir(subdir) if '.' not in x]
        elif not isinstance(sessions, list):  
            logger.error(f"'sessions' must be a list of session ids or 'all', not {sessions}.")
            return
        
        # Loop through
        for s, ses in enumerate(sessions):
            
            # Make session output directory
            if not path.exists(f'{xml_dir}/{sub}/{ses}'):
                mkdir(f'{xml_dir}/{sub}/{ses}')
                
            # Load BIDS stage event file
            datadir = f'{in_dir}/{sub}/{ses}/eeg'
            stagefile = [x for x in listdir(datadir) if 'acq-PSGScoring_events.tsv' in x][0]
            stagedf = read_csv(f'{datadir}/{stagefile}', sep ='\t' ) 
            
            # Create new annotations for staging
            xml_file = f'{xml_dir}/{sub}/{ses}/{sub}_{ses}_staging.xml'
            if not path.exists(xml_file):
                edf_file = [x for x in listdir(datadir) if x.endswith(filetype)][0]
                dset = Dataset(f'{datadir}/{edf_file}')
                create_empty_annotations(xml_file, dset)
                logger.debug(f'Creating annotations file for {sub}, {ses}')
            else:
                logger.warning(f"Annotations file exists for {sub}, {ses}, staging "
                               "will be overwritten.")
                flag += 1
            annot = Annotations(xml_file)
            
            # Save staging to annotations
            annot.add_rater('manual')
            for i in range(0, stagedf.shape[0]):
                epoch_beg = stagedf.loc[i, 'onset']
                one_stage = stage_key[stagedf.loc[i, 'staging']]
                annot.set_stage_for_epoch(epoch_beg, one_stage,
                                         attr='stage',
                                         save=False)
            annot.save()
    
    return flag

def check_chans(datapath, chan, ref_chan, logger):
    if chan is None:
        chan = read_tracking_sheet(f'{datapath}', logger)
        if not isinstance(chan, DataFrame) and chan == 'error':
            logger.error("Channels haven't been defined, and there was an error "
                         "reading the tracking file.")
            logger.info('')
            logger.info("Check documentation for how to set up channel data: "
                        "https://seapipe.readthedocs.io/en/latest/index.html")
            logger.info('-' * 10)
            return chan, ref_chan
        
    if ref_chan is None:
        ref_chan = read_tracking_sheet(f'{datapath}', logger)
        if not isinstance(ref_chan, DataFrame) and ref_chan == 'error':
            logger.warning("Reference channels haven't been defined, and there "
                           "was an error reading the tracking file.")
            logger.info('')
            logger.info("Check documentation for how to set up channel data: "
                        "https://seapipe.readthedocs.io/en/latest/index.html")
            logger.info('-' * 10)
            logger.warning('No re-referencing will be performed prior to analysis.')
            ref_chan = None
    
    if ref_chan is False:
        return chan
    else:
        return chan, ref_chan

def load_sessions(sub, ses, rec_dir = None, flag = 0, 
                  logger = create_logger("Load Sessions"), verbose = 2):
    
    ''' Function to load in sessions available for a subject from a tracking 
        sheet.
    
        Parameters
        ----------
                sub : str
                      The subject id to load channel information for
                ses : str
                      The session id to load channel information for. Can be 
                      set as 'all' if all sessions for that subject are 
                      desired.
                rec_dir : str or None (optional)
                      The path to the BIDS dataset. If ses = 'all', this MUST
                      be specified.
                flag : int
                       the status of error logging
                logger : Class of logging.Logger
                         To log output (of course)
                         see: https://docs.python.org/3/library/logging.html
                verbose : int
                        The verbosity of logging.
                        verbose = 0 (error only)
                        verbose = 1 (warning only)
                        verbose = 2 (debug mode)
        
        Returns
        -------
                flag : int
                       whether an error occured in loading channel names
                
                chanset : dict
                         Format 'chan_name':'ref_chan' specifying the 
                         channel and it's associated reference channels.
    '''
    
    if type(ses) == type(DataFrame()):
        if verbose==2:
            logger.debug("Reading channel names from tracking file")
        # Search participant
        sub_row = ses[ses['sub']==sub]
        if sub_row.size == 0:
            if verbose>0:
                logger.warning(f"Participant {sub} not found in column 'sub' in "
                               "tracking file.")
                flag+=1
                return flag, None
        ses = [x for x in sub_row['ses']]
    elif type(ses) == str and ses == 'all':
        if not rec_dir:
            logger.error("If loading sessions from tracking sheet with ses = 'all', "
                         "rec_dir MUST be specified.")
            ses = None
            flag +=1
        else:   
            ses = listdir(rec_dir + '/' + sub)
            ses = [x for x in ses if not '.' in x]
    elif not type(ses) == list:
        logger.error("'sessions' must be set to None,'all' or a list of sub ids. "
                     "For session setup options, refer to documentation:")
        logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
        flag+=1
        return flag, None
    
    return flag, ses


def load_channels(sub, ses, chan, ref_chan, flag = 0, 
                  logger = create_logger("Load Channels"), verbose = 2):
    
    ''' Function to load in channels and their references from a tracking sheet.
    
        Parameters
        ----------
                sub : str
                      The subject id to load channel information for
                ses : str
                      The session id to load channel information for
                chan : list or DataFrame
                      The channel names to either load (DataFrame) or 
                      associate with a reference (list). If a DataFrame, 
                      the columns 'sub', 'ses', and 'chanset' must exist.
                ref_chan : list or DataFrame
                      The ref channel names to either load (DataFrame) or 
                      associate with a reference (list). If a DataFrame, 
                      the columns 'sub', 'ses' and 'refset' must exist.
                flag : int
                       the status of error logging
                logger : Class of logging.Logger
                         To log output (of course)
                         see: https://docs.python.org/3/library/logging.html
                verbose : int
                        The verbosity of logging.
                        verbose = 0 (error only)
                        verbose = 1 (warning only)
                        verbose = 2 (debug mode)
        
        Returns
        -------
                flag : int
                       whether an error occured in loading channel names
                
                chanset : dict
                         Format 'chan_name':'ref_chan' specifying the 
                         channel and it's associated reference channels.
    '''
    
    # verbose = 0 (error only)
    # verbose = 1 (warning only)
    # verbose = 2 (debug)
    
    # 1. Extract relevant cells from tracking sheet
    if isinstance(chan, DataFrame):
        if verbose==2:
            logger.debug("Reading channel names from tracking file")
        # Search participant
        chans = chan[chan['sub']==sub]
        if chans.size == 0:
            if verbose>0:
                logger.warning(f"Participant {sub} not found in column 'sub' in"
                               f" tracking file for {sub}, {ses}.")
            flag+=1
            return flag, None
        # Search session
        chans = chans[chans['ses']==ses]
        if chans.size == 0:
            if verbose>0:
                logger.warning(f"Session {ses} not found in column 'ses' in "
                               f"tracking file for {sub}, {ses}.")
            flag+=1
            return flag, None
        # Extract channels
        chans = search_chans(chans)
        chans = chans.dropna(axis=1, how='all')
        if len(chans.columns) == 0:
            if verbose>0:
                logger.warning("No channel set found in tracking file for "
                               f"{sub}, {ses}, skipping...")
            flag+=1
            return flag, None
    else:
        chans = chan
    
    if isinstance(ref_chan, DataFrame):
        if verbose==2:
            logger.debug("Reading reference channel names from tracking file ")
        ref_chans = ref_chan[ref_chan['sub']==sub]
        if ref_chans.size == 0:
            if verbose>0:
                logger.warning("Participant not found in column 'sub' in "
                               f"tracking file for {sub}, {ses}.")
            flag+=1
            return flag, None
        ref_chans = ref_chans[ref_chans['ses']==ses]
        if ref_chans.size == 0:
            if verbose>0:
                logger.warning("Session not found in column 'ses' in tracking "
                               f"file for {sub}, {ses}.")
            flag+=1
            return flag, None
        ref_chans = ref_chans.filter(regex='refset')
        ref_chans = ref_chans.dropna(axis=1, how='all')
        if len(ref_chans.columns) == 0:
            if verbose>0:
                logger.warning("No reference channel set found in tracking "
                               f"file for {sub}, {ses}. Progressing without "
                               "re-referencing...")
            ref_chans = []
    elif isinstance(ref_chan, list):
        ref_chans = ref_chan
    else:
        ref_chans = []
    
    if isinstance(chans, list):
        if isinstance(ref_chans, DataFrame) and len(ref_chans.columns) >1:
            chan = ref_chan[ref_chan['sub']==sub]
            chan = chan[chan['ses']==ses]
            chan = search_chans(chan)
            ref_chan=[]
            for c in chans:
                ref_chan.append([ref_chans[ref_chans.columns[x]].iloc[0] for 
                                 x, y in enumerate(chan) if c in chan[y].iloc[0]][0])
            ref_chan = [char.split(x, sep=', ').tolist() for x in ref_chan]
            chanset = {chn:[ref_chan[i]] if isinstance(ref_chan[i],str) else 
                       ref_chan[i] for i,chn in enumerate(chans)}
            
        elif isinstance(ref_chans, DataFrame):
            ref_chans = ref_chans.to_numpy()[0]
            ref_chans = ref_chans.astype(str)
            ref_chans_all = []
            for cell in ref_chans:
                cell = cell.split(', ')
                for x in cell:
                    if ',' in x: 
                        x = x.split(',')
                        ref_chans_all = ref_chans_all + x
                    else:
                        ref_chans_all.append(x)
            ref_chans = [x for x in ref_chans_all if not x=='']
            chanset = {chn:ref_chans for chn in chans} 
        
        elif isinstance(ref_chans, list):
            chanset = {chn:ref_chans for chn in chans}
        else:
            chanset = {chn:[] for chn in chans}
    
    elif isinstance(chans, DataFrame):
        if type(ref_chans) == DataFrame and len(ref_chans.columns) != len(chans.columns):
            logger.error("There must be the same number of channel sets and "
                         "reference channel sets in 'tracking file, but for "
                         f"{sub}, {ses}, there were {len(chans.columns)} channel "
                         f"sets and {len(ref_chans.columns)} reference channel sets. "
                         "For channel setup options, refer to documentation:")
            logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
            flag+=1
            return flag, None
        elif isinstance(ref_chans, DataFrame):
            ref_chans = ref_chans.to_numpy()[0]
            ref_chans = ref_chans.astype(str)
            ref_chans_all = []
            for cell in ref_chans:
                cell = cell.split(', ')
                refcell = []
                for x in cell:
                    if ',' in x: 
                        x = x.split(',')
                        refcell = refcell + x
                    else:
                        refcell.append(x)
                refcell = [x for x in refcell if not x=='']
                ref_chans_all.append(refcell)
            ref_chans = [x for x in ref_chans_all if not x=='']
        
        chans = chans.to_numpy()[0]
        chans = chans.astype(str)
        chans_all = []
        for cell in chans:
            cell = cell.split(', ')
            chancell = []
            for x in cell:
                if ',' in x: 
                    x = x.split(',')
                    chancell = chancell + x
                else:
                    chancell.append(x)
            chancell = [x for x in chancell if not x=='']
            chans_all.append(chancell)
        chans = chans_all

        if len(ref_chans)>0:
            chanset = {key:ref_chans[i] for i,chn in enumerate(chans) for key in chn}
        else:
            chanset = {key:[] for i,chn in enumerate(chans) for key in chn}
        
    else:
        logger.error("The variable 'chan' should be a [list] or definied in the "
                     "'chanset' column of tracking file - NOT a string.")
        flag+=1
        return flag, None
    
    return  flag, chanset


def rename_channels(sub, ses, chan, logger):
    
    if type(chan) == type(DataFrame()):
        # Search participant
        chans = chan[chan['sub']==sub]
        if len(chans.columns) == 0:
            return None
        # Search session
        chans = chans[chans['ses']==ses]
        if len(chans.columns) == 0:
            return None
        # Search channels
        oldchans = search_chans(chans)
        chans = chans.filter(regex='^((?!invert).)*$')
        oldchans = oldchans.dropna(axis=1, how='all')
        if len(oldchans.columns) == 0:
            return None
        newchans = chans.filter(regex='rename')
        newchans = newchans.dropna(axis=1, how='all')
        if len(newchans.columns) == 0:
            return None
    else:
        return None  
    
    if isinstance(oldchans,DataFrame):
        if isinstance(newchans,DataFrame) and len(newchans.columns) != len(oldchans.columns):
            try:
                oldchans_to_be_renamed = oldchans[list({i for i in oldchans if any(i in j for j in newchans)})]
                oldchans_to_be_kept = oldchans[list({i for i in oldchans if not any(i in j for j in newchans)})]
            except:
                logger.warning(f"There must be the same number of channel sets "
                               "and channel rename sets in tracking file, but "
                               f"for {sub}, {ses}, there were {len(oldchans.columns)} "
                               f"channel sets and {len(newchans.columns)} channel "
                               "rename sets. For info on how to rename channels, "
                               "refer to documentation:")
                logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
                logger.warning(f"Using original channel names for {sub}, {ses}...")
                return None
        else:
            oldchans_to_be_renamed = oldchans
            oldchans_to_be_kept = oldchans
        
        # Split cells in tracking
        # OLD Channels to be renamed
        oldchans_to_be_renamed = oldchans_to_be_renamed.to_numpy() 
        oldchans_to_be_renamed = oldchans_to_be_renamed[0].astype(str)
        oldchans_all = []
        for cell in oldchans_to_be_renamed:
            cell = cell.split(', ')
            chancell = []
            for x in cell:
                if ',' in x: 
                    x = x.split(',')
                    chancell = chancell + x
                else:
                    chancell.append(x)
            chancell = [x for x in chancell if not x=='']
            oldchans_all.append(chancell)       
        oldchans_to_be_renamed = list(chain(*oldchans_all))
        
        # OLD Channels to be kept
        oldchans_to_be_kept = oldchans_to_be_kept.to_numpy() 
        oldchans_to_be_kept = oldchans_to_be_kept[0].astype(str)
        oldchans_all = []
        for cell in oldchans_to_be_kept:
            cell = cell.split(', ')
            chancell = []
            for x in cell:
                if ',' in x: 
                    x = x.split(',')
                    chancell = chancell + x
                else:
                    chancell.append(x)
            chancell = [x for x in chancell if not x=='']
            oldchans_all.append(chancell)       
        oldchans_to_be_kept = list(chain(*oldchans_all))
        
        if type(newchans) == DataFrame:
            newchans = newchans.to_numpy() 
            newchans = newchans[0].astype(str)
            newchans_all = []
            for cell in newchans:
                cell = cell.split(', ')
                chancell = []
                for x in cell:
                    if ',' in x: 
                        x = x.split(',')
                        chancell = chancell + x
                    else:
                        chancell.append(x)
                chancell = [x for x in chancell if not x=='']
                newchans_all.append(chancell)       
            newchans = list(chain(*newchans_all))
        
        if len(oldchans_to_be_renamed) == len(newchans):
            newchans = {chn:newchans[i] for i,chn in enumerate(oldchans_to_be_renamed)}
            n = {x:x for x in oldchans_to_be_kept}
            newchans = n | newchans
        else:
            logger.warning(f"There must be the same number of original channel "
                           f"names and new renamed channels in tracking file, "
                           f"but for {sub}, {ses}, there were {len(oldchans)} old "
                           f"channel and {len(newchans)} new channel names. For "
                           "info on how to rename channels, refer to documentation:")
            logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
            logger.warning(f"Using original channel names for {sub}, {ses}...")
            return None
    else:
        return None
    
    return newchans


def load_stagechan(sub, ses, chan, ref_chan, flag, logger, verbose=2):
    # verbose = 0 (error only)
    # verbose = 1 (warning only)
    # verbose = 2 (debug)
    
    if type(chan) == type(DataFrame()):
        if verbose==2:
            logger.debug("Reading channel names from tracking file")
        # Search participant
        chans = chan[chan['sub']==sub]
        if chans.size == 0:
            if verbose>0:
                logger.warning(f"Participant {sub} not found in column 'sub' in "
                               f"tracking file for {sub}, {ses}.")
            flag+=1
            return flag, None
        # Search session
        chans = chans[chans['ses']==ses]
        if chans.size == 0:
            if verbose>0:
                logger.warning(f"Session {ses} not found in column 'ses' in "
                               f"tracking file for {sub}, {ses}.")
            flag+=1
            return flag, None
        # Search channel
        chans = chans.filter(regex='stagechan')
        chans = chans.dropna(axis=1, how='all')
        if len(chans.columns) == 0:
            if verbose>0:
                logger.warning(f"No stagechan found in tracking file for {sub}, "
                               f"{ses}, skipping...")
            flag+=1
            return flag, None
    else:
        chans = chan
    
    if type(ref_chan) == type(DataFrame()):
        if verbose==2:
            logger.debug("Reading reference channel names from tracking file")
        ref_chans = ref_chan[ref_chan['sub']==sub]
        if ref_chans.size == 0:
            if verbose>0:
                logger.warning(f"Participant not found in column 'sub' in tracking "
                               f"file for {sub}, {ses}.")
            flag+=1
            return flag, None
        ref_chans = ref_chans[ref_chans['ses']==ses]
        if ref_chans.size == 0:
            if verbose>0:
                logger.warning(f"Session not found in column 'ses' in tracking "
                               f"file for {sub}, {ses}.")
            flag+=1
            return flag, None
        ref_chans = ref_chans.filter(regex='refset')
        ref_chans = ref_chans.dropna(axis=1, how='all')
        if len(ref_chans.columns) == 0:
            if verbose>0:
                logger.warning(f"No reference channel set found in tracking file "
                               f"for {sub}, {ses}. Progressing without re-referencing...")
            ref_chans = []
    elif ref_chan:
        ref_chans = ref_chan
    else:
        ref_chans = []
    
    if type(chans) == list:
        if type(ref_chans) == DataFrame and len(ref_chans.columns) >1:
            chan = ref_chan[ref_chan['sub']==sub]
            chan = chan[chan['ses']==ses]
            chan = search_chans(chan)
            ref_chan=[]
            for c in chans:
                ref_link = [ref_chans[ref_chans.columns[x]].iloc[0] for x, y 
                                           in enumerate(chan) if c in chan[y].iloc[0]]
                if len(ref_link) > 1:
                    ref_chan.append(ref_link[0])
            ref_chan = [char.split(x, sep=', ').tolist() for x in ref_chan]
            
            if len(ref_chan) > 0:
                chanset = {chn:[ref_chan[i]] if isinstance(ref_chan[i],str) 
                           else ref_chan[i] for i,chn in enumerate(chans)}
            else:
                chanset = {chn:ref_chan for chn in chans}
            
        elif type(ref_chans) == DataFrame:
            ref_chans = ref_chans.to_numpy()[0]
            ref_chans = ref_chans.astype(str)
            ref_chans = char.split(ref_chans, sep=', ')
            ref_chans = [x for y in ref_chans for x in y]
            chanset = {chn:ref_chans for chn in chans}    
        else:
            chanset = {chn:[] for chn in chans}
    
    elif type(chans) == type(DataFrame()):
        
        if type(ref_chans) == DataFrame:
            if len(ref_chans.columns) != len(chans.columns):
                logger.warning(f"There were >2 reference channel sets in 'tracking' "
                               f"file for {sub}, {ses}, we will just use the "
                               "first set for automatic staging.")
                ref_chans = ref_chans.iloc[:,0]
                
            ref_chans = ref_chans.to_numpy()[0]
            if not isinstance(ref_chans, str):
                ref_chans = ref_chans.astype(str)
            ref_chans = char.split(ref_chans, sep=', ')
            if ref_chans.size < 2:
                ref_chans = reshape(ref_chans, (1,1))
                ref_chans = [x for x in ref_chans[0][0]]
            else:
                ref_chans = [x for x in ref_chans]

        chans = chans.to_numpy()[0]
        chans = chans.astype(str)
        chans = char.split(chans, sep=', ')
        chans = [x for x in chans]
        if len(ref_chans)>0:
            chanset = {key:ref_chans[i] for i,chn in enumerate(chans) for key in chn}
        else:
            chanset = {key:[] for i,chn in enumerate(chans) for key in chn}
        
    else:
        logger.error("The variable 'chan' should be definied in the 'chanset' "
                     "column of tracking file or a [list] - NOT a string.")
        flag+=1
        return flag, None
    
    return  flag, chanset


def load_eog(sub, ses, chan, flag, logger, verbose=2):
    # verbose = 0 (error only)
    # verbose = 1 (warning only)
    # verbose = 2 (debug)
    
    if type(chan) == type(DataFrame()):
        if verbose==2:
            logger.debug("Reading eog channel names from tracking file")
        # Search participant
        chans = chan[chan['sub']==sub]
        if chans.size == 0:
            if verbose>0:
                logger.warning(f"Participant not found in column 'sub' in "
                               f"the tracking file for {sub}, {ses}.")
            flag+=1
            return flag, None
        # Search session
        chans = chans[chans['ses']==ses]
        if chans.size == 0:
            if verbose>0:
                logger.warning(f"Session not found in column 'ses' in the tracking "
                               f"file for {sub}, {ses}.")
            flag+=1
            return flag, None
        # Search channel
        chans = chans.filter(regex='eog')
        chans = chans.dropna(axis=1, how='all')
        if len(chans.columns) == 0:
            if verbose>0:
                logger.warning("No stagechan found in tracking file for "
                               f"{sub}, {ses}...")
            flag+=1
            return flag, []
        chans = [x for x in chans['eog']]
    else:
        chans = chan
    return flag, chans


def load_emg(sub, ses, chan, flag, logger, verbose=2):
    # verbose = 0 (error only)
    # verbose = 1 (warning only)
    # verbose = 2 (debug)
    
    if type(chan) == type(DataFrame()):
        if verbose==2:
            logger.debug("Reading emg channel names from tracking file")
        # Search participant
        chans = chan[chan['sub']==sub]
        if chans.size == 0:
            if verbose>0:
                logger.warning("Participant not found in column 'sub' in the "
                               f"tracking file for {sub}, {ses}.")
            flag+=1
            return flag, None
        # Search session
        chans = chans[chans['ses']==ses]
        if chans.size == 0:
            if verbose>0:
                logger.warning("Session not found in column 'ses' in the "
                               f"tracking file for {sub}, {ses}.")
            flag+=1
            return flag, None
        # Search channel
        chans = chans.filter(regex='emg')
        chans = chans.dropna(axis=1, how='all')
        if len(chans.columns) == 0:
            if verbose>0:
                logger.warning(f"No emg found in tracking file for {sub}, {ses}...")
            flag+=1
            return flag, []
        chans = [x for x in chans['emg']]
    else:
        chans = chan
    return flag, chans

def check_adap_bands(rootpath, subs, sessions, chan, logger):
    
    try:
        track = read_tracking_sheet(rootpath, logger)
    except:
        logger.error("Error reading tracking sheet. Check that it isn't open.")
        logger.info("For info how to use adap_bands = 'Manual' in detections, "
                    "refer to documentation:")
        logger.info(" https://seapipe.readthedocs.io/en/latest/index.html")
        logger.info('-' * 10)
        return 'error'
    
    chans = search_chans(track)
    chans = chans.dropna(axis=1, how='all')
    peaks = track.filter(regex='peaks')
    peaks = peaks.dropna(axis=1, how='all')
    
    if len(peaks.columns) == 0:
        logger.error("No spectral peaks have been provided in tracking file. "
                     "Peaks will need to be detected.")
        logger.info("Check documentation for how to use adap_bands = 'Manual' "
                    "in detections: https://seapipe.readthedocs.io/en/latest/index.html")
        logger.info('-' * 10)
        return 'error'
    elif len(peaks.columns) != len(chans.columns):
        logger.error("There must be the same number of channel sets and spectral "
                     "peaks sets in tracking file")
        logger.info("Check documentation for how to use adap_bands = 'Manual' "
                    "in detections: https://seapipe.readthedocs.io/en/latest/index.html")
        return 'error'
    
    sub = {}
    for c, col in enumerate(chans.columns):
        for r, row in enumerate(chans[col]):
            chs = reshape(char.split(str(row), sep=', '), (1,1))[0][0]
            pks = reshape(char.split(str(peaks.iloc[r,c]), sep=', '), (1,1))[0][0]  
            if len(chs) != len(pks) and 'nan' not in (pks):
                logger.warning(f"For {track['sub'][r]}, {track['ses'][r]} the "
                               f"number of channels provided ({len(chs)}) != the "
                               f"number of spectral peaks ({len(pks)}).")
                if not track['sub'][r] in sub.keys():
                    sub[track['sub'][r]] = [track['ses'][r]]
                else:
                    sub[track['sub'][r]].append(track['ses'][r])
            elif 'nan' in (pks) and 'nan' not in (chs):
                logger.warning(f"For {track['sub'][r]}, {track['ses'][r]} no "
                               "peaks have been provided.")
                if not track['sub'][r] in sub.keys():
                    sub[track['sub'][r]] = [track['ses'][r]]
                else:
                    sub[track['sub'][r]].append(track['ses'][r])
    
    if len(sub) == 0:
        flag = 'approved'
        sub = 'all'
    else:
        flag = 'review'

    return flag
    

def read_manual_peaks(rootpath, sub, ses, chan, adap_bw, logger):
    
    try:
        track = read_tracking_sheet(rootpath, logger)
    except:
        logger.error("Error reading tracking sheet. Check that it isn't open.")
        logger.info("For info how to use adap_bands = 'Manual' in detections, "
                    "refer to documentation:")
        logger.info(" https://seapipe.readthedocs.io/en/latest/index.html")
        logger.info('-' * 10)
        return 'error'

    track = track[track['sub']==sub]
    if len(track.columns) == 0:
        logger.warning(f"Participant not found in column 'sub' in tracking file "
                       "for {sub}, {ses}.")
        return None
    # Search session
    track = track[track['ses']==ses]
    if len(track.columns) == 0:
        logger.warning("Session not found in column 'ses' in tracking file for "
                       f"{sub}, {ses}.")
        return None
    
    # Search channel
    chans = search_chans(track)
    chans = chans.dropna(axis=1, how='all')
    peaks = track.filter(regex='peaks')
    peaks = peaks.dropna(axis=1, how='all')
    
    if len(peaks.columns) == 0:
        logger.warning(f"No spectral peaks found in tracking file for {sub}, {ses}.")
        return None

    chans = chans.to_numpy()[0]
    chans = chans.astype(str)
    chans_all = []
    for cell in chans:
        cell = cell.split(', ')
        for x in cell:
            if ',' in x: 
                x = x.split(',')
                chans_all = chans_all + x
            else:
                chans_all.append(x)
    chans = [x for x in chans_all if not x=='']
    
    
    peaks = peaks.to_numpy()[0]
    peaks = peaks.astype(str)
    peaks_all = []
    for cell in peaks:
        cell = cell.split(', ')
        for x in cell:
            if ',' in x: 
                x = x.split(',')
                peaks_all = peaks_all + x
            else:
                peaks_all.append(x)
    peaks = [float(x) for x in peaks_all if not x=='']
    
    try:
        freq = (peaks[chans.index(chan)] - adap_bw/2, 
            peaks[chans.index(chan)] + adap_bw/2)
    except:
        logger.warning('Inconsistent number of peaks and number of channels '
                       'listed in tracking sheet for {sub}, {ses}. Will use '
                       'Fixed frequency bands instead...')
        freq = None
    
    return freq


def load_adap_bands(tracking, sub, ses, ch, stage, band_limits, adap_bw, logger):
    
    logger.debug(f'Searching for spectral peaks for {sub}, {ses}, {ch}.')
    
    try:
        files = tracking[sub][ses][ch]
    except:
        logger.warning(f'No specparams export file found for {sub}, {ses}, {ch}.')
        return None
    
    files = [x for x in files if stage in x['Stage']]
    files = [x for x in files if band_limits in x['Bandwidth']]

    if len(files) == 0:
        logger.warning(f"No specparams export file found for {sub}, {ses}, {ch}," 
                       f"{stage}, {band_limits}.")
        return None
    elif len(files) > 1:
        logger.warning(f">1 specparams export files found for {sub}, {ses}, {ch}, "
                       f"{stage}, {band_limits} ?")
        return None
    else:
        file = files[0]['File']
    
    
    # Read file and extract peak
    df = read_csv(file)
    df = df.filter(regex='peak')
    df = df.dropna(axis=1, how='all')
    
    if len(df.columns) == 3:
        peak = df.filter(regex='CF').values[0][0]
    elif len(df.columns) == 0: 
        logger.warning(f"No peaks found in export file for {sub}, {ses}, {ch}, "
                       f"{stage}, {band_limits}.")
        return None
    else:
        BW = df.filter(regex='BW')
        maxcol = BW.idxmax(axis='columns')[0].split('_')[1]
        df = df.filter(regex=maxcol)
        peak = df.filter(regex='CF').values[0][0]
            
    freq = (peak - adap_bw/2, 
            peak + adap_bw/2)
    
    return freq
    

def read_inversion(sub, ses, invert, chan, logger):
    
    if type(invert) == type(DataFrame()):
        # Search participant
        chans = invert[invert['sub']==sub]
        if len(chans.columns) == 0:
            logger.warning(f"Participant not found in column 'sub' in tracking "
                           f"file for {sub}, {ses}.")
            return None
        # Search session
        chans = chans[chans['ses']==ses]
        if len(chans.columns) == 0:
            logger.warning(f"Session not found in column 'ses' in tracking file "
                           f"for {sub}, {ses}.")
            return None
        
        # Search channel
        inversion = chans.filter(regex='invert')
        inversion = inversion.dropna(axis=1, how='all')
        chans = search_chans(chans)
        chans = chans.dropna(axis=1, how='all')
        
        if len(inversion.columns) == 0:
            logger.warning(f"No inversion info found in tracking file for "
                           f"{sub}, {ses}.")
            return None

        chans = chans.to_numpy()[0]
        chans = chans.astype(str)
        chans = char.split(chans, sep=', ')
        chans = [x for y in chans for x in y]
        
        inversion = inversion.to_numpy()[0]
        inversion = inversion.astype(str)
        inversion = char.split(inversion, sep=', ')
        inversion = [x for y in inversion for x in y]
        
        if len(inversion) == len(chans):
            inversion = inversion[chans.index(chan)]
            return inversion
        else:
            logger.warning(f"Error reading inversion info for {sub}, {ses}, {chan} "
                           "- check documentation for how to provide information "
                           "for inversion:")
            logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
            return None
    

def infer_ref(sub, ses, chan, logger, verbose=0):
    # verbose = 0 (error only)
    # verbose = 1 (warning only)
    # verbose = 2 (debug)
    logger.debug(f"Attempting to infer reference channel for {sub}, {ses} from "
                 "tracking sheet...")
    if isinstance(chan, DataFrame):
        # Search participant
        chans = chan[chan['sub']==sub]
        if len(chans.columns) == 0:
            if verbose>1:
                logger.debug(f"'{sub}' not found in tracking sheet.")
            return None
        # Search session
        chans = chans[chans['ses']==ses]
        if len(chans.columns) == 0:
            if verbose>1:
                logger.debug(f"{ses} for {sub} not found in tracking sheet.")
            return None
        # Search chansets (names in edf)
        oldchans = search_chans(chans)
        oldchans = oldchans.dropna(axis=1, how='all')
        if len(oldchans.columns) == 0:
            if verbose>1:
                logger.debug(f"'Chanset' for {sub} {ses} empty or not found in "
                             "tracking sheet (reference cannot be inferred).")
            return None
        # Search chansets (new names for output)
        newchans = chans.filter(regex='rename')
        newchans = newchans.dropna(axis=1, how='all')
        if len(newchans.columns) == 0:
            if verbose>1:
                logger.debug(f"Column 'newchans' empty or not found in tracking "
                             "sheet, so cannot infer reference channels true name.")
            return None
    else:
        return None
    
    # Infer ref chan name from oldchans:newchans mapping
    if isinstance(oldchans, DataFrame):
        if isinstance(newchans, DataFrame) and len(newchans.columns) != len(oldchans.columns):
            if verbose>0:
                logger.warning(f"There must be the same number of channel sets "
                               "and channel rename sets in tracking file, but "
                               f"for {sub}, {ses}, there were {len(oldchans.columns)} "
                               f"channel sets and {len(newchans.columns)} channel rename sets.")
                if verbose > 1:
                    logger.debug("For info on how to rename channels, refer to documentation:")
                    logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
            return None
        
        oldchans = oldchans.to_numpy() 
        oldchans = oldchans.astype(str)
        oldchans = char.split(oldchans[0], sep=', ')
        oldchans = [x for y in oldchans for x in y]
        
        if isinstance(newchans, DataFrame):
            newchans = newchans.to_numpy()
            newchans = newchans.astype(str)
            newchans = char.split(newchans[0], sep=', ')
            newchans = [x for y in newchans for x in y]
        
        if len(oldchans) == len(newchans):
            ref_chan = [newchans[i] for i,chn in enumerate(oldchans) if chn == '_REF'][0]
            if len(ref_chan) < 1:
                if verbose>1:
                    logger.debug(f"No channels named '_REF' in tracking sheet, "
                                 "so cannot infer reference channel.")
                return None
        else:
            if verbose>0:
                logger.warning(f"There must be the same number of original channel "
                               "names and new renamed channels in tracking file, "
                               f"but for {sub}, {ses}, there were {len(oldchans)} "
                               f"old channel and {len(newchans)} new channel names.")
                if verbose > 1:
                    logger.debug("For info on how to rename channels, refer to documentation:")
                    logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
            return None
    else:
        return None
    
    return ref_chan
    
    
def search_chans(chans):

    chans = chans.filter(regex='chanset')
    chans = chans.filter(regex='^((?!rename).)*$')
    chans = chans.filter(regex='^((?!peaks).)*$')
    chans = chans.filter(regex='^((?!invert).)*$')    
    
    return chans
    
    

