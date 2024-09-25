#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 29 18:02:41 2024

@author: ncro8394
"""
from os import listdir, mkdir, path
from numpy import char, reshape
from pandas import DataFrame, read_csv, read_excel


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

def select_input_dirs(outpath, xml_dir, evt_name=None):
    if not xml_dir:
        if evt_name in ['spindle', 'Ferrarelli2007', 'Nir2011', 'Martin2013', 
                        'Moelle2011', 'Wamsley2012', 'Ray2015', 'Lacourse2018', 
                        'FASST', 'FASST2', 'Concordia','UCSD','spindle_adap', 
                        'Ferrarelli2007_adap', 'Nir2011_adap', 'Martin2013_adap', 
                        'Moelle2011_adap', 'Wamsley2012_adap', 'Ray2015_adap', 'Lacourse2018_adap', 
                        'FASST_adap', 'FASST2_adap', 'Concordia_adap','UCSD_adap']:
            xml_dir = f'{outpath}/spindle'
        elif evt_name in ['Ngo2015','Staresina2015','Massimini2004','slowwave',
                          'slowosc','SO']:
            xml_dir = f'{outpath}/slowwave'
        elif evt_name in ['macro']:
            xml_dir = f'{outpath}/staging'
        elif evt_name is None:
            xml_dir = f'{outpath}/staging'
        else:
            xml_dir = f'{outpath}/{evt_name}'
        
    return xml_dir

def select_output_dirs(outpath, out_dir, evt_name=None):
            
    if not out_dir:
        if evt_name in ['spindle', 'Ferrarelli2007', 'Nir2011', 'Martin2013', 
                        'Moelle2011', 'Wamsley2012', 'Ray2015', 'Lacourse2018', 
                        'FASST', 'FASST2', 'Concordia','UCSD','spindle_adap', 
                        'Ferrarelli2007_adap', 'Nir2011_adap', 'Martin2013_adap', 
                        'Moelle2011_adap', 'Wamsley2012_adap', 'Ray2015_adap', 'Lacourse2018_adap', 
                        'FASST_adap', 'FASST2_adap', 'Concordia_adap','UCSD_adap']:
            out_dir = f'{outpath}/spindle'
        elif evt_name in ['Ngo2015','Staresina2015','Massimini2004','slowwave',
                          'slowosc','SO']:
            out_dir = f'{outpath}/slowwave'
        elif evt_name in ['macro']:
            out_dir = f'{outpath}/staging'
        elif evt_name in ['pac']:
            out_dir = f'{outpath}/pac'
        else:
            out_dir = f'{outpath}/{evt_name}'
    
    if not path.exists(out_dir):
        mkdir(out_dir)
        
    return out_dir

def check_chans(rootpath, chan, ref_chan, logger):
    if chan is None:
        chan = read_tracking_sheet(f'{rootpath}', logger)
        if not isinstance(chan, DataFrame) and chan == 'error':
            logger.error("Channels haven't been defined, and there was an error reading the tracking file.")
            logger.info('')
            logger.info('Check documentation for how to set up channel data: https://seapipe.readthedocs.io/en/latest/index.html')
            logger.info('-' * 10)
        
    if ref_chan is None:
        ref_chan = read_tracking_sheet(f'{rootpath}', logger)
        if not isinstance(ref_chan, DataFrame) and ref_chan == 'error':
            logger.warning("Reference channels haven't been defined, and there was an error reading the tracking file.")
            logger.info('')
            logger.info('Check documentation for how to set up channel data: https://seapipe.readthedocs.io/en/latest/index.html')
            logger.info('-' * 10)
            logger.warning('No re-referencing will be performed prior to analysis.')
            ref_chan = None
    
    if ref_chan is False:
        return chan
    else:
        return chan, ref_chan

def load_sessions(sub, ses, rec_dir, flag, logger, verbose=2):
    # verbose = 0 (error only)
    # verbose = 1 (warning only)
    # verbose = 2 (debug)
    
    if type(ses) == type(DataFrame()):
        if verbose==2:
            logger.debug("Reading channel names from tracking file")
        # Search participant
        sub_row = ses[ses['sub']==sub]
        if sub_row.size == 0:
            if verbose>0:
                logger.warning(f"Participant {sub} not found in column 'sub' in tracking file.")
                flag+=1
                return flag, None
        ses = [x for x in sub_row['ses']]
    elif type(ses) == str and ses == 'all':
            ses = listdir(rec_dir + '/' + sub)
            ses = [x for x in ses if not '.' in x]
    elif not type(ses) == list:
        logger.error("'sessions' must be set to None,'all' or a list of sub ids. For session setup options, refer to documentation:")
        logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
        flag+=1
        return flag, None
    
    return flag, ses


def load_channels(sub, ses, chan, ref_chan, flag, logger, verbose=2):
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
                logger.warning(f"Participant {sub} not found in column 'sub' in tracking file for {sub}, {ses}.")
            flag+=1
            return flag, None
        # Search session
        chans = chans[chans['ses']==ses]
        if chans.size == 0:
            if verbose>0:
                logger.warning(f"Session {ses} not found in column 'ses' in tracking file for {sub}, {ses}.")
            flag+=1
            return flag, None
        # Search channel
        chans = chans.filter(regex='chanset')
        chans = chans.filter(regex='^((?!rename).)*$')
        chans = chans.filter(regex='^((?!peaks).)*$')
        chans = chans.filter(regex='^((?!invert).)*$')
        
        chans = chans.dropna(axis=1, how='all')
        if len(chans.columns) == 0:
            if verbose>0:
                logger.warning(f"No channel set found in tracking file for {sub}, {ses}, skipping...")
            flag+=1
            return flag, None
    else:
        chans = chan
    
    if type(ref_chan) == type(DataFrame()):
        if verbose==2:
            logger.debug("Reading reference channel names from tracking file ")
        ref_chans = ref_chan[ref_chan['sub']==sub]
        if ref_chans.size == 0:
            if verbose>0:
                logger.warning(f"Participant not found in column 'sub' in tracking file for {sub}, {ses}.")
            flag+=1
            return flag, None
        ref_chans = ref_chans[ref_chans['ses']==ses]
        if ref_chans.size == 0:
            if verbose>0:
                logger.warning(f"Session not found in column 'ses' in tracking file for {sub}, {ses}.")
            flag+=1
            return flag, None
        ref_chans = ref_chans.filter(regex='refset')
        ref_chans = ref_chans.dropna(axis=1, how='all')
        if len(ref_chans.columns) == 0:
            if verbose>0:
                logger.warning(f"No reference channel set found in tracking file for {sub}, {ses}. Progressing without re-referencing...")
            ref_chans = []
    elif ref_chan:
        ref_chans = ref_chan
    else:
        ref_chans = []
    
    if type(chans) == list:
        if type(ref_chans) == DataFrame and len(ref_chans.columns) >1:
            chan = ref_chan[ref_chan['sub']==sub]
            chan = chan[chan['ses']==ses]
            chan = chan.filter(regex='chanset')
            chan = chan.filter(regex='^((?!rename).)*$')
            ref_chan=[]
            for c in chans:
                ref_chan.append([ref_chans[ref_chans.columns[x]].iloc[0] for x, y in enumerate(chan) if c in chan[y].iloc[0]][0])
            ref_chan = [char.split(x, sep=', ').tolist() for x in ref_chan]
            
            chanset = {chn:[ref_chan[i]] if isinstance(ref_chan[i],str) else ref_chan[i] for i,chn in enumerate(chans)}
            
        elif type(ref_chans) == DataFrame:
            ref_chans = ref_chans.to_numpy()[0]
            ref_chans = ref_chans.astype(str)
            ref_chans = char.split(ref_chans, sep=', ')
            ref_chans = [x for y in ref_chans for x in y]
            chanset = {chn:ref_chans for chn in chans}    
        else:
            chanset = {chn:[] for chn in chans}
    
    elif type(chans) == type(DataFrame()):
        if type(ref_chans) == DataFrame and len(ref_chans.columns) != len(chans.columns):
            logger.error(f"There must be the same number of channel sets and reference channel sets in 'tracking file, but for {sub}, {ses}, there were {len(chans.columns)} channel sets and {len(ref_chans.columns)} reference channel sets. For channel setup options, refer to documentation:")
            logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
            flag+=1
            return flag, None
        elif type(ref_chans) == DataFrame:
            ref_chans = ref_chans.to_numpy()[0]
            ref_chans = ref_chans.astype(str)
            ref_chans = char.split(ref_chans, sep=', ')
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
        logger.error("The variable 'chan' should be a [list] or definied in the 'chanset' column of tracking file - NOT a string.")
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
        oldchans = chans.filter(regex='chanset')
        oldchans = oldchans.filter(regex='^((?!rename).)*$')
        oldchans = oldchans.filter(regex='^((?!peaks).)*$')
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
    
    if type(oldchans) == type(DataFrame()):
        if type(newchans) == DataFrame and len(newchans.columns) != len(oldchans.columns):
            logger.warning(f"There must be the same number of channel sets and channel rename sets in tracking file, but for {sub}, {ses}, there were {len(oldchans.columns)} channel sets and {len(newchans.columns)} channel rename sets. For info on how to rename channels, refer to documentation:")
            logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
            logger.warning(f"Using original channel names for {sub}, {ses}...")
            return None
        
        oldchans=oldchans.to_numpy() 
        oldchans = oldchans.astype(str)
        oldchans = char.split(oldchans[0], sep=', ')
        oldchans = [x for y in oldchans for x in y]
        
        if type(newchans) == DataFrame:
            newchans = newchans.to_numpy()
            newchans = newchans.astype(str)
            newchans = char.split(newchans[0], sep=', ')
            newchans = [x for y in newchans for x in y]
        
        if len(oldchans) == len(newchans):
            newchans = {chn:newchans[i] for i,chn in enumerate(oldchans)}
        else:
            logger.warning(f"There must be the same number of original channel names and new renamed channels in tracking file, but for {sub}, {ses}, there were {len(oldchans)} old channel and {len(newchans)} new channel names. For info on how to rename channels, refer to documentation:")
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
                logger.warning(f"Participant {sub} not found in column 'sub' in tracking file for {sub}, {ses}.")
            flag+=1
            return flag, None
        # Search session
        chans = chans[chans['ses']==ses]
        if chans.size == 0:
            if verbose>0:
                logger.warning(f"Session {ses} not found in column 'ses' in tracking file for {sub}, {ses}.")
            flag+=1
            return flag, None
        # Search channel
        chans = chans.filter(regex='stagechan')
        chans = chans.dropna(axis=1, how='all')
        if len(chans.columns) == 0:
            if verbose>0:
                logger.warning(f"No stagechan found in tracking file for {sub}, {ses}, skipping...")
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
                logger.warning(f"Participant not found in column 'sub' in tracking file for {sub}, {ses}.")
            flag+=1
            return flag, None
        ref_chans = ref_chans[ref_chans['ses']==ses]
        if ref_chans.size == 0:
            if verbose>0:
                logger.warning(f"Session not found in column 'ses' in tracking file for {sub}, {ses}.")
            flag+=1
            return flag, None
        ref_chans = ref_chans.filter(regex='refset')
        ref_chans = ref_chans.dropna(axis=1, how='all')
        if len(ref_chans.columns) == 0:
            if verbose>0:
                logger.warning(f"No reference channel set found in tracking file for {sub}, {ses}. Progressing without re-referencing...")
            ref_chans = []
    elif ref_chan:
        ref_chans = ref_chan
    else:
        ref_chans = []
    
    if type(chans) == list:
        if type(ref_chans) == DataFrame and len(ref_chans.columns) >1:
            chan = ref_chan[ref_chan['sub']==sub]
            chan = chan[chan['ses']==ses]
            chan = chan.filter(regex='chanset')
            chan = chan.filter(regex='^((?!rename).)*$')
            ref_chan=[]
            for c in chans:
                ref_chan.append([ref_chans[ref_chans.columns[x]].iloc[0] for x, y in enumerate(chan) if c in chan[y].iloc[0]][0])
            ref_chan = [char.split(x, sep=', ').tolist() for x in ref_chan]
            
            chanset = {chn:[ref_chan[i]] if isinstance(ref_chan[i],str) else ref_chan[i] for i,chn in enumerate(chans)}
            
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
                logger.warning(f"There were >2 reference channel sets in 'tracking' file for {sub}, {ses}, we will just use the first set for automatic staging.")
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
        logger.error("The variable 'chan' should be a [list] or definied in the 'chanset' column of tracking file - NOT a string.")
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
                logger.warning(f"Participant not found in column 'sub' in tracking file for {sub}, {ses}.")
            flag+=1
            return flag, None
        # Search session
        chans = chans[chans['ses']==ses]
        if chans.size == 0:
            if verbose>0:
                logger.warning(f"Session not found in column 'ses' in tracking file for {sub}, {ses}.")
            flag+=1
            return flag, None
        # Search channel
        chans = chans.filter(regex='eog')
        chans = chans.dropna(axis=1, how='all')
        if len(chans.columns) == 0:
            if verbose>0:
                logger.warning(f"No stagechan found in tracking file for {sub}, {ses}...")
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
                logger.warning(f"Participant not found in column 'sub' in tracking file for {sub}, {ses}.")
            flag+=1
            return flag, None
        # Search session
        chans = chans[chans['ses']==ses]
        if chans.size == 0:
            if verbose>0:
                logger.warning(f"Session not found in column 'ses' in tracking file for {sub}, {ses}.")
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
        logger.info("For info how to use adap_bands = 'Manual' in detections, refer to documentation:")
        logger.info(" https://seapipe.readthedocs.io/en/latest/index.html")
        logger.info('-' * 10)
        return 'error'
    
    chans = track.filter(regex='chanset')
    chans = chans.filter(regex='^((?!rename).)*$')
    chans = chans.filter(regex='^((?!peaks).)*$')
    chans = chans.filter(regex='^((?!invert).)*$')
    chans = chans.dropna(axis=1, how='all')
    peaks = track.filter(regex='peaks')
    peaks = peaks.dropna(axis=1, how='all')
    
    if len(peaks.columns) == 0:
        logger.error("No spectral peaks have been provided in tracking file. Peaks will need to be detected.")
        logger.info("Check documentation for how to use adap_bands = 'Manual' in detections: https://seapipe.readthedocs.io/en/latest/index.html")
        logger.info('-' * 10)
        return 'error'
    elif len(peaks.columns) != len(chans.columns):
        logger.error("There must be the same number of channel sets and spectral peaks sets in tracking file")
        logger.info("Check documentation for how to use adap_bands = 'Manual' in detections: https://seapipe.readthedocs.io/en/latest/index.html")
        return 'error'
    
    sub = {}
    for c, col in enumerate(chans.columns):
        for r, row in enumerate(chans[col]):
            chs = reshape(char.split(str(row), sep=', '), (1,1))[0][0]
            pks = reshape(char.split(str(peaks.iloc[r,c]), sep=', '), (1,1))[0][0]  
            if len(chs) != len(pks) and 'nan' not in (pks):
                logger.warning(f"For {track['sub'][r]}, {track['ses'][r]} the number of channels provided ({len(chs)}) != the number of spectral peaks ({len(pks)}).")
                if not track['sub'][r] in sub.keys():
                    sub[track['sub'][r]] = [track['ses'][r]]
                else:
                    sub[track['sub'][r]].append(track['ses'][r])
            elif 'nan' in (pks) and 'nan' not in (chs):
                logger.warning(f"For {track['sub'][r]}, {track['ses'][r]} no peaks have been provided.")
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
        logger.info("For info how to use adap_bands = 'Manual' in detections, refer to documentation:")
        logger.info(" https://seapipe.readthedocs.io/en/latest/index.html")
        logger.info('-' * 10)
        return 'error'

    track = track[track['sub']==sub]
    if len(track.columns) == 0:
        logger.warning(f"Participant not found in column 'sub' in tracking file for {sub}, {ses}.")
        return None
    # Search session
    track = track[track['ses']==ses]
    if len(track.columns) == 0:
        logger.warning(f"Session not found in column 'ses' in tracking file for {sub}, {ses}.")
        return None
    
    # Search channel
    chans = track.filter(regex='chanset')
    chans = chans.filter(regex='^((?!rename).)*$')
    chans = chans.filter(regex='^((?!peaks).)*$')
    chans = chans.filter(regex='^((?!invert).)*$')
    chans = chans.dropna(axis=1, how='all')
    peaks = track.filter(regex='peaks')
    peaks = peaks.dropna(axis=1, how='all')
    
    if len(peaks.columns) == 0:
        logger.warning(f"No spectral peaks found in tracking file for {sub}, {ses}.")
        return None

    chans = chans.to_numpy()[0]
    chans = chans.astype(str)
    chans = char.split(chans, sep=', ')
    chans = [x for y in chans for x in y]
    
    peaks = peaks.to_numpy()[0]
    peaks = peaks.astype(str)
    peaks = char.split(peaks, sep=', ')
    peaks = [float(x) for y in peaks for x in y]
    
    freq = (peaks[chans.index(chan)] - adap_bw/2, 
            peaks[chans.index(chan)] + adap_bw/2)

    
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
        logger.warning(f'No specparams export file found for {sub}, {ses}, {ch}, {stage}, {band_limits}.')
        return None
    elif len(files) > 1:
        logger.warning(f'>1 specparams export files found for {sub}, {ses}, {ch}, {stage}, {band_limits} ?')
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
        logger.warning(f'No peaks found in export file for {sub}, {ses}, {ch}, {stage}, {band_limits}.')
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
            logger.warning(f"Participant not found in column 'sub' in tracking file for {sub}, {ses}.")
            return None
        # Search session
        chans = chans[chans['ses']==ses]
        if len(chans.columns) == 0:
            logger.warning(f"Session not found in column 'ses' in tracking file for {sub}, {ses}.")
            return None
        
        # Search channel
        chans = chans.filter(regex='chanset')
        inversion = chans.filter(regex='invert')
        inversion = inversion.dropna(axis=1, how='all')
        chans = chans.filter(regex='^((?!rename).)*$')
        chans = chans.filter(regex='^((?!peaks).)*$')
        chans = chans.filter(regex='^((?!invert).)*$')
        chans = chans.dropna(axis=1, how='all')
        
        if len(inversion.columns) == 0:
            logger.warning(f"No inversion info found in tracking file for {sub}, {ses}.")
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
            logger.warning(f"Error reading inversion info for {sub}, {ses}, {chan} - check documentation for how to provide information for inversion:")
            logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
            return None
    

def infer_ref(sub, ses, chan, logger, verbose=0):
    # verbose = 0 (error only)
    # verbose = 1 (warning only)
    # verbose = 2 (debug)
    
    if isinstance(chan, DataFrame):
        # Search participant
        chans = chan[chan['sub']==sub]
        if len(chans.columns) == 0:
            return None
        # Search session
        chans = chans[chans['ses']==ses]
        if len(chans.columns) == 0:
            return None
        # Search channels
        oldchans = chans.filter(regex='chanset')
        oldchans = oldchans.filter(regex='^((?!rename).)*$')
        oldchans = oldchans.filter(regex='^((?!peaks).)*$')
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
    
    if isinstance(oldchans, DataFrame):
        if isinstance(newchans, DataFrame) and len(newchans.columns) != len(oldchans.columns):
            if verbose>1:
                logger.warning(f"There must be the same number of channel sets and channel rename sets in tracking file, but for {sub}, {ses}, there were {len(oldchans.columns)} channel sets and {len(newchans.columns)} channel rename sets. For info on how to rename channels, refer to documentation:")
                logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
                logger.warning(f"Using original channel names for {sub}, {ses}...")
            return None
        
        oldchans=oldchans.to_numpy() 
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
                return None
        else:
            logger.warning(f"There must be the same number of original channel names and new renamed channels in tracking file, but for {sub}, {ses}, there were {len(oldchans)} old channel and {len(newchans)} new channel names. For info on how to rename channels, refer to documentation:")
            logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
            logger.warning(f"Using original channel names for {sub}, {ses}...")
            return None
    else:
        return None
    
    return ref_chan
        
    
    
    

