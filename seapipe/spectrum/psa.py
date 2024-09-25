# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 10:29:11 2021

@author: Nathan Cross
"""
from copy import deepcopy
from csv import reader
from datetime import datetime
from fooof import FOOOF
from fooof.analysis import get_band_peak_fm
from glob import glob
from itertools import product
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from numpy import (append, arange, asarray, ceil, concatenate, empty, floor, 
                   mean, nan, ones, pi, reshape, sqrt, stack, sum, where, zeros)
from openpyxl import Workbook
from os import listdir, mkdir, path, walk
from pandas import DataFrame, ExcelWriter, read_csv
from scipy.fftpack import next_fast_len
from wonambi import ChanFreq, Dataset
from wonambi.attr import Annotations
from wonambi.trans import (fetch, frequency, get_descriptives, export_freq, 
                           export_freq_band)
from ..utils.misc import bandpass_mne, laplacian_mne, notch_mne, notch_mne2
from ..utils.logs import create_logger
from ..utils.load import infer_ref, load_channels, load_sessions, rename_channels


def default_general_opts():
    general_opts = {'freq_full':True, 
                    'freq_band':True, 
                    'freq_plot':False,
                    'max_freq_plot':35,
                    'suffix':None,
                    'chan_grp_name':'eeg'}
    return general_opts

def default_frequency_opts():
    headers = ['SO', 
               'Delta', 
               'Theta', 
               'Alpha', 
               'Sigma', 
               'Low Beta', 
               'High Beta']
    bands = [(0.25, 1.25),  #SO
             (0.5, 4),      #delta
             (4, 7.75),     #theta
             (8, 11),       #alpha
             (11.25, 16),   #sigma
             (16.25, 19),   #low beta
             (19.25, 35)]   #high beta
    frequency_opts = {'bands': bands,
                      'headers': headers,
                      'frequency': (11, 16),
                      'output':'spectraldensity', 
                      'scaling':'power', 
                      'sides':'one', 
                      'taper':'hann', 
                      'detrend':'linear', 
                      'n_fft': None, 
                      'fast_fft': True, 
                      'duration': 4,
                      'overlap': 0.5,
                      'step': None,
                      'centend':'mean',
                      'log_trans': True,
                      'halfbandwidth': 3,
                      'NW': None}
    return frequency_opts

def default_filter_opts():
    filter_opts = {'laplacian': False, 
                   'oREF': None, 
                   'lapchan': None, 
                   'laplacian_rename': False, 
                   'renames': None, 
                   'montage': 'standard_alphabetic',
                   'notch': True, 
                   'notch_freq': 50, 
                   'notch_harmonics': False,
                   'bandpass': True,
                   'highpass': 0.25,
                   'lowpass': 40,
                   'dcomplex': 'hilbert',
                   'filtcycle': (3, 6),
                   'width': 7}
    return filter_opts

def default_epoch_opts():
    epoch_opts = {'epoch': None,
                  'reject_epoch': True,
                  'reject_artf': True,
                  'min_dur':0,
                  'epoch_dur': 30,
                  'epoch_overlap': 0.5,
                  'epoch_step': None,
                  'concat_signal':True}
    
    return epoch_opts

def default_event_opts():
    event_opts = {'evt_type': None,
                  'event_chan': None,
                  'buffer': 2,
                  'concat_events':True}
    return event_opts

def default_fooof_opts():
    fooof_opts = {'psd_dur': 2, 
                  'peak_width_limits': [1, 12], 
                  'max_n_peaks': 8,
                  'min_peak_amplitude': 0.0,
                  'peak_threshold': 2.0,
                  'freq_range': [0.5, 35],
                  'select_highest': True,
                  'thresh_param': 'PW',
                  'bands_fooof': None, 
                  'thresh_select': None}
    return fooof_opts

def default_norm_opts():
    norm_opts = {'norm_cat': (1, 1, 1, 0),
                 'norm_evt_type': ['norm'], 
                 'norm_stage': None, 
                 'norm_epoch': None}
    return norm_opts


class Spectrum:
    
    """Design power spectral analyses on segmented data.

    Parameters
    ----------
    method : str
        one of the predefined methods
    frequency : tuple of float
        low and high frequency of frequency band
    duration : tuple of float
        min and max duration of spindles
    merge : bool
        if True, then after events are detected on every channel, events on 
        different channels that are separated by less than min_interval will be 
        merged into a single event, with 'chan' = the chan of the earlier-onset 
        event.
        
    Functions
    ----------
    fooof_it : 
        calculation of spectral parameters.
    powerspec_it : 
        Call to calculate power spectral analysis for each of:
            <sub>, <ses>, <chan>, <segment>
    powerspec_summary_full : 
        for power spectrum statistics. Returns value for whole spectrum, in 
        N/fs bins.
    powerspec_summary_bands : 
        for power spectrum statistics. Returns value for whole spectrum, 
        for pre-defined bands of interest.
        
    Notes
    -----
    See individual functions for other attribute descriptions.
    """ 
    
    def __init__(self, rec_dir, xml_dir, out_dir, log_dir, chan, ref_chan, 
                 grp_name, stage, cat, rater = None, cycle_idx = None,
                 subs = 'all', sessions = 'all', tracking = None):
        
        self.rec_dir = rec_dir
        self.xml_dir = xml_dir
        self.out_dir = out_dir
        self.log_dir = log_dir
        
        self.subs = subs
        self.sessions = sessions
        
        self.chan = chan
        self.ref_chan = ref_chan
        self.grp_name = grp_name
        self.rater = rater
        
        self.stage = stage
        self.cycle_idx = cycle_idx
        self.cat = cat
        
        if tracking == None:
            tracking = {}
        self.tracking = tracking
            

    
    def fooof_it(self, general_opts, frequency_opts, filter_opts, epoch_opts, 
                       event_opts, fooof_opts, filetype = '.edf'):
        

        '''
        FOOOF is a fast, efficient, and physiologically-informed tool to 
        parameterize neural power spectra.
        
        Outputs
        -------
        ap_params: 
                  Parameters for the aperiodic component, reflecting 1/f like 
                  characteristics, including: 
                                    - offset
                                    - slope
        pk_params:
                  Parameters for periodic components (putative oscillations), 
                  as peaks rising above the aperiodic component, including:
                                    - peak frequency (in Hz)
                                    - peak bandwidth
                                    - peak amplitude (in µV^2)
        
        Notes
        -----
        See: https://fooof-tools.github.io/fooof/index.html for more info.
        
        '''
        
        ### 0.a. Set up logging
        logger = create_logger('Specparam')
        logger.info('')
        tracking = self.tracking
        flag = 0
        
        ### 0.b. Set up organisation of export
        suffix = general_opts['suffix']
        freq_bw = frequency_opts['frequency']
        if self.cat[0] + self.cat[1] == 2:
            model = 'whole_night'
            logger.debug(f'Parameterizing power spectrum in range {freq_bw[0]}-{freq_bw[1]} Hz for the whole night.')
        elif self.cat[0] + self.cat[1] == 0:
            model = 'stage*cycle'
            logger.debug(f'Parameterizing power spectrum in range {freq_bw[0]}-{freq_bw[1]} Hz per stage and cycle separately.')
        elif self.cat[0] == 0:
            model = 'per_cycle'
            logger.debug(f'Parameterizing power spectrum in range {freq_bw[0]}-{freq_bw[1]} Hz per cycle separately.')
        elif self.cat[1] == 0:
            model = 'per_stage'  
            logger.debug(f'Parameterizing power spectrum in range {freq_bw[0]}-{freq_bw[1]} Hz per stage separately.')
        if 'cycle' in model and self.cycle_idx == None:
            logger.info('')
            logger.critical("To run cycles separately (i.e. cat[0] = 0), cycle_idx cannot be 'None'")
            return
        cat = tuple((self.cat[0],self.cat[1],1,1)) # force concatenation of discontinuous & events
        logger.info('')
        logger.debug(r"""
                     
                                  |
                                  |
                                  |  
                                  |  .
                                  |   `~.
                        (µV2)     |      `~.       FOOOF !
                                  |         `^.
                                  |            `~.
                                  |               `•._
                                  |                   `~¬.…_
                                  |_____________________________
                                              (Hz)
                                  
                                                    """,)
        
        # 1.a. Define model parameters
        if fooof_opts is None:
            logger.critical('Options not set for FOOOF')
            return
        else:
            fm = FOOOF(fooof_opts['peak_width_limits'], fooof_opts['max_n_peaks'], 
                       fooof_opts['min_peak_amplitude'], fooof_opts['peak_threshold'])
            def gaussian_integral(a, c):
                """ Returns definite integral of a gaussian function with height a and
                standard deviation c."""
                return sqrt(2) * a * abs(c) * sqrt(pi)
        
        # 1.b. Prepare Bands
        if fooof_opts['bands_fooof'] is not None:
            bands = fooof_opts['bands_fooof']
        else:
            stp = min(fooof_opts['peak_width_limits'])
            low = int(floor(fooof_opts['freq_range'][0]))
            hi = int(ceil(fooof_opts['freq_range'][1]))
            bands = [(x,x + stp) for x in range(low,hi)]
        
        # 2.a. Check for output folder, if doesn't exist, create
        if path.exists(self.out_dir):
                logger.debug(f"Output directory: {self.out_dir} exists")
        else:
            mkdir(self.out_dir)
        
        # 2.b. Get subjects
        subs = self.subs
        if isinstance(subs, list):
            None
        elif subs == 'all':
                subs = next(walk(self.xml_dir))[1]
        else:
            logger.info('')
            logger.error("'subs' must either be an array of Participant IDs or = 'all' ")
            return
        subs.sort()
        
        # 3.a. Begin loop through participants
        for p, sub in enumerate(subs):
            
            # b. Begin loop through sessions
            flag, sessions = load_sessions(sub, self.sessions, self.rec_dir, flag, 
                                     logger, verbose=2)

            for v, ses in enumerate(sessions):
                logger.info('')
                logger.debug(f'Commencing {sub}, {ses}')
                
                ## Define files
                rdir = f'{self.rec_dir}/{sub}/{ses}/eeg/'
                xdir = f'{self.xml_dir}/{sub}/{ses}/'
                
                try:
                    edf_file = [x for x in listdir(rdir) if x.endswith(filetype)]
                    dset = Dataset(rdir + edf_file[0])
                except:
                    logger.warning(f' No input {filetype} file in {rdir}. Skipping...')
                    flag+=1
                    break
                
                try:
                    xml_file = [x for x in listdir(xdir) if x.endswith('.xml')]
                except:
                    logger.warning(f"No input annotations file in {xdir} or path doesn't exist. Skipping...")
                    flag+=1
                    break
                
                ### Define output path
                if not path.exists(self.out_dir):
                    mkdir(self.out_dir)
                if not path.exists(f'{self.out_dir}/{sub}'):
                    mkdir(f'{self.out_dir}/{sub}')
                if not path.exists(f'{self.out_dir}/{sub}/{ses}'):
                    mkdir(f'{self.out_dir}/{sub}/{ses}')
                outpath = f'{self.out_dir}/{sub}/{ses}'
                
                ## Now import data
                dset = Dataset(rdir + edf_file[0])
                annot = Annotations(xdir + xml_file[0], rater_name=self.rater)
                 
                ### get cycles
                if self.cycle_idx is not None:
                    all_cycles = annot.get_cycles()
                    cycle = [all_cycles[i - 1] for i in self.cycle_idx if i <= len(all_cycles)]
                else:
                    cycle = None
                    
                ### if event channel only, specify event channels
                # 4.d. Channel setup
                flag, chanset = load_channels(sub, ses, self.chan, 
                                              self.ref_chan, flag, logger)
                if not chanset:
                    flag+=1
                    break
                newchans = rename_channels(sub, ses, self.chan, logger)

                # get segments
                for c, ch in enumerate(chanset):
                    logger.debug(f"Reading data for {ch}:{'/'.join(chanset[ch])}")
                    segments = fetch(dset, annot, cat = cat, 
                                     evt_type = event_opts['evt_type'], 
                                     stage = self.stage, cycle=cycle,  
                                     epoch = epoch_opts['epoch'], 
                                     epoch_dur = epoch_opts['epoch_dur'], 
                                     epoch_overlap = epoch_opts['epoch_overlap'], 
                                     epoch_step = epoch_opts['epoch_step'], 
                                     reject_epoch = epoch_opts['reject_epoch'], 
                                     reject_artf = epoch_opts['reject_artf'],
                                     min_dur = epoch_opts['min_dur'])
                    
                    # 5.b Rename channel for output file (if required)
                    if newchans:
                        fnamechan = newchans[ch]
                    else:
                        fnamechan = ch
                    
                    # Read data
                    segments.read_data(chan = ch, ref_chan = chanset[ch])
                    if len(segments)==0:
                        logger.warning(f"No valid data found for {ch}:{'/'.join(chanset[ch])}.")
                        flag+=1
                        continue
                    
                    for sg, seg in enumerate(segments):
                        logger.debug(f'Analysing segment {sg+1} of {len(segments)}')
                        out = dict(seg)
                        data = seg['data']
                        timeline = data.axis['time'][0]
                        out['start'] = timeline[0]
                        out['end'] = timeline[-1]
                        out['duration'] = len(timeline) / data.s_freq
                        
                        if frequency_opts['fast_fft']:
                            n_fft = frequency_opts['n_fft']
                        else:
                            n_fft = None
                    
                        Fooofxx = frequency(data, output=frequency_opts['output'], 
                                        scaling=frequency_opts['scaling'], 
                                        sides=frequency_opts['sides'], 
                                        taper=frequency_opts['taper'],
                                        halfbandwidth=frequency_opts['halfbandwidth'], 
                                        NW=frequency_opts['NW'],
                                        duration=frequency_opts['duration'], 
                                        overlap=frequency_opts['overlap'], 
                                        step=frequency_opts['step'],
                                        detrend=frequency_opts['detrend'], 
                                        n_fft=n_fft, 
                                        log_trans=False, 
                                        centend=frequency_opts['centend'])
                        
                        freqs = Fooofxx.axis['freq'][0]     
                        fooof_powers = zeros((len(bands)))
                        fooof_ap_params = zeros((len(bands), 2))
                        fooof_pk_params = ones((len(bands), 9)) * nan
                        fm.fit(freqs, Fooofxx.data[0][0], fooof_opts['freq_range'])
                        
                        for j, band in enumerate(bands):
                            fp = get_band_peak_fm(fm, band, fooof_opts['select_highest'],
                                                  threshold=fooof_opts['thresh_select'],
                                                  thresh_param=fooof_opts['thresh_param'])
                            if fp.ndim == 1:
                                fooof_powers[j] = gaussian_integral(fp[1], fp[2])
                            else:
                                pwr = asarray([gaussian_integral(fp[x, 1], fp[x, 2]) \
                                               for x in range(fp.shape[0])]).sum()
                                fooof_powers[j] = pwr
                                
                            # get fooof aperiodic parameters
                            fooof_ap_params[j, :] = fm.aperiodic_params_
                            
                            # get fooof peak parameters
                            fp = get_band_peak_fm(fm, band, False,
                                                  threshold=fooof_opts['thresh_select'],
                                                  thresh_param=fooof_opts['thresh_param'])
                            if fp.ndim == 1:
                                fooof_pk_params[j, :3] = fp
                            else:
                                n_peaks = min(fp.shape[0], 3)
                                fooof_pk_params[j, :n_peaks * 3] = fp[:n_peaks, 
                                                                      :].ravel()
                        out['fooof_powers'] = fooof_powers
                        out['fooof_ap_params'] = fooof_ap_params
                        out['fooof_pk_params'] = fooof_pk_params
                                
                        
    
                        seg_info = ['Start time', 'End time', 'Duration', 'Stitches', 
                                    'Stage', 'Cycle', 'Event type', 'Channel']
                        band_hdr = [f'{b1}-{b2} Hz' for b1, b2 in bands]
                        pk_params_hdr = ['peak1_CF', 'peak1_PW', 'peak1_BW', 
                                         'peak2_CF', 'peak2_PW', 'peak2_BW', 
                                         'peak3_CF', 'peak3_PW', 'peak3_BW', ]
                        ap_params_hdr = ['Offset', 'Exponent']
                        band_pk_params_hdr = ['_'.join((b, p)) for b in band_hdr 
                                              for p in pk_params_hdr]
                        band_ap_params_hdr = ['_'.join((b, p)) for b in band_hdr 
                                              for p in ap_params_hdr]
                        one_record = zeros((1, 
                                            (len(seg_info) + len(bands) + 
                                            len(band_pk_params_hdr) + 
                                            len(band_ap_params_hdr))),
                                            dtype='O')
                        one_record[0, :] = concatenate((asarray([
                                                                out['start'],
                                                                out['end'],
                                                                out['duration'], 
                                                                out['n_stitch'], # number of concatenated segments minus 1
                                                                out['stage'],
                                                                out['cycle'],
                                                                out['name'], # event type
                                                                ch,
                                                                ]),
                                                                out['fooof_powers'],
                                                                out['fooof_pk_params'].ravel(),
                                                                out['fooof_ap_params'].ravel(),
                                                                ))
                        
                    ### Output filename ###
                    if model == 'whole_night':
                        stagename = '-'.join(self.stage)
                        outputfile = f'{outpath}/{sub}_{ses}_{fnamechan}_{stagename}_specparams_{suffix}.csv'
                    elif model == 'stage*cycle':    
                        outputfile = f'{outpath}/{sub}_{ses}_{fnamechan}_{self.stage[sg]}_cycle{self.cycle_idx[sg]}_specparams_{suffix}.csv'
                    elif model == 'per_stage':
                        outputfile = f'{outpath}/{sub}_{ses}_{fnamechan}_{self.stage[sg]}_specparams_{suffix}.csv'
                    elif model == 'per_cycle':
                        stagename = '-'.join(self.stage)
                        outputfile = f'{outpath}/{sub}_{ses}_{fnamechan}_{stagename}_cycle{self.cycle_idx[sg]}_specparams_{suffix}.csv'

                    logger.debug(f'Saving {outputfile}')
                    df = DataFrame(data = one_record, 
                                   columns = (seg_info + band_hdr + band_pk_params_hdr
                                              + band_ap_params_hdr))
                    df.to_csv(outputfile)
                    
                    # Update tracking TO DO
                    #update = datetime.fromtimestamp(path.getmtime(outputfile)).strftime("%m-%d-%Y, %H:%M:%S")
        ### 3. Check completion status and print
        if flag == 0:
            logger.info('')
            logger.debug('Specparam finished without ERROR.')  
        else:
            logger.info('')
            logger.warning('Specparam finished with WARNINGS. See log for details.')
        return 

    
    
    def powerspec_it(self, general_opts, frequency_opts, filter_opts, epoch_opts, 
                           event_opts, norm, norm_opts, filetype = '.edf'):
        
        ### 0.a. Set up logging
        logger = create_logger('Power Spectral Analysis')
        tracking = self.tracking
        logger.info('')
        flag = 0
        
        ### 0.b. Set up organisation of export
        suffix = general_opts['suffix']
        freq_bw = frequency_opts['frequency']
        if self.cat[0] + self.cat[1] == 2:
            model = 'whole_night'
            logger.debug(f'Calculating power spectrum in range {freq_bw[0]}-{freq_bw[1]} Hz for the whole night.')
        elif self.cat[0] + self.cat[1] == 0:
            model = 'stage*cycle'
            logger.debug(f'Calculating power spectrum in range {freq_bw[0]}-{freq_bw[1]} Hz per stage and cycle separately.')
        elif self.cat[0] == 0:
            model = 'per_cycle'
            logger.debug(f'Calculating power spectrum in range {freq_bw[0]}-{freq_bw[1]} Hz per cycle separately.')
        elif self.cat[1] == 0:
            model = 'per_stage'  
            logger.debug(f'Calculating power spectrum in range {freq_bw[0]}-{freq_bw[1]} Hz per stage separately.')
        
        if self.cat[3] == 0:
            ev_model = True
            logger.debug(f"Calculating power spectrum in range {freq_bw[0]}-{freq_bw[1]} Hz per event: {event_opts['evt_type']} separately.")
        else:
            ev_model = False
        if 'cycle' in model and self.cycle_idx == None:
            logger.info('')
            logger.critical("To run cycles separately (i.e. cat[0] = 0), cycle_idx cannot be 'None'")
            return
        cat = self.cat
        #cat = tuple((self.cat[0],self.cat[1],1,1)) # force concatenation of discontinuous & events
        
        
        logger.info('')
        logger.debug(r"""
                     
                                  |
                                  |
                                  |  ,^\
                                  | /   \           Power Spectral
                                  |     `.             Analysis
                        (µV2)     |      `~.       
                                  |         `./\
                                  |             `~.
                                  |                `•._
                                  |                    `~¬.…_._._
                                  |_______________________________
                                              (Hz)
                                  
                                                    """,)
        
        
        # 1.a. Check for output folder, if doesn't exist, create
        if path.exists(self.out_dir):
                logger.debug(f"Output directory: {self.out_dir} exists")
        else:
            mkdir(self.out_dir)
        
        # 2.b. Get subjects
        subs = self.subs
        if isinstance(subs, list):
            None
        elif subs == 'all':
                subs = next(walk(self.rec_dir))[1]
        else:
            logger.info('')
            logger.critical("'subs' must either be an array of Participant IDs or = 'all' ")
            return
        subs.sort()
        
        # 3.a. Begin loop through participants
        for p, sub in enumerate(subs):
            tracking[f'{sub}'] = {}
            
            # b. Begin loop through sessions
            sessions = self.sessions
            if sessions == 'all':
                sessions = listdir(self.rec_dir + '/' + sub)
                sessions = [x for x in sessions if not '.' in x] 

            for v, ses in enumerate(sessions):
                logger.info('')
                logger.debug(f'Commencing {sub}, {ses}')
                tracking[f'{sub}'][f'{ses}'] = {'powerspec':{}} 
                
                ## Define files
                rdir = f'{self.rec_dir}/{sub}/{ses}/eeg/'
                xdir = f'{self.xml_dir}/{sub}/{ses}/'
                
                try:
                    edf_file = [x for x in listdir(rdir) if x.endswith(filetype)]
                    dset = Dataset(rdir + edf_file[0])
                except:
                    logger.warning(f' No input {filetype} file in {rdir}. Skipping...')
                    break
                
                xml_file = [x for x in listdir(xdir) if x.endswith('.xml')]            
                
                ### Define output path
                if not path.exists(self.out_dir):
                    mkdir(self.out_dir)
                if not path.exists(f'{self.out_dir}/{sub}'):
                    mkdir(f'{self.out_dir}/{sub}')
                if not path.exists(f'{self.out_dir}/{sub}/{ses}'):
                    mkdir(f'{self.out_dir}/{sub}/{ses}')
                outpath = f'{self.out_dir}/{sub}/{ses}'
                
                ## Now import data
                dset = Dataset(rdir + edf_file[0])
                annot = Annotations(xdir + xml_file[0], rater_name=self.rater)
                 
                ### get cycles
                if self.cycle_idx is not None:
                    all_cycles = annot.get_cycles()
                    cycle = [all_cycles[i - 1] for i in self.cycle_idx if i <= len(all_cycles)]
                else:
                    cycle = None
                    
                # 4.d. Channel setup          
                flag, chanset = load_channels(sub, ses, self.chan, 
                                              self.ref_chan, flag, logger)
                if not chanset:
                    logger.error('Problem loading channels from tracking sheet')
                    break
            
                newchans = rename_channels(sub, ses, self.chan, logger)
                filter_opts['oREF'] = infer_ref(sub, ses, self.chan, logger)    
                
                # Loop through channels
                for c, ch in enumerate(chanset):
                    
                    chan_full = [ch + ' (' + general_opts['chan_grp_name'] + ')']

                    # 5.b Rename channel for output file (if required)
                    if newchans:
                        fnamechan = newchans[ch]
                    else:
                        fnamechan = ch
                    
                    # Normalisation (if requested)
                    if norm == 'baseline':
                        logger.debug('Extracting data for baseline normalisation...')
                        norm_seg = fetch(dset, annot, cat=norm_opts['norm_cat'], 
                                         evt_type=norm_opts['norm_evt_type'],
                                         stage=norm_opts['norm_stage'], 
                                         epoch=norm_opts['norm_epoch'])
                        if not norm_seg:
                            logger.warning('No valid baseline data found. Skipping...')
                            continue
                        if filter_opts['laplacian']:
                            try:
                                norm_seg.read_data(filter_opts['lapchan'], chanset[ch]) 
                                logger.debug("Applying Laplacian filtering to baseline data...")
                                laplace_flag = True
                            except:
                                logger.error(f"Channels listed in filter_opts['lapchan']: {filter_opts['lapchan']} are not found in recording.")
                                logger.warning("Laplacian filtering will NOT be run for BASELINE data. Check parameters under: filter_opts['lapchan']")
                                norm_seg.read_data(ch, chanset[ch])
                                laplace_flag = False
                                flag += 1
                        else:
                            norm_seg.read_data(ch, chanset[ch])            
                        all_nSxx = []
                        for seg in norm_seg:
                            normdata = seg['data']
                            if laplace_flag:
                                normdata.data[0] = laplacian_mne(normdata, 
                                                         filter_opts['oREF'], 
                                                         channel=chan_full, 
                                                         ref_chan=chanset[ch], 
                                                         laplacian_rename=filter_opts['laplacian_rename'], 
                                                         renames=filter_opts['renames'])
                            Sxx = frequency(normdata, output=frequency_opts['output'], 
                                            scaling=frequency_opts['scaling'],
                                            sides=frequency_opts['sides'], 
                                            taper=frequency_opts['taper'],
                                            halfbandwidth=frequency_opts['halfbandwidth'], 
                                            NW=frequency_opts['NW'],
                                            duration=frequency_opts['duration'], 
                                            overlap=frequency_opts['overlap'], 
                                            step=frequency_opts['step'],
                                            detrend=frequency_opts['detrend'], 
                                            n_fft=frequency_opts['n_fft'], 
                                            log_trans=frequency_opts['log_trans'], 
                                            centend=frequency_opts['centend'])
                            all_nSxx.append(Sxx)
                            nSxx = ChanFreq()
                            nSxx.s_freq = Sxx.s_freq
                            nSxx.axis['freq'] = Sxx.axis['freq']
                            nSxx.axis['chan'] = Sxx.axis['chan']
                            nSxx.data = empty(1, dtype='O')
                            nSxx.data[0] = empty((Sxx.number_of('chan')[0],
                                     Sxx.number_of('freq')[0]), dtype='f')
                            nSxx.data[0] = mean(
                                    stack([x()[0] for x in all_nSxx], axis=2), axis=2)
                    
                    
                    
                    ## Get segments of data
                    logger.debug(f"Reading data for {ch}:{'/'.join(chanset[ch])}")
                    
                    segments = fetch(dset, annot, cat = cat, 
                                     evt_type = event_opts['evt_type'], 
                                     stage = self.stage, cycle=cycle,   
                                     epoch = epoch_opts['epoch'], 
                                     epoch_dur = epoch_opts['epoch_dur'], 
                                     epoch_overlap = epoch_opts['epoch_overlap'], 
                                     epoch_step = epoch_opts['epoch_step'], 
                                     reject_epoch = epoch_opts['reject_epoch'], 
                                     reject_artf = epoch_opts['reject_artf'],
                                     min_dur = epoch_opts['min_dur'],
                                     buffer = event_opts['buffer'])
                    if len(segments)==0:
                        logger.warning(f"No valid data found for {sub}, {ses}, {self.stage}, {cycle}.")
                        continue
                    
                    # Read data
                    if filter_opts['laplacian']:
                        try:
                            segments.read_data(filter_opts['lapchan'], chanset[ch]) 
                            logger.debug("Applying Laplacian filtering to data...")
                            laplace_flag = True
                        except:
                            logger.error(f"Channels listed in filter_opts['lapchan']: {filter_opts['lapchan']} are not found in recording.")
                            logger.warning("Laplacian filtering will NOT be run. Check parameters under: filter_opts['lapchan']")
                            segments.read_data(ch, chanset[ch])
                            laplace_flag = False
                            flag += 1
                    else:
                        segments.read_data(ch, chanset[ch])
                    
                    # Loop over segments and apply transforms
                    for sg, seg in enumerate(segments):
                        logger.debug(f'Analysing segment {sg} of {len(segments)}')
                        out = dict(seg)
                        data = seg['data']
                        timeline = data.axis['time'][0]
                        out['start'] = timeline[0]
                        out['end'] = timeline[-1]
                        out['duration'] = len(timeline) / data.s_freq
                        
                        if not frequency_opts['fast_fft']:
                            n_fft = n_fft = 256*8
                        else:
                            n_fft = frequency_opts['n_fft']
                        
                        if filter_opts['laplacian']:
                            selectchans = filter_opts['lapchan']
                        else:
                            selectchans = ch
                        
                        # Notch filters
                        if filter_opts['notch']:
                            logger.debug(f"Applying notch filtering: {filter_opts['notch_freq']}")
                            data.data[0] = notch_mne(data, oREF=filter_opts['oREF'], 
                                                        channel=selectchans, 
                                                        freq=filter_opts['notch_freq'],
                                                        rename=filter_opts['laplacian_rename'],
                                                        renames=filter_opts['renames'],
                                                        montage=filter_opts['montage'])
                            
                        if filter_opts['notch_harmonics']: 
                            logger.debug('Applying notch harmonics filtering.')
                            logger.debug(f'{selectchans}')
                            data.data[0] = notch_mne2(data, oREF=filter_opts['oREF'], 
                                                      channel=selectchans, 
                                                      rename=filter_opts['laplacian_rename'],
                                                      renames=filter_opts['renames'],
                                                      montage=filter_opts['montage'])    
                        
                        # Bandpass filters
                        if filter_opts['bandpass']:
                            logger.debug(f"Applying bandpass filtering: {filter_opts['highpass']} - {filter_opts['lowpass']} Hz")
                            data.data[0] = bandpass_mne(data, oREF=filter_opts['oREF'], 
                                                      channel=selectchans,
                                                      highpass=filter_opts['highpass'], 
                                                      lowpass=filter_opts['lowpass'], 
                                                      rename=filter_opts['laplacian_rename'],
                                                      renames=filter_opts['renames'],
                                                      montage=filter_opts['montage'])
                        
                        # Laplacian transform
                        if filter_opts['laplacian'] and laplace_flag:
                            logger.debug('Applying Laplacian filtering.')
                            data.data[0] = laplacian_mne(data, 
                                                 filter_opts['oREF'], 
                                                 channel=selectchans, 
                                                 ref_chan=chanset[ch], 
                                                 laplacian_rename=filter_opts['laplacian_rename'], 
                                                 renames=filter_opts['renames'],
                                                 montage=filter_opts['montage'])
                            data.axis['chan'][0] = asarray([x for x in chanset])
                            selectchans = ch
                        
                        ### Frequency transformation
                        Sxx = frequency(data, output=frequency_opts['output'], 
                                        scaling=frequency_opts['scaling'], 
                                        sides=frequency_opts['sides'], 
                                        taper=frequency_opts['taper'],
                                        halfbandwidth=frequency_opts['halfbandwidth'], 
                                        NW=frequency_opts['NW'],
                                        duration=frequency_opts['duration'], 
                                        overlap=frequency_opts['overlap'], 
                                        step=frequency_opts['step'],
                                        detrend=frequency_opts['detrend'], 
                                        n_fft=n_fft, 
                                        log_trans=frequency_opts['log_trans'], 
                                        centend=frequency_opts['centend'])
                        
                        # Baseline normalisation
                        if norm:
                            logger.debug('Applying baseline normalisation')
                            for j, ch in enumerate(Sxx.axis['chan'][0]):
                                dat = Sxx.data[0][j,:]
                                sf = Sxx.axis['freq'][0]
                                f_res = sf[1] - sf[0] # frequency resolution
                                if norm == 'integral':
                                    norm_dat = sum(dat) * f_res # integral by midpoint rule
                                else:
                                    norm_dat = nSxx(chan=ch)[0]
                                Sxx.data[0][j,:] = dat / norm_dat
                                
                        out['data'] = Sxx

                        # Export and save data
                        ### Output filename ###
                        if not ev_model:
                            if model == 'whole_night':
                                stagename = '-'.join(self.stage)
                                outputfile = f'{outpath}/{sub}_{ses}_{fnamechan}_{stagename}**{suffix}.csv'
                            elif model == 'stage*cycle':    
                                outputfile = f'{outpath}/{sub}_{ses}_{fnamechan}_{self.stage[sg]}_cycle{self.cycle_idx[sg]}**{suffix}.csv'
                            elif model == 'per_stage':
                                outputfile = f'{outpath}/{sub}_{ses}_{fnamechan}_{self.stage[sg]}**{suffix}.csv'
                            elif model == 'per_cycle':
                                stagename = '-'.join(self.stage)
                                outputfile = f'{outpath}/{sub}_{ses}_{fnamechan}_{stagename}_cycle{self.cycle_idx[sg]}**{suffix}.csv'
                            if general_opts['freq_full']:
                                save_freq([out], freq_bw, general_opts, outputfile, logger)
                            if general_opts['freq_band']:
                                save_bands([out], freq_bw, general_opts, frequency_opts, outputfile, logger)
                            if general_opts['freq_plot']:
                                save_plot([out], general_opts, outputfile, logger)
                        elif sg == 0:
                            out_full = [deepcopy(out)]
                        else:
                            out_full.append(out)
                    
                    # If analyses set to be per event        
                    if ev_model:
                        ev_out = []
                        outputfile = [] #f"{outputfile.split('*')[0]}_{event_opts['evt_type'][0]}_**{outputfile.split('*')[2]}"
                        if model == 'whole_night':
                            outputfile.append(f"{outpath}/{sub}_{ses}_{fnamechan}_{stagename}_{event_opts['evt_type'][0]}**{suffix}.csv")
                            ev_out.append(out_full)
                        if model == 'stage*cycle':
                            for (stage,cycle) in product(self.stage,self.cycle_idx):
                                outputfile.append(f"{outpath}/{sub}_{ses}_{fnamechan}_{stage}_cycle{cycle}_{event_opts['evt_type'][0]}**{suffix}.csv")
                                out = [x for x in out_full if stage in x['stage']]
                                out = [x for x in out if cycle in x['cycle']]
                                ev_out.append(out)
                        if model == 'per_stage':
                            for stage in self.stage:
                                outputfile.append(f'{outpath}/{sub}_{ses}_{fnamechan}_{stage}**{suffix}.csv')
                                out = [x for x in out_full if stage in x['stage']]
                                ev_out.append(out)
                        if model == 'per_cycle':
                            for cycle in self.cycle_idx:
                                outputfile.append(f'{outpath}/{sub}_{ses}_{fnamechan}_{stagename}_cycle{cycle}**{suffix}.csv') 
                                out = [x for x in out if cycle in x['cycle']]
                                ev_out.append(out)
                        for out, outfile in zip(ev_out, outputfile):
                            if general_opts['freq_full']:
                                save_freq(out, freq_bw, general_opts, outfile, logger)
                            if general_opts['freq_band']:
                                save_bands(out, freq_bw, general_opts, frequency_opts, outfile, logger)
                            if general_opts['freq_plot']:
                                save_plot(out, general_opts, outfile, logger)     
                        
        ### 3. Check completion status and print
        if flag == 0:
            logger.debug('Power spectral analyses finished without ERROR.')  
        else:
            logger.warning('Power spectral analyses finished with WARNINGS. See log for details.')
        
        #self.tracking = tracking   ## TO UPDATE - FIX TRACKING
        
        return 
           
        
    def powerspec_summary(self, chan, general_opts, frequency_opts, filter_opts, 
                                epoch_opts, event_opts, logger):
        
        ### 0.a. Set up logging
        logger = create_logger('Power Spectral Summary')
        logger.info('')
        
        logger.info('')
        logger.debug(r""" Summarising Power Spectrum Data...
                     
                                     .^.
                                    / |\\_
                                   / | \ ¯;
                                   `|  \ /
                                    |  \ \
                                    |  \ \
                                    \__\_|
                                    (    )
                                    (0 ) )
                                   // ||\\
                                 /(( // ||
                               // \\))||\\
                              )) //|| ||))
                                (( )) |//
                                  //  ((
                                  
                    Spectral Quantification Unified In a Dataset
                    (S.Q.U.I.D)
                    
                                                    """,)
        
        ### 0.b. Set up output file params
        suffix = general_opts['suffix']
        frequency = frequency_opts['frequency']
        flag = 0
        if general_opts['freq_full']:
            variables_full = [str(x) for x in (arange(0.25,filter_opts['lowpass']*4,0.25))]
            idx_data_full = list(range(9,8+filter_opts['lowpass']*4))
        if general_opts['freq_band']:
            col_headers = frequency_opts['headers'] 
            variables_band = col_headers
            idx_data_band = list(range(9,9+len(col_headers)))
        
        ### 0.c. Set up organisation of export
        if self.cat[0] == 0 and self.cycle_idx == None:
            logger.info('')
            logger.critical("To run cycles separately (i.e. cat[0] = 0), cycle_idx cannot be 'None'")
            return
        if self.cat[0] + self.cat[1] == 2:
            model = 'whole_night'
            logger.debug(f'Summarising power spectrum in range {frequency[0]}-{frequency[1]} Hz for the whole night.')
        elif self.cat[0] + self.cat[1] == 0:
            model = 'stage*cycle'
            logger.debug(f'Calculating power spectrum in range {frequency[0]}-{frequency[1]} Hz per stage and cycle separately.')
        elif self.cat[0] == 0:
            model = 'per_cycle'
            logger.debug(f'Calculating power spectrum in range {frequency[0]}-{frequency[1]} Hz per cycle separately.')
        elif self.cat[1] == 0:
            model = 'per_stage'  
            logger.debug(f'Calculating power spectrum in range {frequency[0]}-{frequency[1]} Hz per stage separately.')
        
        ### 1. First check the directories
        # a. Check for output folder, if doesn't exist, create
        if not path.exists(self.out_dir):
            mkdir(self.out_dir)
        outdir = f'{self.out_dir}/powerspec_{model}'
        if path.exists(outdir):
            logger.debug(f"Output directory: {outdir} exists")
        else:
            logger.debug(f"Creating output directory: {outdir}")
            mkdir(outdir)
        
        # b. Get subject IDs
        subs = self.subs
        if isinstance(subs, list):
            None
        elif subs == 'all':
                subs = next(walk(self.xml_dir))[1]
        else:
            logger.error("'subs' must either be an array of participant ids or = 'all' ")
            
        # c. Take a look through directories to get sessions
        subs.sort()
        sessions = {}
        for s, sub in enumerate(subs):
            if self.sessions == 'all':
                sessions[sub] = next(walk(f'{self.xml_dir}/{sub}'))[1]
        sessions = list(set([y for x in sessions.values() for y in x]))           
        sessions.sort() 
  
        # 4. Begin data extraction
        for c, ch in enumerate(chan):
            logger.debug(f'Creating a summary dataset for channel: {ch}')
            
            # a. Create column names (append chan and ses names)
            for v, ses in enumerate(sessions):
                if general_opts['freq_full']:
                    sesvar = []
                    for pair in product(variables_full, [ses]):
                        sesvar.append('_'.join(pair))
                    columns_f = []
                    for pair in product(sesvar, [ch]):
                        columns_f.append('_'.join(pair))
                if general_opts['freq_band']:
                    sesvar = []
                    for pair in product(variables_band, [ses]):
                        sesvar.append('_'.join(pair))
                    columns_b = []
                    for pair in product(sesvar, [ch]):
                        columns_b.append('_'.join(pair))
                        
                # b. Extract data based on cycle and stage setup
                if model == 'whole_night':
                    stagename = '-'.join(self.stage)
                    if general_opts['freq_full']:
                        st_columns = [x + f'_{stagename}' for x in columns_f]
                        freq_full = DataFrame(index=subs, columns=st_columns, dtype=float)
                    if general_opts['freq_full']:    
                        st_columns = [x + f'_{stagename}' for x in columns_b]
                        freq_band = DataFrame(index=subs, columns=st_columns, dtype=float)
                    if not path.exists(f'{self.out_dir}/spectrum_{model}'):
                        mkdir(f'{self.out_dir}/spectrum_{model}')
                    for sub in subs: 
                        logger.debug(f'Extracting from {sub}, {ses}')
                        paramsfile = f'{self.xml_dir}/{sub}/{ses}/{sub}/{ses}/{sub}_{ses}_{ch}_{stagename}**{suffix}.csv'
                        if general_opts['freq_full']:
                            filename = f"{paramsfile.split('*')[0]}_freq_full_{paramsfile.split('*')[2]}"
                            if not path.exists(filename):
                                logger.warning(f"No (exported) spectrum data found for {sub}, {ses}. Has 'power_spectrum' by 'freq_full' been run for {model}? Skipping...")
                                flag += 1
                                continue
                            else:
                                try:
                                    df = read_csv(filename, skiprows=1)
                                    freq_full.loc[sub] = df.iloc[-1,idx_data_full]
                                except:
                                    extract_psa_error(logger)
                                    flag +=1
                                    continue
                        if general_opts['freq_band']:
                            filename = f"{paramsfile.split('*')[0]}_freq_band_{paramsfile.split('*')[2]}" 
                            if not path.exists(filename):
                                logger.warning(f"No (exported) spectrum data found for {sub}, {ses}. Has 'power_spectrum' by 'freq_band' been run for {model}? Skipping...")
                                flag += 1
                                continue
                            else:
                                try:
                                    df = read_csv(filename, skiprows=1)
                                    freq_band.loc[sub] = df.iloc[-1,idx_data_band]
                                except:
                                    extract_psa_error(logger)
                                    flag +=1
                                    continue
                    if general_opts['freq_full']:
                        freq_full.to_csv(f"{self.out_dir}/spectrum_{ses}_{ch}_{stagename}.csv")
                    if general_opts['freq_band']:
                        freq_band.to_csv(f"{self.out_dir}/spectrum_{ses}_{ch}_{stagename}.csv")
 
                elif model == 'stage*cycle':
                    for cyc in self.cycle_idx:
                        cycle = f'cycle{cyc}'
                        for st in self.stage:
                            if general_opts['freq_full']:
                                st_columns = [x + f'_{st}_{cycle}' for x in columns_f]
                                freq_full = DataFrame(index=subs, columns=st_columns, dtype=float)
                            if general_opts['freq_full']:    
                                st_columns = [x + f'_{st}_{cycle}' for x in columns_b]
                                freq_band = DataFrame(index=subs, columns=st_columns, dtype=float) 
                            for sub in subs: 
                                logger.debug(f'Extracting from {sub}, {ses}, stage {st}, {cycle}')
                                paramsfile = f'{self.xml_dir}/{sub}/{ses}/{sub}_{ses}_{ch}_{st}_cycle{cycle}**{suffix}.csv'
                                if general_opts['freq_full']:
                                    filename = f"{paramsfile.split('*')[0]}_freq_full_{paramsfile.split('*')[2]}"
                                    if not path.exists(filename):
                                        logger.warning(f"No (exported) spectrum data found for {sub}, {ses}. Has 'power_spectrum' by 'freq_full' been run for {model}? Skipping...")
                                        flag += 1
                                        continue
                                    else:
                                        try:
                                            df = read_csv(filename, skiprows=1)
                                            freq_full.loc[sub] = df.iloc[-1,idx_data_full]
                                        except:
                                            extract_psa_error(logger)
                                            flag +=1
                                            continue
                                if general_opts['freq_band']:
                                    filename = f"{paramsfile.split('*')[0]}_freq_band_{paramsfile.split('*')[2]}" 
                                    if not path.exists(filename):
                                        logger.warning(f"No (exported) spectrum data found for {sub}, {ses}. Has 'power_spectrum' by 'freq_band' been run for {model}? Skipping...")
                                        flag += 1
                                        continue
                                    else:
                                        try:
                                            df = read_csv(filename, skiprows=1)
                                            freq_band.loc[sub] = df.iloc[-1,idx_data_band]
                                        except:
                                            extract_psa_error(logger)
                                            flag +=1
                                            continue
                            if general_opts['freq_full']:
                                freq_full.to_csv(f"{self.out_dir}/pectrum_{ses}_{ch}_{st}_{cycle}.csv")
                            if general_opts['freq_band']:
                                freq_band.to_csv(f"{self.out_dir}/spectrum_{ses}_{ch}_{st}_{cycle}.csv")

                elif model == 'per_cycle':
                    for cyc in self.cycle_idx:
                        cycle = f'cycle{cyc}'
                        stagename = '-'.join(self.stage)
                        if general_opts['freq_full']:
                            st_columns = [x + f'_{stagename}_{cycle}' for x in columns_f]
                            freq_full = DataFrame(index=subs, columns=st_columns, dtype=float)
                        if general_opts['freq_full']:    
                            st_columns = [x + f'_{stagename}_{cycle}' for x in columns_b]
                            freq_band = DataFrame(index=subs, columns=st_columns, dtype=float) 
                        for sub in subs: 
                            logger.debug(f'Extracting from {sub}, {ses}, {cycle}')
                            paramsfile = f'{self.xml_dir}/{sub}/{ses}/{sub}_{ses}_{ch}_{stagename}_{cycle}**{suffix}.csv'
                            if general_opts['freq_full']:
                                filename = f"{paramsfile.split('*')[0]}_freq_full_{paramsfile.split('*')[2]}"
                                if not path.exists(filename):
                                    logger.warning(f"No (exported) spectrum data found for {sub}, {ses}. Has 'power_spectrum' by 'freq_full' been run for {model}? Skipping...")
                                    flag += 1
                                    continue
                                else:
                                    try:
                                        df = read_csv(filename, skiprows=1)
                                        freq_full.loc[sub] = df.iloc[-1,idx_data_full]
                                    except:
                                        extract_psa_error(logger)
                                        flag +=1
                                        continue
                            if general_opts['freq_band']:
                                filename = f"{paramsfile.split('*')[0]}_freq_band_{paramsfile.split('*')[2]}" 
                                if not path.exists(filename):
                                    logger.warning(f"No (exported) spectrum data found for {sub}, {ses}. Has 'power_spectrum' by 'freq_band' been run for {model}? Skipping...")
                                    flag += 1
                                    continue
                                else:
                                    try:
                                        df = read_csv(filename, skiprows=1)
                                        freq_band.loc[sub] = df.iloc[-1,idx_data_band]
                                    except:
                                        extract_psa_error(logger)
                                        flag +=1
                                        continue
                        if general_opts['freq_full']:
                            freq_full.to_csv(f"{self.out_dir}/spectrum_{ses}_{ch}_{stagename}_{cycle}.csv")
                        if general_opts['freq_band']:
                            freq_band.to_csv(f"{self.out_dir}/spectrum_{ses}_{ch}_{stagename}_{cycle}.csv")

                elif model == 'per_stage':
                    for st in self.stage:
                        if general_opts['freq_full']:
                            st_columns = [x + f'_{st}' for x in columns_f]
                            freq_full = DataFrame(index=subs, columns=st_columns, dtype=float)
                        if general_opts['freq_full']:    
                            st_columns = [x + f'_{st}' for x in columns_b]
                            freq_band = DataFrame(index=subs, columns=st_columns, dtype=float) 
                        for sub in subs: 
                            logger.debug(f'Extracting from {sub}, {ses}, stage {st}')
                            paramsfile = f'{self.xml_dir}/{sub}/{ses}/{sub}_{ses}_{ch}_{st}**{suffix}.csv'
                            if general_opts['freq_full']:
                                filename = f"{paramsfile.split('*')[0]}_freq_full_{paramsfile.split('*')[2]}"
                                if not path.exists(filename):
                                    logger.warning(f"No (exported) spectrum data found for {sub}, {ses}. Has 'power_spectrum' by 'freq_full' been run for {model}? Skipping...")
                                    flag += 1
                                    continue
                                else:
                                    try:
                                        df = read_csv(filename, skiprows=1)
                                        freq_full.loc[sub] = df.iloc[-1,idx_data_full]
                                    except:
                                        extract_psa_error(logger)
                                        flag +=1
                                        continue
                            if general_opts['freq_band']:
                                filename = f"{paramsfile.split('*')[0]}_freq_band_{paramsfile.split('*')[2]}" 
                                if not path.exists(filename):
                                    logger.warning(f"No (exported) spectrum data found for {sub}, {ses}. Has 'power_spectrum' by 'freq_band' been run for {model}? Skipping...")
                                    flag += 1
                                    continue
                                else:
                                    try:
                                        df = read_csv(filename, skiprows=1)
                                        freq_band.loc[sub] = df.iloc[-1,idx_data_band]
                                    except:
                                        extract_psa_error(logger)
                                        flag +=1
                                        continue
                        if general_opts['freq_full']:
                            freq_full.to_csv(f"{outdir}/spectrum_{ses}_{ch}_{st}.csv")
                        if general_opts['freq_band']:
                            freq_band.to_csv(f"{outdir}/spectrum_{ses}_{ch}_{st}.csv")
                            
         
        ### 3. Check completion status and print
        if flag == 0:
            logger.info('')
            logger.debug('Create powerspec dataset finished without ERROR.')  
        else:
            logger.info('')
            logger.warning('Create powerspec dataset finished with WARNINGS. See log for details.')
        return 
    
                        
# MISC FUNCTIONS        
def extract_psa_error(logger):
    logger.critical("Data extraction error: Check that all 'params' are written correctly.")
    logger.info('                      Check documentation for how event parameters need to be written:')
    logger.info('                      https://seapipe.readthedocs.io/en/latest/index.html')            
        
        
def save_freq(out, freq_bw, general_opts, outputfile, logger):

    low = where(out[0]['data'].axis['freq'][0]==freq_bw[0])[0][0]
    hi = where(out[0]['data'].axis['freq'][0]==freq_bw[1])[0][0]
    
    #[x[low:hi] for x['data'].axis['freq'][0] in out]
    for i, x in enumerate(out):
        out[i]['data'].axis['freq'][0] = x['data'].axis['freq'][0][low:hi]
        out[i]['data'].data[0] = reshape(out[i]['data'].data[0][0][low:hi], 
                                         (1, len(out[i]['data'].data[0][0][low:hi])))
    as_matrix = asarray([x['data'].data[0][0] for x in out])
    # out['data'].axis['freq'][0] = out['data'].axis['freq'][0][low:hi]
    # out['data'].data[0] = asarray([out['data'].data[0][0][low:hi]])
    #as_matrix = asarray([y for y in out[0]['data'].data[0]])
    desc = get_descriptives(as_matrix)
    filename = f"{outputfile.split('*')[0]}_freq_full_{outputfile.split('*')[2]}"
    logger.debug(f'Writing to {filename}')  
    export_freq(out, filename, desc=desc)
        
def save_bands(out, freq_bw, general_opts, frequency_opts, outputfile, logger):
    filename = f"{outputfile.split('*')[0]}_freq_band_{outputfile.split('*')[2]}"
    logger.debug(f'Writing to {filename}') 
    try:
        export_freq_band(out, frequency_opts['bands'], filename)
    except:
        logger.error('Cannot export power in user-defined frequency bands. Check bands. For info on how to define frequency bands, refer to documentation:')
        logger.info('https://seapipe.readthedocs.io/en/latest/index.html')

def save_plot(out, general_opts, outputfile, logger):
    
    ''' TO COMPLETE '''
    
    # fig = Figure()
    # FigureCanvas(fig)
    # idx_plot = 1
    # data = asarray([x['data'].data[0][0] for x in out])
    # seg_chan = data[0].axis['chan'][0]
    # sf = out[0]['data'].axis['freq'][0]
    # if general_opts['max_freq_plot']:
    #     idx_max = asarray(
    #             [abs(x - general_opts['max_freq_plot']) for x in sf]).argmin()
    # for ch in seg_chan:
    #     Sxx = data(chan=ch)[0]
    #     ax = fig.add_subplot(len([out]), len(seg_chan), idx_plot)
    #     ax.semilogy(sf[1:idx_max], Sxx[1:idx_max])
    #     ax.set_xlabel('Frequency (Hz)')
    #     idx_plot += 1 
    # filename = f"{outputfile.split('*')[0]}_powerspectrum.svg"
    # logger.debug(f'Saving figure to {filename}')
    # fig.savefig(f'{filename}')  
        
        