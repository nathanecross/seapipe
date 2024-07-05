# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 10:29:11 2021

@author: Nathan Cross
"""
from copy import deepcopy
from csv import reader
from fooof import FOOOF
from fooof.analysis import get_band_peak_fm
from glob import glob
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from numpy import (asarray, ceil, concatenate, empty, floor, mean, nan, 
                   ones, pi, sqrt, stack, sum, zeros)
from openpyxl import Workbook
from os import listdir, mkdir, path, walk
from pandas import DataFrame, ExcelWriter
from scipy.fftpack import next_fast_len
from wonambi import ChanFreq, Dataset
from wonambi.attr import Annotations
from wonambi.trans import (fetch, frequency, get_descriptives, export_freq, 
                           export_freq_band)
from ..utils.misc import bandpass_mne, laplacian_mne, notch_mne, notch_mne2
from ..utils.logs import create_logger
from ..utils.load import load_channels, rename_channels


def default_general_opts():
    general_opts = {'freq_full':True, 
                    'freq_band':True, 
                    'freq_plot':False,
                    'max_freq_plot':35,
                    'suffix':None,
                    'chan_grp_name':'eeg'}
    return general_opts

def default_frequency_opts():
    bands = [(0.25, 1.25),  #SO
             (0.5, 4),      #delta
             (4, 7.75),     #theta
             (8, 11),       #alpha
             (11.25, 16),   #sigma
             (16.25, 19),   #low beta
             (19.25, 35)]   #high beta
    frequency_opts = {'bands': bands,
                      'output':'spectraldensity', 
                      'scaling':'power', 
                      'sides':'one', 
                      'taper':'hann', 
                      'detrend':'linear', 
                      'n_fft': None, 
                      'fast_fft': False, 
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
                   'bandpass': False,
                   'highpass': None,
                   'lowpass': 50}
    return filter_opts

def default_epoch_opts():
    epoch_opts = {'epoch': None,
                  'reject_epoch': True,
                  'reject_artf': True,
                  'min_dur':0,
                  'epoch_dur': 30,
                  'epoch_overlap': 0.5,
                  'epoch_step': None}
    
    return epoch_opts

def default_event_opts():
    event_opts = {'evt_type': None,
                  'event_chan': None}
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
                 grp_name, stage, frequency = (11,16), rater = None, 
                 subs = 'all', sessions = 'all'):
        
        self.rec_dir = rec_dir
        self.xml_dir = xml_dir
        self.out_dir = out_dir
        self.log_dir = log_dir
        
        self.chan = chan
        self.ref_chan = ref_chan
        self.grp_name = grp_name
        self.stage = stage
        self.frequency = frequency
        self.rater = rater
        
        self.subs = subs
        self.sessions = sessions
        
        self.tracking = {}

    
    def fooof_it(self, general_opts, frequency_opts, filter_opts, 
                 epoch_opts, event_opts, fooof_opts, chan = None, rater = None, 
                 grp_name = 'eeg', cat = (1,1,1,1), cycle_idx = None,
                 filetype = '.edf', suffix = 'fooof'):
        

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
        tracking = self.tracking
        flag = 0
        
        logger.info('')
        
        ### 0.b. Set up organisation of export
        if cat[0] + cat[1] == 2:
            model = 'whole_night'
            logger.debug(f'Parameterizing power spectrum in range {self.frequency[0]}-{self.frequency[1]} Hz for the whole night.')
        elif cat[0] + cat[1] == 0:
            model = 'stage*cycle'
            logger.debug(f'Parameterizing power spectrum in range {self.frequency[0]}-{self.frequency[1]} Hz per stage and cycle separately.')
        elif cat[0] == 0:
            model = 'per_cycle'
            logger.debug(f'Parameterizing power spectrum in range {self.frequency[0]}-{self.frequency[1]} Hz per cycle separately.')
        elif cat[1] == 0:
            model = 'per_stage'  
            logger.debug(f'Parameterizing power spectrum in range {self.frequency[0]}-{self.frequency[1]} Hz per stage separately.')
        if 'cycle' in model and cycle_idx == None:
            logger.info('')
            logger.critical("To run cycles separately (i.e. cat[0] = 0), cycle_idx cannot be 'None'")
            return
        cat = tuple((cat[0],cat[1],1,1)) # force concatenation of discontinuous & events
        logger.info('')
        logger.debug(r"""
                     
                                  |
                                  |
                                  |  .
                                  |  \
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
            logger.critical("'subs' must either be an array of Participant IDs or = 'all' ")
            return
        subs.sort()
        
        # 3.a. Begin loop through participants
        for p, sub in enumerate(subs):
            tracking[f'{sub}'] = {}
            
            # b. Begin loop through sessions
            sessions = self.sessions
            if sessions == 'all':
                sessions = listdir(f'{self.rec_dir}/{sub}')
                sessions = [x for x in sessions if not '.' in x] 

            for v, ses in enumerate(sessions):
                logger.info('')
                logger.debug(f'Commencing {sub}, {ses}')
                tracking[f'{sub}'][f'{ses}'] = {'Spec_peaks':{}} 
                
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
                annot = Annotations(xdir + xml_file[0], rater_name=rater)
                 
                ### get cycles
                if cycle_idx is not None:
                    all_cycles = annot.get_cycles()
                    cycle = [all_cycles[i - 1] for i in cycle_idx if i <= len(all_cycles)]
                else:
                    cycle = None
                    
                ### if event channel only, specify event channels
                # 4.d. Channel setup
                if not chan:
                    flag, chanset = load_channels(sub, ses, self.chan, 
                                                  self.ref_chan, flag, logger)
                    if not chanset:
                        break
                
                    newchans = rename_channels(sub, ses, self.chan, logger)

                
                for c, ch in enumerate(chanset):
                    logger.debug(f"Reading data for {ch}:{'/'.join(chanset[ch])}")
                    segments = fetch(dset, annot, cat = cat, evt_type = None, 
                                     stage = self.stage, cycle=cycle,  
                                     epoch = epoch_opts['epoch'], 
                                     epoch_dur = epoch_opts['epoch_dur'], 
                                     epoch_overlap = epoch_opts['epoch_overlap'], 
                                     epoch_step = epoch_opts['epoch_step'], 
                                     reject_epoch = epoch_opts['reject_epoch'], 
                                     reject_artf = epoch_opts['reject_artf'],
                                     min_dur = epoch_opts['min_dur'])
                    
                    if len(segments)==0:
                        logger.warning(f"No valid data found for {ch}:{'/'.join(chanset[ch])}.")
                        continue
                    
                    # 5.b Rename channel for output file (if required)
                    if newchans:
                        fnamechan = newchans[ch]
                    else:
                        fnamechan = ch
                    
                    segments.read_data(chan = ch, ref_chan = chanset[ch])
                    
                    for sg, seg in enumerate(segments):
                        logger.debug(f'Analysing segment {sg+1} of {len(segments)}')
                        out = dict(seg)
                        data = seg['data']
                        timeline = data.axis['time'][0]
                        out['start'] = timeline[0]
                        out['end'] = timeline[-1]
                        out['duration'] = len(timeline) / data.s_freq
                        
                        if frequency_opts['fast_fft']:
                            n_fft = next_fast_len(data.number_of('time')[0])
                        else:
                            n_fft = frequency_opts['n_fft']
                    
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
                        
                        ### WHOLE NIGHT ###
                        if model == 'whole_night':
                            stagename = '-'.join(self.stage)
                            outputfile = f'{outpath}/{sub}_{ses}_{fnamechan}_{stagename}_specparams_{suffix}.csv'
                        elif model == 'stage*cycle':    
                            outputfile = f'{outpath}/{sub}_{ses}_{fnamechan}_{self.stage[sg]}_cycle{cycle_idx[sg]}_specparams_{suffix}.csv'
                        elif model == 'per_stage':
                            outputfile = f'{outpath}/{sub}_{ses}_{fnamechan}_{self.stage[sg]}_specparams_{suffix}.csv'
                        elif model == 'per_cycle':
                            stagename = '-'.join(self.stage)
                            outputfile = f'{outpath}/{sub}_{ses}_{fnamechan}_{stagename}_cycle{cycle_idx[sg]}_specparams_{suffix}.csv'

                        logger.debug(f'Saving {outputfile}')
                        df = DataFrame(data = one_record, 
                                       columns = (seg_info + band_hdr + band_pk_params_hdr
                                                  + band_ap_params_hdr))
                        df.to_csv(outputfile)
        
        return

    def powerspec_it(self, cat, stage, chan, ref_chan, 
                     general_opts, frequency_opts, rater = None, cycle_idx = None, 
                     subs = 'all', sessions = 'all', filter_opts = None, 
                     epoch_opts = None, event_opts = None, norm = None, 
                     norm_opts = None, fooof_it = False,
                     fooof_opts = None):
        
        print(r"""    Calculating power spectrum...
              
                  |
                  | /\ 
                  |/  \
              uV2 |    ^-_
                  |       ^-___-__
                  |________________
                         (Hz)
              
              """)
            
        if path.exists(self.out_dir):
                print(self.out_dir + " already exists")
        else:
            mkdir(self.out_dir)
        
        
        suffix = general_opts['suffix']
            
        if fooof_it:
            if fooof_opts is None:
                print('Error: Options not set for FOOOF')
            else:
                # Define model parameters
                fm = FOOOF(fooof_opts['peak_width_limits'], fooof_opts['max_n_peaks'], 
                           fooof_opts['min_peak_amplitude'], fooof_opts['peak_threshold'])
            
                def gaussian_integral(a, c):
                    """ Returns definite integral of a gaussian function with height a and
                    standard deviation c."""
                    return sqrt(2) * a * abs(c) * sqrt(pi)
            
            # Prepare Bands
            if fooof_opts['bands_fooof'] is not None:
                bands = fooof_opts['bands_fooof']
            else:
                stp = min(fooof_opts['peak_width_limits'])
                low = int(floor(fooof_opts['freq_range'][0]))
                hi = int(ceil(fooof_opts['freq_range'][1]))
                bands = [(x,x + stp) for x in range(low,hi)]
            
        
        
        ### loop through records
        if isinstance(subs, list):
            None
        elif subs == 'all':
                subs = listdir(self.rec_dir)
                subs = [ p for p in subs if not '.' in p]
        else:
            print("ERROR: 'part' must either be an array of subject ids or = 'all' ")  
    
        if norm:
            if norm not in ('integral', 'baseline'):
                exit('Invalid value for norm: ' + str(norm))
                
        for i, p in enumerate(subs):
            # loop through visits
            if sessions == 'all':
                sessions = listdir(self.rec_dir + '/' + p)
                sessions = [x for x in sessions if not '.' in x]    
            
            for v, ses in enumerate(sessions):
                ## Define files
                rdir = self.rec_dir + p + '/' + ses + '/'
                xdir = self.xml_dir + p + '/' + ses + '/'
                edf_file = [x for x in listdir(rdir) if x.endswith('.edf') or x.endswith('.rec') 
                            or x.endswith('.eeg') or x.endswith('.set')]
                xml_file = [x for x in listdir(xdir) if x.endswith('.xml')]            
                
                ### Define output path
                if not path.exists(self.out_dir):
                    mkdir(self.out_dir)
                if not path.exists(self.out_dir + p ):
                    mkdir(self.out_dir + p)
                if not path.exists(self.out_dir + p + '/' + ses):
                    mkdir(self.out_dir + p + '/' + ses)
                outpath = self.out_dir + p + '/' + ses + '/'
                
                
                ## Now import data
                dset = Dataset(rdir + edf_file[0])
                annot = Annotations(xdir + xml_file[0], rater_name=rater)
                 
                ### get cycles
                if cycle_idx is not None:
                    all_cycles = annot.get_cycles()
                    cycle = [all_cycles[i - 1] for i in cycle_idx if i <= len(all_cycles)]
                else:
                    cycle = None
                    
                ### if event channel only, specify event channels
                chan_full = None
                if event_opts['event_chan']:
                    chan_full = [i + ' (' + general_opts['chan_grp_name'] + ')' for i in event_opts['event_chan']]
                
                ### select and read data
                print('Reading data for ' + p + ', Visit: ' + ses)
                segments = fetch(dset, annot, cat=cat, evt_type=event_opts['evt_type'], 
                                 stage=stage, cycle=cycle, chan_full=chan_full, 
                                 epoch=epoch_opts['epoch'], 
                                 epoch_dur=epoch_opts['epoch_dur'], 
                                 epoch_overlap=epoch_opts['epoch_overlap'], 
                                 epoch_step=epoch_opts['epoch_step'], 
                                 reject_epoch=general_opts['reject_epoch'], 
                                 reject_artf=general_opts['reject_artf'],
                                 min_dur=general_opts['min_dur'])
                
                if not segments:
                    print('No valid data found. Skipping to next record.')
                    continue
                
                # BASELINE NORMALISATION
                if norm == 'baseline':
                    norm_seg = fetch(dset, annot, cat=norm_opts['norm_cat'], 
                                     evt_type=norm_opts['norm_evt_type'],
                                     stage=norm_opts['norm_stage'], 
                                     epoch=norm_opts['norm_epoch'])
                    
                    if not norm_seg:
                        print('No valid normalization data found. '
                              'Skipping to next record.')
                        continue
                    
                    if filter_opts['laplacian']:
                        norm_seg.read_data(filter_opts['lapchan'], ref_chan) 
                    else:
                        norm_seg.read_data(chan, ref_chan)            
                    all_nSxx = []
                    
                    for seg in norm_seg:
                        
                        normdata = seg['data']
                        
                        if filter_opts['laplacian']:
                            normdata.data[0] = laplacian_mne(normdata, 
                                                     filter_opts['oREF'], 
                                                     channel=chan, 
                                                     ref_chan=ref_chan, 
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
                
                
                if filter_opts['laplacian']:
                    segments.read_data(filter_opts['lapchan'], ref_chan)
                else:
                    segments.read_data(chan, ref_chan)
                xfreq = []
                
                for sg, seg in enumerate(segments):
                    print(f'Analysing segment {sg} of {len(segments)}')
                    out = dict(seg)
                    data = seg['data']
                    timeline = data.axis['time'][0]
                    out['start'] = timeline[0]
                    out['end'] = timeline[-1]
                    out['duration'] = len(timeline) / data.s_freq
                    
                    if frequency_opts['fast_fft']:
                        n_fft = next_fast_len(data.number_of('time')[0])
                    else:
                        n_fft = frequency_opts['n_fft']
                    
                    if filter_opts['laplacian']:
                        selectchans = filter_opts['lapchan']
                    else:
                        selectchans = chan
                    
                    if filter_opts['notch']:
                        print('Applying notch filtering.')
                        data.data[0] = notch_mne(data, oREF=filter_opts['oREF'], 
                                                    channel=selectchans, 
                                                    freq=filter_opts['notch_freq'],
                                                    rename=filter_opts['laplacian_rename'],
                                                    renames=filter_opts['renames'],
                                                    montage=filter_opts['montage'])
                        
                    if filter_opts['notch_harmonics']: 
                        print('Applying notch harmonics filtering.')
                        print(f'{selectchans}')
                        data.data[0] = notch_mne2(data, oREF=filter_opts['oREF'], 
                                                  channel=selectchans, 
                                                  rename=filter_opts['laplacian_rename'],
                                                  renames=filter_opts['renames'],
                                                  montage=filter_opts['montage'])    
                    
                    if filter_opts['bandpass']:
                        print('Applying bandpass filtering.')
                        data.data[0] = bandpass_mne(data, oREF=filter_opts['oREF'], 
                                                  channel=selectchans,
                                                  highpass=filter_opts['highpass'], 
                                                  lowpass=filter_opts['lowpass'], 
                                                  rename=filter_opts['laplacian_rename'],
                                                  renames=filter_opts['renames'],
                                                  montage=filter_opts['montage'])
                    
                    if filter_opts['laplacian']:
                        print('Applying Laplacian filtering.')
                        data.data[0] = laplacian_mne(data, 
                                             filter_opts['oREF'], 
                                             channel=chan, 
                                             ref_chan=ref_chan, 
                                             laplacian_rename=filter_opts['laplacian_rename'], 
                                             renames=filter_opts['renames'],
                                             montage=filter_opts['montage'])
                        data.axis['chan'][0] = asarray([x for x in chan])
                        selectchans = chan
                    
                    
                    ### frequency transformation
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
                    
                    
                    
                    if norm:
            
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
                    
                    
                    if fooof_it:
                        
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
                        
    
                                
                        fooof_powers = zeros((len(chan), len(bands)))
                        fooof_ap_params = zeros((len(chan), len(bands), 2))
                        fooof_pk_params = ones((len(chan), len(bands), 9)) * nan
                        
                        for i in range(len(chan)):
                            fm.fit(freqs, Fooofxx.data[0][i], fooof_opts['freq_range'])
                            
                            for j, band in enumerate(bands):
                                fp = get_band_peak_fm(fm, band, fooof_opts['select_highest'],
                                                      threshold=fooof_opts['thresh_select'],
                                                      thresh_param=fooof_opts['thresh_param'])
                                if fp.ndim == 1:
                                    fooof_powers[i, j] = gaussian_integral(fp[1], fp[2])
                                else:
                                    pwr = asarray([gaussian_integral(fp[x, 1], fp[x, 2]) \
                                                   for x in range(fp.shape[0])]).sum()
                                    fooof_powers[i, j] = pwr
                                    
                                # get fooof aperiodic parameters
                                fooof_ap_params[i, j, :] = fm.aperiodic_params_
                                
                                # get fooof peak parameters
                                fp = get_band_peak_fm(fm, band, False,
                                                      threshold=fooof_opts['thresh_select'],
                                                      thresh_param=fooof_opts['thresh_param'])
                                if fp.ndim == 1:
                                    fooof_pk_params[i, j, :3] = fp
                                else:
                                    n_peaks = min(fp.shape[0], 3)
                                    fooof_pk_params[i, j, :n_peaks * 3] = fp[:n_peaks, 
                                                                          :].ravel()
                        
                        out['fooof_powers'] = fooof_powers
                        out['fooof_ap_params'] = fooof_ap_params
                        out['fooof_pk_params'] = fooof_pk_params
                                
                    xfreq.append(out)
                    
                if general_opts['freq_full']:                        
                    if len(xfreq) == 1:
                        desc = None
                    else:
                        as_matrix = asarray([y for x in xfreq for y in x['data']()[0]])
                        desc = get_descriptives(as_matrix)
                    
                    filename = outpath + p + '_' + ses + f'_freq_full{suffix}.csv'
                    print('Writing to ' + filename)  
                    export_freq(xfreq, filename, desc=desc)
                    
                if general_opts['freq_band']:
                    filename = outpath + p + '_' + ses + f'_freq_band{suffix}.csv'
                    print('Writing to ' + filename)  
                    export_freq_band(xfreq, frequency_opts['bands'], filename)
                    
                if general_opts['freq_plot']:
                    fig = Figure()
                    FigureCanvas(fig)
                    idx_plot = 1
                    
                    for seg in xfreq:
                        data = seg['data']
                        seg_chan = data.axis['chan'][0]
                        sf = data.axis['freq'][0]
                        if general_opts['max_freq_plot']:
                            idx_max = asarray(
                                    [abs(x - general_opts['max_freq_plot']) for x in sf]).argmin()
                        for ch in seg_chan:
                            Sxx = data(chan=ch)[0]
                            ax = fig.add_subplot(len(xfreq), len(seg_chan), idx_plot)
                            ax.semilogy(sf[1:idx_max], Sxx[1:idx_max])
                            ax.set_xlabel('Frequency (Hz)')
                            
                            idx_plot += 1
                        
                    print('Saving figure to ' + outpath + p + '.png')
                    fig.savefig(outpath + p + suffix + '.png')
                    
                if fooof_it:
                    seg_info = ['Start time', 'End time', 'Duration', 'Stitches', 
                        'Stage', 'Cycle', 'Event type', 'Channel']
                    #band_hdr = [str(b1) + '-' + str(b2) for b1, b2 in bands_fooof]
                    band_hdr = [f'{b1}-{b2} Hz' for b1, b2 in bands]
                    pk_params_hdr = ['peak1_CF', 'peak1_PW', 'peak1_BW', 
                                  'peak2_CF', 'peak2_PW', 'peak2_BW', 
                                  'peak3_CF', 'peak3_PW', 'peak3_BW', ]
                    ap_params_hdr = ['Offset', 'Exponent']
                    band_pk_params_hdr = ['_'.join((b, p)) for b in band_hdr 
                                          for p in pk_params_hdr]
                    band_ap_params_hdr = ['_'.join((b, p)) for b in band_hdr 
                                          for p in ap_params_hdr]
                    one_record = zeros((len(xfreq) * len(chan), 
                                        (len(seg_info) + len(bands) + 
                                        len(band_pk_params_hdr) + 
                                        len(band_ap_params_hdr))),
                                        dtype='O')
                    for i, seg in enumerate(xfreq):
                        for j, ch in enumerate(chan):
                            cyc = None
                            if seg['cycle'] is not None:
                                cyc = seg['cycle'][2]
                            
                            one_record[i * len(chan) + j, :] = concatenate((asarray([
                                seg['start'],
                                seg['end'],
                                seg['duration'], 
                                seg['n_stitch'], # number of concatenated segments minus 1
                                seg['stage'],
                                cyc,
                                seg['name'], # event type
                                ch,
                                ]),
                                seg['fooof_powers'][j, :],
                                seg['fooof_pk_params'][j, ...].ravel(),
                                seg['fooof_ap_params'][j, ...].ravel(),
                                ))
                            
                    outpath_fooof = outpath + p + '_' + ses + f'_fooofbands{suffix}.csv'
                    print(f'Saving {outpath_fooof}')
                    df = DataFrame(data=one_record, 
                                   columns=(seg_info + band_hdr + band_pk_params_hdr
                                            + band_ap_params_hdr))
                    df.to_csv(outpath_fooof)
        
        print("Job's done.")
        
        
        
    def powerspec_summary_bands(root_dir, out_dir, part, visit, chan, stage, col_headers, 
                          excel_output_file):
    
        ## Script will create a file with multiple tab: once per STAGE*CHANNEL
        
        """ SCRIPT """
        csvs = []
        which_files = '*_freq_band.csv' # Name extension of output files to search for
        tabs = [' '.join((c,s)) for c in chan for s in stage] # Create matrices for each stage*channel combo
        
        # Loop through participants and visits to extract individual CSV output files
        if isinstance(part, list):
                None
        elif part == 'all':
                part = listdir(root_dir)
                part = [x for x in part if not('.' in x)]
        else:
            print("ERROR: 'part' must either be an array of subject ids or = 'all' ")      
                
        for i, p in enumerate(part):
                print(p)
                # loop through visits
                if visit == 'all':
                    visit = listdir(root_dir + '/' + p)
                    visit = [x for x in visit if not('.' in x)]
                
                for v, vis in enumerate(visit):
                    ## Define files
                    rdir = root_dir + p + '/' + vis + '/'
                    
                    try:
                        csvs.append(glob(f'{rdir}/{which_files}')[0]) #glob= va chercher tous les path du filename
                    except Exception as e:
                        print(f'No files found for {p}, visit {vis}')
                    
        ids = [x + '_' + y for x in part for y in visit] 
        idx_data_col = list(range(9,9+len(col_headers)))
        data = ones((len(chan) * len(stage), len(ids), len(idx_data_col))) * nan #create matrice - *nan to make sure that missing value will be blank
            
        for i, one_csv in enumerate(csvs): #on va passer dans la liste 1 a 1 (commence a compter a 0)    
            with open(one_csv, newline='') as f: #with open = to open a text/csv file
                csv_reader = reader(f, delimiter=',') #f = file, delimiter = what separate value in xcl (;, ou, ou .)
                
                for row in range(2): #skip first line
                    row = next(csv_reader) #next = move forward
                    
                if col_headers is None:
                    col_headers = [row[x] for x in idx_data_col] #if None in col_header then will use the name of the csv
                
                for row in csv_reader:
                    try: #try below, but if valueerror continue to except
                        int(row[0]) #find integer (chiffre)
                        idx_tab = tabs.index(' '.join((row[8], row[5]))) #join column 8 (channel) and 5(stage) to know which tab to go to
                        
                        for j, idx_col in enumerate(idx_data_col):
                            data[idx_tab, i, j] = row[idx_col] #create data, 1 column at a time and write it in correct tab
                                        
                    except ValueError:
                        continue
        
        wb = Workbook()
        wb.save(excel_output_file)
                       
        with ExcelWriter(excel_output_file, engine="openpyxl", mode='a') as writer: #create xcl
        
            for i, tab_label in enumerate(tabs):
                df = DataFrame(data[i, ...], index=ids, columns=col_headers) #create dataframe and write data in it
                df.to_excel(writer, sheet_name=tab_label) #method of df to create xcl
            
        print(f'Saved to {excel_output_file}')
        
        
        
    def powerspec_summary_full(root_dir, out_dir, part, visit, chan, stage, lowpass,
                          excel_output_file):
    
        ## Script will create a file with multiple tab: once per STAGE*CHANNEL
        
        """ SCRIPT """
        csvs = []
        which_files = '*_freq_full.csv' # Name extension of output files to search for
        tabs = [' '.join((c,s)) for c in chan for s in stage] # Create matrices for each stage*channel combo
        
        # Loop through participants and visits to extract individual CSV output files
        if isinstance(part, list):
                None
        elif part == 'all':
                part = listdir(root_dir)
                part = [x for x in part if not('.' in x)]
        else:
            print("ERROR: 'part' must either be an array of subject ids or = 'all' ")      
                
        for i, p in enumerate(part):
                print(p)
                # loop through visits
                if visit == 'all':
                    visit = listdir(root_dir + '/' + p)
                    visit = [x for x in visit if not('.' in x)]
                
                for v, vis in enumerate(visit):
                    ## Define files
                    rdir = root_dir + p + '/' + vis + '/'
                    
                    try:
                        csvs.append(glob(f'{rdir}/{which_files}')[0]) #glob= va chercher tous les path du filename
                    except Exception as e:
                        print(f'No files found for {p}, visit {vis}')
                    
        ids = [x + '_' + y for x in part for y in visit] 
        idx_data_col = list(range(9,9+lowpass*4))
        data = ones((len(chan) * len(stage), len(ids), len(idx_data_col))) * nan #create matrice - *nan to make sure that missing value will be blank
            
        for i, one_csv in enumerate(csvs): #on va passer dans la liste 1 a 1 (commence a compter a 0)    
            with open(one_csv, newline='') as f: #with open = to open a text/csv file
                csv_reader = reader(f, delimiter=',') #f = file, delimiter = what separate value in xcl (;, ou, ou .)
                
                for row in range(2): #skip first line
                    row = next(csv_reader) #next = move forward
                    
                
                col_headers = [row[x] for x in idx_data_col] #if None in col_header then will use the name of the csv
                
                for row in csv_reader:
                    try: #try below, but if valueerror continue to except
                        int(row[0]) #find integer (chiffre)
                        idx_tab = tabs.index(' '.join((row[8], row[5]))) #join column 8 (channel) and 5(stage) to know which tab to go to
                        
                        for j, idx_col in enumerate(idx_data_col):
                            data[idx_tab, i, j] = row[idx_col] #create data, 1 column at a time and write it in correct tab
                                        
                    except ValueError:
                        continue
        
        wb = Workbook()
        wb.save(excel_output_file)
                       
        with ExcelWriter(excel_output_file, engine="openpyxl", mode='a') as writer: #create xcl
        
            for i, tab_label in enumerate(tabs):
                df = DataFrame(data[i, ...], index=ids, columns=col_headers) #create dataframe and write data in it
                df.to_excel(writer, sheet_name=tab_label) #method of df to create xcl
            
        print(f'Saved to {excel_output_file}')