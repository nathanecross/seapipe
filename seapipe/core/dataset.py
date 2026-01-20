#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 12:07:36 2023

@author: nathancross
"""
from datetime import datetime, date
from os import listdir, mkdir, path, remove, walk
from pathlib import Path
import numpy as np
import csv
from pandas import DataFrame
from seapipe.events.fish import FISH
from seapipe.events.whales import cluster_peaks, whales
from seapipe.spectrum.psa import (Spectrum, default_epoch_opts, default_event_opts,
                     default_fooof_opts, default_filter_opts, default_frequency_opts, 
                     default_general_opts,default_norm_opts)
from seapipe.utils.squid import SQUID, gather_qc_reports
from seapipe.stats import sleepstats
from seapipe.utils.audit import (check_dataset, extract_channels, make_bids,
                        track_processing)
from seapipe.utils.logs import create_logger, create_logger_outfile
from seapipe.utils.load import (check_chans, read_tracking_sheet, 
                                select_input_dirs, select_output_dirs, load_stages)
from seapipe.utils.misc import adap_bands_setup

## TO DO:
#   - adapt load channels to be flexible for non-equivalent refsets and chansets
#   - related to the above, fix logging when using adapted bands but chan or ref_chan
#               are set to something other than None - e.g. psa.py, line 367
#   - add in log for detection whether auto, fixed or adapted bands was run
#   - add logging to save to output file (not implemented for all functions)
#   - update adapted bands in tracking.tsv
#   - fix discrepency between track WARNINGS and output in dataframe 
#   - update initial tracking to include spindles, slow_oscillation, cfc, power_spectrum
#   - update export sleepstats to export by stage/cycle separately
#   - possibility of cycle_idx = 'all'
#   - enable macro_dataset per sleep cycle
#   - enable downsampling of data

## FOR DOCUMENTATION:
#   - Clearly describe how chanset & refset works, ie. chanset = per chan, refset = ALL chans

class pipeline:
        
    """Contains specific information and allows the application of methods of 
    analysis, associated with a dataset. 

    Parameters
    ----------
    indir : str 
        name of the root level directory containing the BIDS organised data
        
    outfile : bool / str
        whether to save log of dataset audit to file. If False (default) - does
        not save any log file. If True - saves a log under the filepath 
        /derivatives/seapipe/audit/audit.csv. Else if a string to a filepath, 
        it will save the log under the filepath indicated in the string.
        
    filetype : str
        extension of file to search for. Default is '.edf' - but according to 
        BIDS convention other filetypes can be '.vhdr', '.vmrk', '.eeg' or '.set'

    Attributes
    ----------
    rootpath : str
        name of the root level directory
    datapath : str
        name of the directory containing the raw data (recordings and annotations)
    outpath : str
        name of the directory containing the output (analysed) data

    """
        
    def __init__(self, indir, tracking = False, outfile = False, 
                       filetype = '.edf'):
        
        self.rootpath = indir
        if path.exists(indir + '/DATA'):
            self.datapath = indir + '/DATA'
        else:
            self.datapath = indir + '/sourcedata'
        self.outpath = indir + '/derivatives'
        if not path.exists(self.outpath):
            mkdir(self.outpath)
        self.outfile = outfile
        if not path.exists(f'{self.outpath}/audit'):
            mkdir(f'{self.outpath}/audit')
        self.log_dir = self.outpath + '/audit/logs/'
        if not path.exists(self.log_dir):
            mkdir(self.log_dir)
        self.tracking = {}
        self.audit_init = check_dataset(self.rootpath, self.datapath, 
                                        self.outfile, filetype, tracking)
        
        if tracking:
            self.track(subs = 'all', ses = 'all', 
                   step = ['staging','spindle','slowwave','pac','sync','psa'],
                   show = False, log = False)
        
    
    #--------------------------------------------------------------------------
    '''
    MISCELLANEOUS FUNCTIONS
    
    audit -> Audits dataset structure for compatibility with seapipe analysis.
    
    list_dataset ->  Intended to walk from root directory through participant 
                        folders and list all participants and their files.
    
    track -> Tracks what seapipe processing or functions have already been applied
                to a dataset, with information on which channels and parameters 
                have been used.
                
    make_bids (beta) -> Transforms data from (some) data structures into the 
                            correct BIDS format compatible with use for seapipe.
                            
    extract_channels -> Extracts and lists which channels exist in the dataset. 

    load_stages -> Extracts stages from the BIDS formatted dataset, in which
                        staging has been listed in a file *acq-PSGScoring_events.tsv
    
    '''    
        
        
    def audit(self, outfile = False, tracking = False, filetype = '.edf'):
        
        ''' Audits the dataset for BIDS compatibility.
            Includes option to save the audit to an output file.
        '''
        
        # Create audit directory
        out_dir = f'{self.outpath}/audit'
        if not path.exists(out_dir):
            mkdir(out_dir)
            
        if not outfile and not self.outfile:
            logger = create_logger("Audit")
            logger.propagate = False
            self.audit_update = check_dataset(self.rootpath, self.datapath,
                                              outfile, filetype, tracking,  
                                              logger)
        else:
            if not outfile:
                outfile = self.outfile
            out = f'{out_dir}/{outfile}'
            if path.exists(out):
                remove(out)
            logger = create_logger_outfile(outfile, name = 'Audit')
            logger.propagate = False
            self.audit_update = check_dataset(self.rootpath, self.datapath,
                                              outfile, filetype, tracking,  
                                              logger)
            
        logger.info('')
        logger.info(self.audit_update)
        
        
    def list_dataset(self, outfile=False): 
        
        """Prints out all the files inside the directory <in_dir> along with the
        directories 1 and 2 levels above containing the files. You can specify 
        an optional output filename that will contain the printout.
        """

        if not outfile and not self.outfile:
            logger = create_logger('Audit')  
        else:
            if not outfile:
                outfile = self.outfile
            out_dir = f'{self.outpath}/audit'
            if not path.exists(out_dir):
                mkdir(out_dir)
            out = f'{out_dir}/{outfile}'
            if path.exists(out):
                remove(out)
            logger = create_logger_outfile(out, name='Audit')

        logger.propagate = False
        
        logger.info("")
        logger.info("")
        for dirPath, dirNames, fileNames in walk(self.datapath):
            try:
                fileNames.remove('.DS_Store')
            except(ValueError):
                pass
            
            if fileNames or dirPath.split('/')[-1]=='eeg':
                dir1 = dirPath.split('/')[-3]
                dir2 = dirPath.split('/')[-2]
                dir3 = dirPath.split('/')[-1]
                logger.info(f"Directory: {dir1}/{dir2}/{dir3}")
                logger.info(f"Files; {fileNames}")
                logger.info('-' * 10)

    
    def track(self, subs = 'all', ses = 'all', step = None, chan = None, 
                    stage = None, outfile = False, show = True, log = True):
        
        ## Set up logging
        logger = create_logger('Tracking')
        logger.info('')
        
        ## Set tracking variable
        if self.tracking:
            tracking = self.tracking
        else:
            tracking = {}
        
        ## Track sessions  
        if not isinstance(subs, list) and subs == 'all':
            try:
                subs = [x for x in listdir(self.datapath) if '.' not in x]
            except:
                logger.critical(f'{self.datapath} does not exist - cannot ascertain '
                             'details of dataset.')
                return
        elif not isinstance(subs, list):
            
            subs = read_tracking_sheet(self.rootpath, logger)
            subs = subs['sub'].drop_duplicates().tolist()
        subs.sort()
        
        # Tracking
        tracking['ses'] = {}
        for sub in subs:
            try:
                tracking['ses'][sub] = [x for x in listdir(f'{self.datapath}/{sub}') 
                                    if '.' not in x]
            except:
                logger.warning(f'No sessions found for {sub}')
                tracking['ses'][sub] = ['-']
            
        # Dataframe
        df = DataFrame(data=None, dtype=object)
        df.index = subs
        df['ses'] = '-'
        for x in df.index:
            df.loc[x,'ses'] = tracking['ses'][x]
        
        # Loop through other steps
        if step: 
            df, tracking = track_processing(self, step, subs, tracking, df, chan, 
                                                  stage, show, log)

        # Update tracking
        try:
            self.tracking = self.tracking | tracking
        except:
            self.tracking = {**self.tracking, **tracking}
        
        if show:
            logger.info('')
            logger.info(df)
        if outfile:
            df.to_csv(f'{self.outpath}/audit/{outfile}')

        return   

    def make_bids(self, subs = 'all', origin = 'SCN', filetype = '.edf',
                  outfile = True):
        # Set up logging
        if outfile == True:
            today = date.today().strftime("%Y%m%d")
            now = datetime.now().strftime("%H%M%S")
            logfile = f'{self.log_dir}/make_bids_{today}_{now}_log.txt'
            logger = create_logger_outfile(logfile=logfile, name='Make bids')
            logger.info('')
            logger.info(f"-------------- New call of 'Load sleep stages' evoked at {now} --------------")
        elif outfile:
            logfile = f'{self.log_dir}/{outfile}'
            logger = create_logger_outfile(logfile=logfile, name='Make bids')
        else:
            logger = create_logger('Make bids')
        logger.info('')
        logger.debug('Formatting dataset into BIDS.')
        
        make_bids(self.datapath, subs, origin, filetype)
        
    def extract_channels(self, exclude = None):
        extract_channels(self.datapath, exclude)
        
    def QC_channels(self, subs = 'all', sessions = 'all', filetype = '.edf',
                    filt = None, chantype = ['eeg', 'eog', 'emg', 'ecg'], 
                    outfile=True):
        
        """
        Run automated quality control (QC) on physiological channels for all subjects and sessions.
    
        This method initializes a SQUID instance and performs QC on the specified channel types
        (e.g., EEG, EOG, EMG, ECG). For each channel, signal quality metrics are computed, such as:
    
        - `time_not_flatline`:        Percentage of time the signal is not flatlined.
        - `time_above_10`:            Percentage of signal above 10 µV.
        - `time_below_200`:           Percentage of signal below 200 µV.
        - `gini_coeff`:               Gini coefficient of signal inequality.
        - `inverse_power_ratio`:      Inverse power ratio of high-to-low frequency content.
        - `ecg_artefact`:             Mean correlation between EEG and ECG signals.
        - `ecg_artefact_perc`:        Percentage of signal with ECG artefact correlation r > 0.5.
    
        Parameters
        ----------
        subs : str or list, optional
            Subjects to include in QC. Can be 'all' or a list of BIDS-compatible subject IDs.
    
        sessions : str or list, optional
            Sessions to include in QC. Can be 'all' or a list of session labels.
    
        filetype : str, default='.edf'
            File format to look for in BIDS dataset (e.g., '.edf', '.set', '.vhdr').
            
        filt : list or None
            Explicit list of channels to process for EEG. If None, then standard 10-20 EEG channel 
            names will be used. 
    
        chantype : list of str, default=['eeg', 'eog', 'emg', 'ecg']
            Channel types to process. Accepted values: ['eeg', 'eog', 'emg', 'ecg'].
    
        outfile : bool or str, default=True
            If True, writes log to auto-generated timestamped file in `self.log_dir`.
            If str, writes log to the specified file path.
            If False, logs only to console.
    
        Returns
        -------
        None
            Results are stored in the SQUID object and/or exported downstream.
        """
        
        # Set up logging
        if outfile == True:
            subs_str, ses_str = out_names(subs, sessions)
            today = date.today().strftime("%Y%m%d")
            now = datetime.now().strftime("%H%M%S")
            logfile = f'{self.log_dir}/load_sleep_stages_subs-{subs_str}_ses-{ses_str}_{today}_{now}_log.txt'
            logger = create_logger_outfile(logfile = logfile, name='Channel QC')
            logger.info('')
            logger.info(f"-------------- New call of 'Channel QC' evoked at {now} --------------")
        elif outfile:
            logfile = f'{self.log_dir}/{outfile}'
            logger = create_logger_outfile(logfile=logfile, name='Channel QC')
        else:
            logger = create_logger('Channel QC')
        logger.info('')
        
        
        # Check chantypes
        valid_chantypes = {'eeg', 'eog', 'emg', 'ecg'}
        if not set(chantype).issubset(valid_chantypes):
            raise ValueError(f"Invalid chantype(s) specified. Allowed: {valid_chantypes}")
        
        # Run detection
        if subs == 'all':
            subs = None
        if sessions == 'all':
            sessions = None
        squid = SQUID(self.rootpath, filetype = filetype, subjects = subs, 
                     sessions = sessions )
        squid.process_all(filt, chantype)
        
    def QC_summary(self, qc_dir = None, chantype = ['eeg', 'eog', 'emg', 'ecg']):
        
        if not qc_dir:
            qc_root = self.outpath + '/QC'
            
        # output directory    
        out_dir = self.outpath + '/datasets/'
        if not path.exists(out_dir ):
            mkdir(out_dir)
        out_dir = self.outpath + '/datasets/QC'
        if not path.exists(out_dir):
            mkdir(out_dir)
        
        # Run and save
        for modality in chantype:
            summary = gather_qc_reports(qc_root, modality)
            summary.to_csv(f'{out_dir}/summary_{modality}.csv')
            
    
    def load_stages(self, xml_dir = None, subs = 'all', sessions = 'all', 
                          filetype = '.edf', stage_key = None, 
                          outfile = True):
        '''
            Extracts stages from the BIDS formatted dataset, in which
            staging has been listed in a file *acq-PSGScoring_events.tsv, and
            saves the information in an annotations (.xml) file
        '''
        
        # Set up logging
        if outfile == True:
            subs_str, ses_str = out_names(subs, sessions)
            today = date.today().strftime("%Y%m%d")
            now = datetime.now().strftime("%H%M%S")
            logfile = f'{self.log_dir}/load_sleep_stages_subs-{subs_str}_ses-{ses_str}_{today}_{now}_log.txt'
            logger = create_logger_outfile(logfile=logfile, name='Load sleep stages')
            logger.info('')
            logger.info(f"-------------- New call of 'Load sleep stages' evoked at {now} --------------")
        elif outfile:
            logfile = f'{self.log_dir}/{outfile}'
            logger = create_logger_outfile(logfile=logfile, name='Load sleep stages')
        else:
            logger = create_logger('Load sleep stages')
        logger.info('')
        
        # Set xml_dir
        if not xml_dir:
            xml_dir = self.outpath + '/staging_manual'
        if not path.exists(xml_dir):
            mkdir(xml_dir)
        
        # Load stages
        flag = load_stages(self.datapath, xml_dir, subs, sessions, filetype, 
                           stage_key)
        
        # Log finish 
        if flag > 0:
            logger.warning(f"'load_stages' finished with {flag} WARNINGS. See log for detail.")
        else:
            logger.debug("'load_stages' finished without error.")
    #--------------------------------------------------------------------------
    '''
    ANALYSIS FUNCTIONS
    
    power_spectrum -> performs power spectral analysis.
                    
    
    '''    
    
    
    def power_spectrum(self, xml_dir = None, out_dir = None, 
                             subs = 'all', sessions = 'all', filetype = '.edf',  
                             chan = None, ref_chan = None, 
                             grp_name = 'eeg', rater = None, 
                             stage = ['NREM1','NREM2','NREM3', 'REM'], 
                             cycle_idx = None, concat_cycle = True, 
                             concat_stage = False, general_opts = None, 
                             frequency_opts = None, filter_opts = None, 
                             epoch_opts = None, event_opts = None, 
                             norm = None, norm_opts = None, 
                             outfile = True):
        
        # Set up logging
        if outfile == True:
            subs_str, ses_str = out_names(subs, sessions)
            today = date.today().strftime("%Y%m%d")
            now = datetime.now().strftime("%H%M%S")
            logfile = f'{self.log_dir}/detect_power_spectrum_subs-{subs_str}_ses-{ses_str}_{today}_{now}_log.txt'
            logger = create_logger_outfile(logfile=logfile, name='Power spectrum')
            logger.info('')
            logger.info(f"-------------- New call of 'Power spectrum' evoked at {now} --------------")
        elif outfile:
            logfile = f'{self.log_dir}/{outfile}'
            logger = create_logger_outfile(logfile=logfile, name='Power spectrum')
        else:
            logger = create_logger('Power spectrum')
        logger.info('')

        
        # Set input/output directories
        in_dir = self.datapath
        
        if not xml_dir:
            xml_dir = select_input_dirs(self.outpath, xml_dir, 'staging') 
            
        if not path.exists(xml_dir):
            logger.info('')
            logger.critical(f"{xml_dir} doesn't exist. Sleep staging has not been "
                            "run or hasn't been converted correctly.")
            logger.info('Check documentation for how to set up staging data:')
            logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
            logger.info('-' * 10)
            return
        else:    
            logger.debug(f'Input annotations being read from: {xml_dir}')

        if not out_dir:
            out_dir = f'{self.outpath}/powerspectrum' 
        if not path.exists(out_dir):
            mkdir(out_dir)
        logger.debug(f'Output being saved to: {out_dir}')
        
        # Set channels
        chan, ref_chan = check_chans(self.rootpath, chan, ref_chan, logger)
        if isinstance(chan, str):
            return
        
        # Set default parameters
        if not general_opts:
            general_opts = default_general_opts()
        if not frequency_opts:
            frequency_opts = default_frequency_opts()
        if not epoch_opts:
            epoch_opts = default_epoch_opts()  
        if not event_opts:
            event_opts = default_event_opts()
        if not norm_opts:
            norm_opts = default_norm_opts()
        if not filter_opts:
            filter_opts = default_filter_opts()    
        frequency_opts['frequency'] = (filter_opts['highpass'], filter_opts['lowpass'])
        
        # Format concatenation
        cat = (int(concat_cycle),int(concat_stage),
               1,
               int(event_opts['concat_events']),
               )
        
        # Set suffix for output filename
        if not general_opts['suffix']:
            general_opts['suffix'] = f"{frequency_opts['frequency'][0]}-{frequency_opts['frequency'][1]}Hz"
        
        # Check annotations directory exists, run detection
        spectrum = Spectrum(in_dir, xml_dir, out_dir, chan, ref_chan, 
                            grp_name, stage, cat, rater, cycle_idx, subs, 
                            sessions, self.tracking)
                
        spectrum.powerspec_it(general_opts, frequency_opts, filter_opts, 
                              epoch_opts, event_opts, norm, norm_opts, 
                              filetype, logger) 
        
        try:
            self.tracking = self.tracking | spectrum.tracking
        except:
            self.tracking = {**self.tracking, **spectrum.tracking}
            
        return 
        

    #--------------------------------------------------------------------------
    '''
    SLEEP EVENTS DETECTIONS
    
    sleep_staging
    detect_artefacts,
    detect_spectral_peaks,
    detect_slow_oscillations,
    detect_spindles,
    
    
    '''
    def detect_sleep_stages(self, out_dir = None, 
                                  subs = 'all', sessions = 'all', filetype = '.edf', 
                                  method = 'Vallat2021', qual_thresh = 0.5,
                                  eeg_chan = None, ref_chan = None, 
                                  eog_chan = None, emg_chan = None, 
                                  rater = None, invert = False, outfile = True):
        
        from seapipe.events.seabass import seabass

        # Set up logging
        if outfile == True:
            subs_str, ses_str = out_names(subs, sessions)
            today = date.today().strftime("%Y%m%d")
            now = datetime.now().strftime("%H:%M:%S")
            logfile = f'{self.log_dir}/detect_sleep_stages_subs-{subs_str}_ses-{ses_str}_{today}_{now}_log.txt'
            logger = create_logger_outfile(logfile=logfile, name='Detect sleep stages')
            logger.info('')
            logger.info(f"-------------- New call of 'Detect sleep stages' evoked at {today}, {now} --------------")
        elif outfile:
            logfile = f'{self.log_dir}/{outfile}'
            logger = create_logger_outfile(logfile=logfile, name='Detect sleep stages')
        else:
            logger = create_logger('Detect sleep stages')
        logger.info('')
        
        # Set input/output directories
        in_dir = self.datapath
        if not out_dir:
            out_dir = f'{self.outpath}/staging_auto'    
        if not path.exists(out_dir):
            mkdir(out_dir)
        logger.debug(f'Output annotations being saved to: {out_dir}')
        
        # Check subs
        track_sheet = read_tracking_sheet(self.rootpath, logger)
        if not subs:
            subs = [x for x in list(set(track_sheet['sub']))]
            subs.sort()
        if not sessions:
            sessions = track_sheet
        
        # Set channels
        eeg_chan, ref_chan = check_chans(self.rootpath, eeg_chan, ref_chan, logger)
        if isinstance(eeg_chan, str) and eeg_chan == 'error':
            return
        
        # Check inversion
        if invert == None:
            invert = check_chans(self.rootpath, None, False, logger)
        elif type(invert) != bool:
            logger.critical("The argument 'invert' must be set to either: "
                            f"'True', 'False' or 'None'; but it was set as {invert}.")
            logger.info('')
            logger.info("Check documentation for how to set up staging data: "
                        "https://seapipe.readthedocs.io/en/latest/index.html")
            logger.info('-' * 10)
            logger.critical('Sleep stage detection finished with ERRORS. See log for details.')
            return
    
        # Check annotations directory exists, run detection
        stages = seabass(in_dir, out_dir, eeg_chan, ref_chan, eog_chan, emg_chan, 
                         rater, subs, sessions, self.tracking) 
        stages.detect_stages(method, qual_thresh, invert, filetype, track_sheet,
                             logger)
        try:
            self.tracking = self.tracking | stages.tracking
        except:
            self.tracking = {**self.tracking, **stages.tracking}
        return
    
    
    def detect_artefacts(self, xml_dir = None, out_dir = None, 
                               subs = 'all', sessions = 'all', filetype = '.edf', 
                               method = 'seapipe', win_size = 5,
                               chan = None, ref_chan = None, 
                               label = 'individual',
                               rater = None, grp_name = 'eeg', 
                               stage = ['NREM1', 'NREM2', 'NREM3', 'REM'],
                               outfile = True):
        
        from seapipe.events.sand import SAND

        # Set up logging
        if outfile == True:
            subs_str, ses_str = out_names(subs, sessions)
            today = date.today().strftime("%Y%m%d")
            now = datetime.now().strftime("%H%M%S")
            logfile = f'{self.log_dir}/detect_artefacts_subs-{subs_str}_ses-{ses_str}_{today}_{now}_log.txt'
            logger = create_logger_outfile(logfile=logfile, name='Detect artefacts')
            logger.info('')
            logger.info(f"-------------- New call of 'Detect artefacts' evoked at {now} --------------")
        elif outfile:
            logfile = f'{self.log_dir}/{outfile}'
            logger = create_logger_outfile(logfile=logfile, name='Detect artefacts')
        else:
            logger = create_logger('Detect artefacts')
        logger.info('')
        
        # Set input/output directories
        in_dir = self.datapath
        if not xml_dir:
            xml_dir = select_input_dirs(self.outpath, xml_dir, 'staging') 
            
        if not path.exists(xml_dir):
            logger.info('')
            logger.critical(f"{xml_dir} doesn't exist. Sleep staging has not been "
                            "run or hasn't been converted correctly.")
            logger.info('Check documentation for how to set up staging data:')
            logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
            logger.info('-' * 10)
            return
        else:    
            logger.debug(f'Input annotations being read from: {xml_dir}')
            
        if not out_dir:
            out_dir = select_input_dirs(self.outpath, xml_dir, 'staging')   
        if not path.exists(out_dir):
            mkdir(out_dir)
        logger.debug(f'Output annotations being saved to: {out_dir}')
        
        # Check subs
        if not subs:
            tracking = read_tracking_sheet(self.rootpath, logger)
            subs = [x for x in list(set(tracking['sub']))]
            subs.sort()
        if not sessions:
            sessions = read_tracking_sheet(self.rootpath, logger)
        
        # Set channels
        chan, ref_chan = check_chans(self.rootpath, chan, ref_chan, 
                                          logger)
        
    
        # Check annotations directory exists, run detection
        artefacts = SAND(in_dir, xml_dir, out_dir, chan, ref_chan, 
                         rater, grp_name, subs, sessions, self.tracking) 
        artefacts.detect_artefacts(method, label,  win_size,  filetype, stage,
                                   logger)
                                   
    
        try:
            self.tracking = self.tracking | artefacts.tracking
        except:
            self.tracking = {**self.tracking, **artefacts.tracking}
        return
        
        
    
    def detect_spectral_peaks(self, xml_dir = None, out_dir = None, 
                                    subs = 'all', sessions = 'all', chan = None, 
                                    ref_chan = None, grp_name = 'eeg', 
                                    rater = None, frequency = (9,16), 
                                    stage = ['NREM2','NREM3'], cycle_idx = None,
                                    concat_cycle = True, concat_stage = False,
                                    general_opts = None, frequency_opts = None,
                                    filter_opts = None, epoch_opts = None, 
                                    event_opts = None, fooof_opts = None, 
                                    filetype = '.edf', suffix = None, 
                                    outfile = True):
        
        # Set up logging
        if outfile == True:
            subs_str, ses_str = out_names(subs, sessions)
            today = date.today().strftime("%Y%m%d")
            now = datetime.now().strftime("%H%M%S")
            logfile = f'{self.log_dir}/detect_specparams_subs-{subs_str}_ses-{ses_str}_{today}_{now}_log.txt'
            logger = create_logger_outfile(logfile = logfile, 
                                           name = 'Detect spectral peaks')
            logger.info('')
            logger.info(f"-------------- New call of 'Detect spectral peaks' evoked at {now} --------------")
        elif outfile:
            logfile = f'{self.log_dir}/{outfile}'
            logger = create_logger_outfile(logfile = logfile, 
                                           name = 'Detect spectral peaks')
        else:
            logger = create_logger('Detect spectral peaks')
        
        # Set input/output directories
        in_dir = self.datapath
        if not xml_dir:
            xml_dir = select_input_dirs(self.outpath, xml_dir, 'staging') 
            
        if not path.exists(xml_dir):
            logger.info('')
            logger.critical(f"{xml_dir} doesn't exist. Sleep staging has not been "
                            "run or hasn't been converted correctly.")
            logger.info('Check documentation for how to set up staging data:')
            logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
            logger.info('-' * 10)
            return
        else:    
            logger.debug(f'Input annotations being read from: {xml_dir}')
            
        if not out_dir:
            out_dir = f'{self.outpath}/fooof' 
        if not path.exists(out_dir):
            mkdir(out_dir)
        logger.debug(f'Output being saved to: {out_dir}')
            
        # Check subs
        if not subs:
            tracking = read_tracking_sheet(self.rootpath, logger)
            subs = [x for x in list(set(tracking['sub']))]
            subs.sort()    
        if not sessions:
            sessions = read_tracking_sheet(self.rootpath, logger)
            
        # Set channels
        chan, ref_chan = check_chans(self.rootpath, chan, ref_chan, logger)
        
        # Format concatenation
        cat = (int(concat_cycle),int(concat_stage),1,1)
        
        # Run detection
        spectrum = Spectrum(in_dir, xml_dir, out_dir, chan, 
                            ref_chan, grp_name, stage, cat, rater, 
                            cycle_idx, subs, sessions, self.tracking)
        
        if not general_opts:
            general_opts = default_general_opts()
        if not frequency_opts:
            frequency_opts = default_frequency_opts()
        if not filter_opts:
            filter_opts = default_filter_opts()
        if not epoch_opts:
            epoch_opts = default_epoch_opts()  
        if not event_opts:
            event_opts = default_event_opts()
        if not fooof_opts:
            fooof_opts = default_fooof_opts() 
            
        fooof_opts['bands_fooof'] = [frequency]
        
        # Set suffix for output filename
        if not suffix:
            general_opts['suffix'] = f'{frequency[0]}-{frequency[1]}Hz'
        
        spectrum.fooof_it(general_opts, frequency_opts, filter_opts, 
                          epoch_opts, event_opts, fooof_opts, 
                          filetype, logger)  
                
        return 
    
    
    def detect_slow_oscillations(self, xml_dir=None, out_dir=None, subs='all', 
                                       sessions='all', filetype='.edf', 
                                       method = ['Staresina2015'], chan=None,
                                       ref_chan=None, rater=None, grp_name='eeg', 
                                       stage = ['NREM2','NREM3'], cycle_idx=None, 
                                       duration=(0.2, 2), invert = None,
                                       reject_artf = ['Artefact', 'Arou', 'Arousal'],
                                       average_channels = False, outfile = True):
        
        from seapipe.events.seasnakes import seasnakes

        # Set up logging
        if outfile == True:
            subs_str, ses_str = out_names(subs, sessions)
            evt_out = '_'.join(method)
            today = date.today().strftime("%Y%m%d")
            now = datetime.now().strftime("%H%M%S")
            logfile = f'{self.log_dir}/detect_slowosc_{evt_out}_subs-{subs_str}_ses-{ses_str}_{today}_{now}_log.txt'
            logger = create_logger_outfile(logfile=logfile, name='Detect slow oscillations')
            logger.info('')
            logger.info(f"-------------- New call of 'Detect slow oscillations' evoked at {now} --------------")
        elif outfile:
            logfile = f'{self.log_dir}/{outfile}'
            logger = create_logger_outfile(logfile=logfile, name='Detect slow oscillations')
        else:
            logger = create_logger('Detect slow oscillations')
        
        logger.info('')
        logger.debug("Commencing SO detection pipeline.")
        logger.info('')
        
        # Set input/output directories
        in_dir = self.datapath
        if not xml_dir:
            xml_dir = select_input_dirs(self.outpath, xml_dir, 'staging') 
            
        if not path.exists(xml_dir):
            logger.info('')
            logger.critical(f"{xml_dir} doesn't exist. Sleep staging has not been "
                            "run or hasn't been converted correctly.")
            logger.info('Check documentation for how to set up staging data:')
            logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
            logger.info('-' * 10)
            return
        else:    
            logger.debug(f'Input annotations being read from: {xml_dir}')
        
        if not out_dir:
            out_dir = f'{self.outpath}/slowwave'    
        if not path.exists(out_dir):
            mkdir(out_dir)
        logger.debug(f'Output annotations being saved to: {out_dir}')
        
        # Check subs
        if not subs:
            tracking = read_tracking_sheet(self.rootpath, logger)
            subs = [x for x in list(set(tracking['sub']))]
            subs.sort()
        if not sessions:
            sessions = read_tracking_sheet(self.rootpath, logger)
            
        # Set channels
        chan, ref_chan = check_chans(self.rootpath, chan, ref_chan, logger)
        
        # Check inversion
        if invert == None:
            invert = check_chans(self.rootpath, None, False, logger)
        elif type(invert) != bool:
            logger.critical(f"The argument 'invert' must be set to either: 'True', "
                            f"'False' or 'None'; but it was set as {invert}.")
            logger.info('Check documentation for how to set up staging data:')
            logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
            logger.info('-' * 10)
            logger.critical('SO detection finished with ERRORS. See log for details.')
            return
            
        # Format concatenation
        cat = (1,1,1,1)
    
        # Run detection
        SO = seasnakes(in_dir, xml_dir, out_dir, chan, ref_chan, 
                         grp_name, stage, rater, subs, sessions, 
                         self.tracking, reject_artf) 
        SO.detect_slowosc(method, cat, cycle_idx, duration, 
                               average_channels, invert, filetype, logger)
        try:
            self.tracking = self.tracking | SO.tracking
        except:
            self.tracking = {**self.tracking, **SO.tracking}
        
        return
    
    
    def detect_spindles(self, xml_dir = None, out_dir = None, subs = 'all', 
                              sessions = 'all', filetype = '.edf', 
                              method = ['Moelle2011'], chan = None, 
                              ref_chan = None, rater = None, 
                              stage = ['NREM2','NREM3'], grp_name = 'eeg', 
                              cycle_idx = None, concat_cycle = True, 
                              frequency = None, adap_bands = 'Fixed', 
                              adap_bw = 4, duration = (0.5, 3),
                              reject_artf = ['Artefact', 'Arou', 'Arousal'], 
                              outfile = True):
        
        # Set up logging
        if outfile == True:
            subs_str, ses_str = out_names(subs, sessions)
            evt_out = '_'.join(method)
            today = date.today().strftime("%Y%m%d")
            now = datetime.now().strftime("%H%M%S")
            logfile = f'{self.log_dir}/detect_spindles_{evt_out}_subs-{subs_str}_ses-{ses_str}_{today}_{now}_log.txt'
            logger = create_logger_outfile(logfile=logfile, name='Detect spindles')
            logger.info('')
            logger.info(f"-------------- New call of 'Detect spindles' evoked at {now} --------------")
        elif outfile:
            logfile = f'{self.log_dir}/{outfile}'
            logger = create_logger_outfile(logfile=logfile, name='Detect spindles')
        else:
            logger = create_logger('Detect spindles')
            

        logger.info('')
        logger.debug("Commencing spindle detection pipeline.")
        logger.info('')
        
        # Set input/output directories
        in_dir = self.datapath   
        if not xml_dir:
            xml_dir = select_input_dirs(self.outpath, xml_dir, 'staging') 
            
        if not path.exists(xml_dir):
            logger.info('')
            logger.critical(f"{xml_dir} doesn't exist. Sleep staging has not been "
                            "run or hasn't been converted correctly.")
            logger.info('Check documentation for how to set up staging data:')
            logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
            logger.info('-' * 10)
            return
        else:    
            logger.debug(f'Input annotations being read from: {xml_dir}')
        
        if not out_dir:
            for met in method:
                out_dir = select_output_dirs(self.outpath, out_dir, met)  
        if not path.exists(out_dir):
            mkdir(out_dir)
        logger.debug(f'Output annotations being saved to: {out_dir}')
        
        # Check subs
        logger.warning(
                            f"[DEBUG whales] entering fallback check: "
                            f"subs={subs!r} (type={type(subs)}), "
                            f"sessions={sessions!r} (type={type(sessions)})"
                        )
        if not subs:
            tracking = read_tracking_sheet(self.rootpath, logger)
            subs = [x for x in list(set(tracking['sub']))]
            subs.sort()
        if not sessions:
            sessions = read_tracking_sheet(self.rootpath, logger)
            
        # Set channels
        chan, ref_chan = check_chans(self.rootpath, chan, ref_chan, logger)
        if not isinstance(chan, DataFrame) and not isinstance(chan, list):
            logger.error('Problem loading channel information')
            return
        elif isinstance(ref_chan, str):
            logger.error('Problem loading ref-channel information')
            return
        
        # Format concatenation
        if concat_cycle == True:
            cat = (1,0,1,1)
        else:
            cat = (0,0,1,1)
            if not cycle_idx:
                logger.error("'concat_cycle' is set to false, but 'cycle_idx' = None. "
                             "Set cycle_idx to a list of integers to use cycles properly.")
                logger.info("Check documentation for how to mark and use sleep cycles:")
                logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
                return
        
        # Check for adapted bands
        logger.debug(f'Detection using {adap_bands} frequency bands has been selected.')
        frequency, flag = adap_bands_setup(self, adap_bands, frequency, subs, 
                                           sessions, chan, ref_chan, stage, cat, 
                                           concat_cycle, cycle_idx, logger)
        if flag == 'error':
            return
       
        # Run detection
        self.track(subs, sessions, step = ['fooof','spindle'], show = False, 
                   log = False)
        
        #self.track(step='fooof', show = False, log = False)
        #self.track(step='spindle', show = False, log = False)
        spindle = whales(self.rootpath, in_dir, xml_dir, out_dir, 
                         chan, ref_chan, grp_name, stage, frequency, rater, 
                         subs, sessions, reject_artf, self.tracking) 
        spindle.whale_it(method, cat, cycle_idx, adap_bands, adap_bw, 
                         duration, filetype, logger)
        try:
            self.tracking = self.tracking | spindle.tracking
        except:
            self.tracking = {**self.tracking, **spindle.tracking}
            
        return
    
    
    def whales(self, xml_dir = None, out_dir = None, subs = 'all', 
                     sessions = 'all', filetype = '.edf', 
                     method = ['Moelle2011', 'Ray2015'], evt_out = 'spindle',
                     merge_type = 'consensus', weights = None,
                     chan = None, ref_chan = None, rater = None, 
                     stage = ['NREM2','NREM3'], grp_name = 'eeg', 
                     cycle_idx = None, s_freq = None, keyword = None, 
                     duration =(0.5, 3),
                     reject_artf = ['Artefact', 'Arou', 'Arousal'], 
                     outfile = True):
        
        # Set up logging
        if outfile == True:
            subs_str, ses_str = out_names(subs, sessions)
            today = date.today().strftime("%Y%m%d")
            now = datetime.now().strftime("%H%M%S")
            logfile = f'{self.log_dir}/detect_spindles_WHALES_subs-{subs_str}_ses-{ses_str}_{today}_{now}_log.txt'
            logger = create_logger_outfile(logfile=logfile, name='Detect spindles (WHALES)')
            logger.info('')
            logger.info(f"-------------- New call of 'Detect spindles (WHALES)' evoked at {now} --------------")
       
        elif outfile:
            logfile = f'{self.log_dir}/{outfile}'
            logger = create_logger_outfile(logfile=logfile, name='Detect spindles (WHALES)')
        else:
            logger = create_logger('Detect spindles (WHALES)')
            
        logger.info('')
        logger.debug("Commencing spindle optimisation pipeline.")
        logger.info('')
        
        # Set input/output directories
        in_dir = self.datapath
        xml_dir = select_input_dirs(self.outpath, xml_dir, 'spindle') 
        logger.debug(f'Input annotations being read from: {xml_dir}')
        
        if not out_dir:
            out_dir = xml_dir  
        if not path.exists(out_dir):
            mkdir(out_dir)
        logger.debug(f'Output annotations being saved to: {out_dir}')
        
        if merge_type == 'consensus':
            cs_thresh = 0.5
        elif merge_type == 'addition':
            cs_thresh = 0.01
        
        # Check subs
        if not subs:
            tracking = read_tracking_sheet(self.rootpath, logger)
            subs = [x for x in list(set(tracking['sub']))]
            subs.sort()
        if not sessions:
            sessions = read_tracking_sheet(self.rootpath, logger)
        
        # Set channels
        chan, ref_chan = check_chans(self.rootpath, chan, ref_chan, logger)
        if not isinstance(chan, DataFrame) and not isinstance(chan, list):
            logger.error('Problem loading channel information')
            return
        elif isinstance(ref_chan, str):
            logger.error('Problem loading ref-channel information')
            return
        
        self.track(step='spindle', show = False, log = False)
        
        logger.debug('Starting Merge Now')
        spindle = whales(self.rootpath, in_dir, xml_dir, out_dir, chan, ref_chan, 
                         grp_name, stage, frequency = None, rater = rater, 
                         subs = subs, sessions = sessions, 
                         reject_artf = reject_artf, tracking = self.tracking) 
        spindle.whales(method, merge_type, chan, rater, stage, ref_chan, grp_name, 
                       keyword, cs_thresh, s_freq, duration, evt_out, weights,
                       filetype, logger)
    
    
    def detect_rems(self, xml_dir = None, out_dir = None, subs = 'all', 
                          sessions = 'all', filetype = '.edf', 
                          method = ['YASA'], chan = None,
                          ref_chan = None, rater = None, grp_name = 'eeg', 
                          stage = ['REM'], cycle_idx = None, 
                          amplitude = (50, 325), duration = (0.3, 1.5),
                          reject_artf = ['Artefact', 'Arou', 'Arousal'],
                          average_channels = False, outfile = True):
        
        from seapipe.events.remora import remora

        # Set up logging
        if outfile == True:
            subs_str, ses_str = out_names(subs, sessions)
            evt_out = '_'.join(method)
            today = date.today().strftime("%Y%m%d")
            now = datetime.now().strftime("%H%M%S")
            logfile = f'{self.log_dir}/detect_rems_{evt_out}_subs-{subs_str}_ses-{ses_str}_{today}_{now}_log.txt'
            logger = create_logger_outfile(logfile=logfile, name='Detect eye movements (REMS)')
            logger.info('')
            logger.info(f"-------------- New call of 'Detect rapid eye movements' evoked at {now} --------------")
            
        elif outfile:
            logfile = f'{self.log_dir}/{outfile}'
            logger = create_logger_outfile(logfile=logfile, name='Detect Detect REMS')
        else:
            logger = create_logger('Detect Detect REMS')
        
        logger.info('')
        logger.debug("Commencing REMS detection pipeline.")
        logger.info('')
        
        # Set input/output directories
        in_dir = self.datapath
        if not xml_dir:
            xml_dir = select_input_dirs(self.outpath, xml_dir, 'staging') 
            
        if not path.exists(xml_dir):
            logger.info('')
            logger.critical(f"{xml_dir} doesn't exist. Sleep staging has not been "
                            "run or hasn't been converted correctly.")
            logger.info('Check documentation for how to set up staging data:')
            logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
            logger.info('-' * 10)
            return
        else:    
            logger.debug(f'Input annotations being read from: {xml_dir}')
        
        if not out_dir:
            out_dir = f'{self.outpath}/rems'    
        if not path.exists(out_dir):
            mkdir(out_dir)
        logger.debug(f'Output annotations being saved to: {out_dir}')
        
        # Check subs
        if not subs:
            tracking = read_tracking_sheet(self.rootpath, logger)
            subs = [x for x in list(set(tracking['sub']))]
            subs.sort()
        if not sessions:
            sessions = read_tracking_sheet(self.rootpath, logger)
            
        # Set channels
        chan, ref_chan = check_chans(self.rootpath, chan, ref_chan, logger)
    
        # Run detection
        REMS = remora(in_dir, xml_dir, out_dir, chan, ref_chan, 
                      rater, grp_name, reject_artf, subs, sessions)
        
        REMS.detect_rems(method, amplitude, duration, stage, cycle_idx, 
                         filetype, logger)
        try:
            self.tracking = self.tracking | REMS.tracking
        except:
            self.tracking = {**self.tracking, **REMS.tracking}
            
        return
    
    #--------------------------------------------------------------------------
    '''
    PLOTTING.
    
    event_spectrogram ->
    
    
    '''
    
    
    def spectrogram(self, xml_dir = None, out_dir = None, subs = 'all', 
                          sessions = 'all', filetype = '.edf', chan = None, 
                          ref_chan = None, rater = None, stage = None, 
                          grp_name = 'eeg', cycle_idx = None, 
                          concat_stage = False, concat_cycle = True, 
                          evt_type = None, buffer = 0, invert = None, 
                          filter_opts = None, progress=True, outfile=False):
        
        from seapipe.spectrum.spectrogram import event_spectrogram

        # Set up logging
        logger = create_logger('Event spectrogram')
        logger.info('')
        logger.debug("Creating spectrogram of events.")
        logger.info('')
        
        # Set input/output directories
        in_dir = self.datapath
            
        if not xml_dir:
            xml_dir = select_input_dirs(self.outpath, xml_dir, 'staging') 
            
        if not path.exists(xml_dir):
            logger.info('')
            logger.critical(f"{xml_dir} doesn't exist. Sleep staging has not been "
                            "run or hasn't been converted correctly.")
            logger.info('Check documentation for how to set up staging data:')
            logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
            logger.info('-' * 10)
            return
        else:    
            logger.debug(f'Input annotations being read from: {xml_dir}')
        
        if not out_dir:
            out_dir = f'{self.outpath}/spindle'    
        if not path.exists(out_dir):
            mkdir(out_dir)
        logger.debug(f'Output annotations being saved to: {out_dir}')
        
        # Check subs
        if not subs:
            tracking = read_tracking_sheet(self.rootpath, logger)
            subs = [x for x in list(set(tracking['sub']))]
            subs.sort()
        if not sessions:
            sessions = read_tracking_sheet(self.rootpath, logger)
            
        # Format concatenation
        cat = (int(concat_cycle),int(concat_stage),1,1)
        
        # Check inversion
        if invert == None:
            invert = check_chans(self.rootpath, chan, False, logger)
        elif type(invert) != bool:
            logger.critical("The argument 'invert' must be set to either: "
                            f"'True', 'False' or 'None'; but it was set as {invert}.")
            logger.info('Check documentation for how to set up staging data:')
            logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
            logger.info('-' * 10)
            return
        
        if not filter_opts:
            filter_opts = default_filter_opts()
        
        
        if not evt_type:
            logger.warning("No event type (evt_type) has been specified. "
                           "Spectrogram will be run on epochs instead. "
                           "This may take some time...")
        
        event_spectrogram(self, in_dir, xml_dir, out_dir, subs, sessions, stage, 
                                cycle_idx, chan, ref_chan, rater, grp_name, 
                                evt_type, buffer, invert, cat, filter_opts,  
                                outfile, filetype, progress, self.tracking)
        
        return
    
    
    #--------------------------------------------------------------------------
    '''
    PHASE AMPLITUDE COUPLING.
    
    mean_amps -> runs Phase Amplitude Coupling analyses on sleep EEG data. 
    
    
    '''    
    def pac(self, xml_dir = None, out_dir = None, subs = 'all', sessions = 'all', 
                  filetype = '.edf', chan = None, ref_chan = None, rater = None, 
                  grp_name = 'eeg', stage = ['NREM2','NREM3'], concat_stage = True, 
                  cycle_idx = None, concat_cycle = True,  
                  method = 'MI', surrogate = 'Time lag', correction = 'Z-score',
                  adap_bands_phase = 'Fixed', frequency_phase = (0.5, 1.25), 
                  adap_bands_amplitude = 'Fixed', frequency_amplitude = (11, 16),
                  adap_bw = 4, min_dur = 1, nbins = 18, invert = None,
                  frequency_opts = None, filter_opts = None, epoch_opts = None, 
                  evt_name = None, event_opts = None, 
                  reject_artf = ['Artefact', 'Arou', 'Arousal'], 
                  progress = True, outfile = True):
        
        from seapipe.pac.octopus import octopus, pac_method
        from seapipe.pac.pacats import pacats

        # Set up logging
        if outfile == True:
            subs_str, ses_str = out_names(subs, sessions)
            today = date.today().strftime("%Y%m%d")
            now = datetime.now().strftime("%H%M%S")
            pha = f'{frequency_phase[0]}-{frequency_phase[1]}'
            amp = f'{frequency_amplitude[0]}-{frequency_amplitude[1]}'
            logfile = f'{self.log_dir}/pac_{pha}_{amp}_subs-{subs_str}_ses-{ses_str}_{today}_{now}_log.txt'
            logger = create_logger_outfile(logfile=logfile, name='Phase-amplitude coupling')
            logger.info('')
            logger.info(f"-------------- New call of 'Phase-amplitude coupling' evoked at {now} --------------")
        elif outfile:
            logfile = f'{self.log_dir}/{outfile}'
            logger = create_logger_outfile(logfile=logfile, name='Phase-amplitude coupling')
        else:
            logger = create_logger('Phase-amplitude coupling')
        logger = create_logger('Phase-amplitude coupling')
        logger.info('')
        
        # Set input/output directories
        in_dir = self.datapath
        
        if not xml_dir:
            xml_dir = select_input_dirs(self.outpath, xml_dir, evt_name) 

        if not path.exists(xml_dir):
            logger.info('')
            logger.critical(f"{xml_dir} doesn't exist. Sleep staging has not been "
                            f"run or Events: {evt_name} haven't been detected.")
            logger.info('Check documentation for how to set up staging data:')
            logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
            logger.info('-' * 10)
            return
        else:    
            logger.debug(f'Input annotations being read from: {xml_dir}')
        
        # Check subs
        if not subs:
            tracking = read_tracking_sheet(self.rootpath, logger)
            subs = [x for x in list(set(tracking['sub']))]
            subs.sort()
        if not sessions:
            sessions = read_tracking_sheet(self.rootpath, logger)
            
        # Set channels
        chan, ref_chan = check_chans(self.rootpath, chan, ref_chan, logger)
        if not isinstance(chan, DataFrame) and not isinstance(chan, list):
            return
        elif isinstance(ref_chan, str):
            return
        
        # Check for adapted bands 
        cat = (int(concat_cycle),int(concat_stage),0,0)
        
        frequency_phase, flag1 = adap_bands_setup(self, adap_bands_phase, frequency_phase, 
                                     subs, sessions, chan, ref_chan, stage, cat, 
                                     concat_cycle, cycle_idx, logger)
        
        frequency_amplitude, flag2 = adap_bands_setup(self, adap_bands_amplitude, 
                                               frequency_amplitude, subs, sessions, 
                                               chan, ref_chan, stage, cat, 
                                               concat_cycle, cycle_idx, logger)
        
        if flag1 == 'error' or flag2 == 'error':
            return
        
        # Set PAC methods
        idpac = pac_method(method, surrogate, correction)
        
        # Set default parameters
        if not frequency_opts:
            frequency_opts = default_frequency_opts()
        if not epoch_opts:
            epoch_opts = default_epoch_opts()  
        if not event_opts:
            event_opts = default_event_opts()
        if not filter_opts:
            filter_opts = default_filter_opts()   
        filter_opts['bandpass'] = False
        
        # Check inversion
        if invert == None:
            invert = check_chans(self.rootpath, None, False, logger)
        elif type(invert) != bool:
            logger.critical("The argument 'invert' must be set to either: "
                            f"'True', 'False' or 'None'; but it was set as {invert}.")
            logger.info('Check documentation for how to set up staging data:')
            logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
            logger.info('-' * 10)
            logger.critical('Phase amplitude coupling finished with ERRORS. See log for details.')
            return
        
        self.track(subs, sessions, step = 'fooof', show = False, 
                   log = False)
        self.tracking['event_pac'] = {}
        
        # Check whether event based or continuous
        if evt_name: #OCTOPUS
            if not out_dir:
                out_dir = f'{self.outpath}/event_pac'  
            if not path.exists(out_dir):
                mkdir(out_dir)
            logger.debug(f'Output being saved to: {out_dir}')
                
            cat = (int(concat_cycle),int(concat_stage),0,0)
            Octopus = octopus(self.rootpath, in_dir, xml_dir, out_dir, 
                              chan, ref_chan, grp_name, stage, rater, 
                              subs, sessions, reject_artf,
                              self.tracking)
            
            Octopus.pac_it(cycle_idx, cat, nbins, filter_opts, epoch_opts, 
                           frequency_opts, event_opts, filetype, idpac, evt_name, 
                           min_dur, adap_bands_phase, frequency_phase, 
                           adap_bands_amplitude, frequency_amplitude, 
                           adap_bw, invert, progress, logger)
        else: #PACATS
            if not out_dir:
                out_dir = f'{self.outpath}/pac'  
            if not path.exists(out_dir):
                mkdir(out_dir)
            logger.debug(f'Output being saved to: {out_dir}')
            cat = (int(concat_cycle),int(concat_stage),1,1)
            Pacats = pacats(self.rootpath, in_dir, xml_dir, out_dir, 
                            chan, ref_chan, grp_name, stage, rater, subs, sessions, 
                            reject_artf, self.tracking)
            Pacats.pac_it(cycle_idx, cat, nbins, filter_opts, epoch_opts, 
                           frequency_opts, filetype, idpac, 
                           min_dur, adap_bands_phase, frequency_phase, 
                           adap_bands_amplitude, frequency_amplitude,
                           adap_bw, invert, progress, logger)

        return
    
    '''
    DYNAMICS
    
    spindle_clustering -> analysis of spindle event temporal clustering 
    
    sigma_fluctuations -> 
    '''
    
    def cluster_flucs(self, evt_name, xml_dir = None, out_dir = None,
                            freq_bands = {'SWA': (0.5, 4), 'Sigma': (10, 15)},
                            subs = 'all', sessions = 'all', 
                            filetype = '.edf',
                            chan = None, ref_chan = None, grp_name = 'eeg', 
                            stage = ['NREM2'], 
                            rater = None, 
                            outfile = True):
        
        from seapipe.events.clam import clam

        # Force evt_name into list, and loop through events    
        if isinstance(evt_name, str):
            evts = [evt_name]
        elif isinstance(evt_name, list):
            evts = evt_name
        else:
            raise TypeError(f"'evt_name' can only be a str or a list, but {type(evt_name)} was passed.")
        
        # Set up logging
        if outfile == True:
            subs_str, ses_str = out_names(subs, sessions)
            evt_out = '_'.join(evt_name)
            today = date.today().strftime("%Y%m%d")
            now = datetime.now().strftime("%H:%M:%S")
            logfile = f'{self.log_dir}/event_clustering_{evt_out}_subs-{subs_str}_ses-{ses_str}_{today}_{now}_log.txt'
            logger = create_logger_outfile(logfile=logfile, name='Event clustering')
            logger.info('')
            logger.debug(f"-------------- New call of 'Event clustering evoked at {now} --------------")
        elif outfile:
            logfile = f'{self.log_dir}/{outfile}'
            logger = create_logger_outfile(logfile=logfile, name='Event clustering')
        else:
            logger = create_logger('Event clustering')
        
        for evt_name in evts:
            
            out_dir = select_output_dirs(self.outpath, out_dir, 'cluster')
            logger.debug(f'Output being save to: {out_dir}')
            
            xml_dir = select_input_dirs(self.outpath, xml_dir, 'cluster')
            logger.debug(f'Input annotations being read from: {xml_dir}')
            
            # Check annotations directory exists
            if not path.exists(xml_dir):
                logger.info('')
                logger.critical(f"{xml_dir} doesn't exist. Event detection has not "
                                "been run or an incorrect event type has been selected.")
                logger.info('Check documentation for how to run a pipeline:')
                logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
                logger.info('-' * 10)
                return
    
            # Set channels
            chan, ref_chan = check_chans(self.rootpath, chan, ref_chan, logger)
            
            self.track(subs, sessions, step = 'spindle', show = False, 
                       log = False)
            self.tracking['cluster'] = self.tracking['spindle']
            
            # Run analysis
            CLAM = clam(self.rootpath, self.datapath, xml_dir, out_dir, 
                        chan, grp_name, stage, 
                        rater, subs, sessions, self.tracking)
            
            CLAM.clustering(evt_name, freq_bands, filetype, grp_name, logger )
        
        return
    
    
    #--------------------------------------------------------------------------
    '''
    DATASET CREATION.
    
    export_macro_stats -> Exports sleep macroarchitecture per participant into 
                            the corresponding folder in output directory 'staging' 
    
    macro_dataset -> Creates a cohort dataset of sleep macroarchitecture and saves
                        it to a single .csv file in output directory 'dataset'
    
    export_eventparams -> Exports descriptives for sleep events per participant into 
                            the corresponding folder in output directory 'staging'
    
    event_dataset -> Creates a cohort dataset of sleep events descriptives and saves
                        it to a single .csv file in output directory 'dataset'
    
    '''    
    
    def export_macro_stats(self, xml_dir = None, out_dir = None, 
                                 subs = 'all', sessions = 'all', 
                                 times = None, rater = None, 
                                 arousal_name = ['Arousal', 'Arou'], 
                                 outfile = True):
        
        # Set up logging
        if outfile == True:
            today = date.today().strftime("%Y%m%d")
            now = datetime.now().strftime("%H:%M:%S")
            logfile = f'{self.log_dir}/export_sleep_macro_stats_{today}_log.txt'
            logger = create_logger_outfile(logfile=logfile, name='Export macro stats')
            logger.info('')
            logger.info(f"-------------- New call of 'Export macro stats' evoked at {now} --------------")
        elif outfile:
            logfile = f'{self.log_dir}/{outfile}'
            logger = create_logger_outfile(logfile=logfile, name='Export macro stats')
        else:
            logger = create_logger('Export macro stats')

        
        # Set input/output directories
        xml_dir = select_input_dirs(self.outpath, xml_dir, 'macro')
        logger.debug(f'Input annotations being read from: {xml_dir}')
        
        out_dir = select_output_dirs(self.outpath, out_dir, 'macro')
        logger.debug(f'Output being saved to: {out_dir}')
        
        # Set channels
        times, ref_chan = check_chans(self.rootpath, None, True, logger)
        
        self.track(subs = subs, ses = sessions, step = ['staging'], show = False, 
                   log = True)
        
        sleepstats.export_sleepstats(xml_dir, out_dir, subs, sessions, 
                                     rater, times, arousal_name, logger)
        return
    
    def macro_dataset(self, xml_dir = None, out_dir = None, 
                      subs = 'all', sessions = 'all', outfile = True):
         
         # Set input/output directories
         if outfile == True:
             today = date.today().strftime("%Y%m%d")
             now = datetime.now().strftime("%H:%M:%S")
             logfile = f'{self.log_dir}/export_sleep_macro_stats_{today}_log.txt'
             logger = create_logger_outfile(logfile=logfile, name='Export macro stats')
             logger.info('')
             logger.info(f"-------------- New call of 'Macro dataset' evoked at {now} --------------")
         elif outfile:
             logfile = f'{self.log_dir}/{outfile}'
             logger = create_logger_outfile(logfile=logfile, name='Export macro stats')
         else:
             logger = create_logger('Export macro stats')
         if not path.exists(self.outpath + '/datasets/'):
             mkdir(self.outpath + '/datasets/')
         out_dir = self.outpath + '/datasets/macro/'
         
         xml_dir = select_input_dirs(self.outpath, xml_dir, 'macro')
         logger.debug(f'Input annotations being read from: {xml_dir}')
         
         out_dir = select_output_dirs(self.outpath, out_dir, 'macro')
         logger.debug(f'Output being save to: {out_dir}')
         
         sleepstats.sleepstats_from_csvs(xml_dir, out_dir,   
                                 subs, sessions, logger)
         return
    
    def export_eventparams(self, evt_name, frequency = None,
                                 xml_dir = None, out_dir = None, subs = 'all', 
                                 sessions = 'all', filetype = '.edf',
                                 chan = None, ref_chan = None, 
                                 stage = ['NREM2','NREM3'], grp_name = 'eeg', 
                                 rater = None, cycle_idx = None, 
                                 concat_cycle = True, concat_stage = False, 
                                 keyword = None, segs = None,  
                                 adap_bands = 'Fixed',  
                                 adap_bw = 4, params = 'all', epoch_dur = 30, 
                                 average_channels = False, outfile = True):
        
        # Force evt_name into list, and loop through events    
        if isinstance(evt_name, str):
            evts = [evt_name]
        elif isinstance(evt_name, list):
            evts = evt_name
        else:
            raise TypeError(f"'evt_name' can only be a str or a list, but {type(evt_name)} was passed.")
        
        # Set up logging
        if outfile == True:
            subs_str, ses_str = out_names(subs, sessions)
            evt_out = '_'.join(evt_name)
            today = date.today().strftime("%Y%m%d")
            now = datetime.now().strftime("%H:%M:%S")
            logfile = f'{self.log_dir}/export_params_{evt_out}_subs-{subs_str}_ses-{ses_str}_{today}_{now}_log.txt'
            logger = create_logger_outfile(logfile=logfile, name='Export params')
            logger.info('')
            logger.debug(f"-------------- New call of 'Export params' evoked at {now} --------------")
        elif outfile:
            logfile = f'{self.log_dir}/{outfile}'
            logger = create_logger_outfile(logfile=logfile, name='Export params')
        else:
            logger = create_logger('Export params')
            
            
        frequency, flag = adap_bands_setup(self, adap_bands, frequency, 
                                                 subs, sessions, 
                                                 chan, ref_chan, stage, None, 
                                                 concat_cycle, cycle_idx, logger)
        if flag == 'error':
            return
        
        # Set input/output directories
        in_dir = self.datapath
        
        for evt_name in evts:
            
            out_dir = select_output_dirs(self.outpath, out_dir, evt_name)
            logger.debug(f'Output being save to: {out_dir}')
            
            xml_dir = select_input_dirs(self.outpath, xml_dir, evt_name)
            logger.debug(f'Input annotations being read from: {xml_dir}')
            
            # Check annotations directory exists
            if not path.exists(xml_dir):
                logger.info('')
                logger.critical(f"{xml_dir} doesn't exist. Event detection has not "
                                "been run or an incorrect event type has been selected.")
                logger.info('Check documentation for how to run a pipeline:')
                logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
                logger.info('-' * 10)
                return
    
            if adap_bands in ['Auto','Manual']:
                evt_name = f'{evt_name}_adap'
                self.track(step='fooof', show = False, log = False)
                peaks = check_chans(self.rootpath, None, False, logger)
            else:
                peaks = None
            
            # Set channels
            chan, ref_chan = check_chans(self.rootpath, chan, ref_chan, logger)
            if average_channels:
                Ngo = {'run':True}
            else:
                Ngo = {'run':False}
            
            # Format concatenation
            cat = (int(concat_cycle),int(concat_stage),1,1)
            
            # Run line
            fish = FISH(self.rootpath, in_dir, xml_dir, out_dir, chan, ref_chan, grp_name, 
                              stage, rater, subs, sessions, self.tracking) 
            fish.line(keyword, evt_name, cat, segs, cycle_idx, frequency, adap_bands, 
                      peaks, adap_bw, params, epoch_dur, Ngo, filetype, logger)
        return
    
    
    def event_dataset(self, chan, evt_name, xml_dir = None, out_dir = None, 
                            subs = 'all', sessions = 'all', 
                            stage = ['NREM2','NREM3'], concat_stage = False, 
                            concat_cycle = True, cycle_idx = None, 
                            grp_name = 'eeg', adap_bands = 'Fixed',  
                            params = 'all', outfile = True):
        
        # Set up logging
        if outfile == True:
            today = date.today().strftime("%Y%m%d")
            now = datetime.now().strftime("%H:%M:%S")
            logfile = f'{self.log_dir}/event_dataset_{evt_name}_subs-{subs}_ses-{sessions}_{today}_{now}_log.txt'
            logger = create_logger_outfile(logfile=logfile, name='Event dataset')
            logger.info('')
            logger.info(f'-------------- New call of Event dataset evoked at {now} --------------')
        elif outfile:
            logfile = f'{self.log_dir}/{outfile}'
            logger = create_logger_outfile(logfile=logfile, name='Event dataset')
        else:
            logger = create_logger('Event dataset')
        logger = create_logger('Event dataset')
        
        # Force evt_name into list, and loop through events    
        if isinstance(evt_name, str):
            evts = [evt_name]
        elif isinstance(evt_name, list):
            evts = evt_name
        else:
            logger.error(TypeError("'evt_name' can only be a str or a list of str, "
                                   f"but {type(evt_name)} was passed."))
            logger.info('Check documentation for how to create an event_dataset:')
            logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
            logger.info('-' * 10)
            return
        
        frequency, flag = adap_bands_setup(self, adap_bands, None, None, 
                                           None, None, None, stage, None, 
                                           concat_cycle, cycle_idx, logger)
        if flag == 'error':
            return
        
        for evt_name in evts:
            # Append 'adap' after event name if adapted bands were used
            if adap_bands in ['Auto', 'Manual']:
                evt_name = f'{evt_name}_adap'
                self.track(step='fooof', show = False, log = False)
            
            # Set input/output directories
            in_dir = self.datapath
            if not out_dir:    
                if not path.exists(self.outpath + '/datasets/'):
                    mkdir(self.outpath + '/datasets/')
                outpath = self.outpath + f'/datasets/{evt_name}'
            else:
                outpath = out_dir
            if not path.exists(outpath):
                mkdir(outpath)
            logger.debug(f'Output being saved to: {outpath}')
            
            xml_dir = select_input_dirs(self.outpath, xml_dir, evt_name)
            logger.debug(f'Input annotations being read from: {xml_dir}')
            if not path.exists(xml_dir):
                logger.info('')
                logger.critical(f"{xml_dir} doesn't exist. Event detection has not "
                                "been run or an incorrect event type has been selected.")
                logger.info('Check documentation for how to run a pipeline:')
                logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
                logger.info('-' * 10)
                return
            
            # Format concatenation
            cat = (int(concat_cycle),int(concat_stage),1,1)
            
            # Format chan
            if isinstance(chan, str):
                chan = [chan]
            
        
            fish = FISH(self.rootpath, in_dir, xml_dir, outpath, chan, None, grp_name, 
                              stage, subs = subs, sessions = sessions) 
            fish.net(chan, evt_name, adap_bands, params,  cat, cycle_idx, logger)
        
        return

    
    def cluster_peak_dataset(self, chan, evt_name, xml_dir = None, out_dir = None,
                                   subs = 'all', sessions = 'all', 
                                   stage = ['NREM2','NREM3'], concat_stage = False, 
                                   concat_cycle = True, cycle_idx = None, 
                                   grp_name = 'eeg', adap_bands = 'Fixed',  
                                   params = 'all', bins = 60, smooth_sigma = 1.5,
                                   cluster_labels = ('low_peak', 'high_peak'),
                                   outfile = True):
        """
        Apply peak clustering to exported event parameter CSVs, then create
        event datasets per peak cluster (mirroring event_dataset outputs).
        """

        # Set up logging
        if outfile == True:
            today = date.today().strftime("%Y%m%d")
            now = datetime.now().strftime("%H:%M:%S")
            logfile = f'{self.log_dir}/cluster_peak_dataset_{evt_name}_subs-{subs}_ses-{sessions}_{today}_{now}_log.txt'
            logger = create_logger_outfile(logfile=logfile, name='Cluster peak dataset')
            logger.info('')
            logger.info(f'-------------- New call of Cluster peak dataset evoked at {now} --------------')
        elif outfile:
            logfile = f'{self.log_dir}/{outfile}'
            logger = create_logger_outfile(logfile=logfile, name='Cluster peak dataset')
        else:
            logger = create_logger('Cluster peak dataset')
        logger = create_logger('Cluster peak dataset')

        # Force evt_name into list
        if isinstance(evt_name, str):
            evts = [evt_name]
        elif isinstance(evt_name, list):
            evts = evt_name
        else:
            logger.error(TypeError("'evt_name' can only be a str or a list of str, "
                                   f"but {type(evt_name)} was passed."))
            return

        # Format concatenation
        cat = (int(concat_cycle), int(concat_stage), 1, 1)

        # Format channels
        if isinstance(chan, str):
            chan = [chan]

        for evt in evts:
            # Determine input/output directories
            xml_dir_evt = select_input_dirs(self.outpath, xml_dir, evt)
            if not path.exists(xml_dir_evt):
                logger.warning(f"{xml_dir_evt} not found for {evt}; skipping.")
                continue

            if not out_dir:
                if not path.exists(self.outpath + '/datasets/'):
                    mkdir(self.outpath + '/datasets/')
                out_dir_evt = self.outpath + f'/datasets/{evt}'
            else:
                out_dir_evt = out_dir
            if not path.exists(out_dir_evt):
                mkdir(out_dir_evt)

            logger.debug(f'Clustering peak frequencies in: {xml_dir_evt}')
            cluster_peaks(xml_dir_evt, bins=bins, smooth_sigma=smooth_sigma)

            cluster_root = Path(self.outpath) / "cluster_peaks" / evt
            cluster_root.mkdir(parents=True, exist_ok=True)

            for cluster_label in cluster_labels:
                logger.debug(f'Building cluster-specific files for {evt} ({cluster_label})')
                cluster_dir = cluster_root / cluster_label
                _build_cluster_param_tree(xml_dir_evt, cluster_dir, evt, cluster_label, logger)

                fish = FISH(self.rootpath, self.datapath, str(cluster_dir), out_dir_evt,
                            chan, None, grp_name, stage, subs=subs, sessions=sessions)
                fish.net(chan, f"{evt}_{cluster_label}", adap_bands, params, cat, cycle_idx, logger)

        return
    
    
    def pac_dataset(self, chan, evt_name = None, subs = 'all', sessions = 'all',
                          xml_dir = None, out_dir = None,  stage = None, 
                          concat_stage = False, concat_cycle = True, 
                          cycle_idx = None, grp_name = 'eeg', 
                          adap_bands_phase = 'Fixed', frequency_phase = (0.5, 1.25), 
                          adap_bands_amplitude = 'Fixed', frequency_amplitude = (11, 16),  
                          params = 'all', outfile=True):
        
        # Set up logging
        if outfile == True:
            today = date.today().strftime("%Y%m%d")
            now = datetime.now().strftime("%H:%M:%S")
            logfile = f'{self.log_dir}/PAC_dataset_{evt_name}_subs-{subs}_ses-{sessions}_{today}_{now}_log.txt'
            logger = create_logger_outfile(logfile=logfile, name='PAC dataset')
            logger.info('')
            logger.info(f'-------------- New call of PAC dataset evoked at {now} --------------')
        elif outfile:
            logfile = f'{self.log_dir}/{outfile}'
            logger = create_logger_outfile(logfile=logfile, name='PAC dataset')
        else:
            logger = create_logger('PAC dataset')
        logger = create_logger('PAC dataset')
        
        # Set input/output directories
        in_dir = self.datapath

        if not out_dir:
            if not path.exists(self.outpath + '/datasets/'):
                mkdir(self.outpath + '/datasets/')
            out_dir = (f'{self.outpath}/datasets/event_pac' if evt_name 
                       else f'{self.outpath}/datasets/pac') 
        if not path.exists(out_dir):
            mkdir(out_dir)
        logger.debug(f'Output being saved to: {out_dir}')
        
        # Check if event-based or continuous PAC
        if isinstance(evt_name, str):
            evt = 'event_pac' 
            xml_dir = select_input_dirs(self.outpath, xml_dir, evt)
            cat = (int(concat_cycle),int(concat_stage),1,1)
        elif evt_name == None:
            evt = 'pac' 
            xml_dir = select_input_dirs(self.outpath, xml_dir, evt)
            cat = (int(concat_cycle),int(concat_stage),0,0)
        else:
            logger.error(TypeError(f"'evt_name' can only be a str or NoneType, but {type(evt_name)} was passed."))
            logger.info('Check documentation for how to create a PAC summary dataset:')
            logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
            logger.info('-' * 10)
            return
        
        logger.debug(f'Input annotations being read from: {xml_dir}')
        if not path.exists(xml_dir):
            logger.info('')
            logger.critical(f"{xml_dir} doesn't exist. PAC detection has not been "
                            "run or an incorrect type has been selected.")
            logger.info('Check documentation for how to run a pipeline:')
            logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
            logger.info('-' * 10)
            return

        # Format chan
        if isinstance(chan, str):
            chan = [chan]
        
        # Default stage
        if stage == None:
            stage = ['NREM2','NREM3']
        
        # Run extraction
        self.track(subs, sessions, step = ['fooof','spindle'], show = False, 
                   log = False)    
        
        fish = FISH(self.rootpath, in_dir, xml_dir, out_dir, chan, None, grp_name, 
                          stage, None, subs, sessions, self.tracking) 
                   
        fish.pac_summary(chan, evt_name, adap_bands_phase, frequency_phase, 
                              adap_bands_amplitude, frequency_amplitude,
                              params = 'all', cat = cat, cycle_idx = None, 
                              logger = logger)
        
        return
    
    
    def cluster_dataset(self, chan, evt_name, xml_dir = None, out_dir = None, 
                                subs = 'all', sessions = 'all', 
                                stage = ['NREM2'], freq_bands = ('SWA', 'Sigma'),
                                params = 'all', outfile=True):
        
        # Set up logging
        if outfile == True:
            today = date.today().strftime("%Y%m%d")
            now = datetime.now().strftime("%H:%M:%S")
            logfile = f'{self.log_dir}/cluster_dataset_{evt_name}_subs-{subs}_ses-{sessions}_{today}_{now}_log.txt'
            logger = create_logger_outfile(logfile=logfile, name='Event dataset')
            logger.info('')
            logger.info(f'-------------- New call of Cluster dataset evoked at {now} --------------')
        elif outfile:
            logfile = f'{self.log_dir}/{outfile}'
            logger = create_logger_outfile(logfile=logfile, name='Event dataset')
        else:
            logger = create_logger('Event dataset')
        logger = create_logger('Event dataset')
        
        # Force evt_name into list, and loop through events    
        if isinstance(evt_name, str):
            evts = [evt_name]
        elif isinstance(evt_name, list):
            evts = evt_name
        else:
            logger.error(TypeError("'evt_name' can only be a str or a list of str, "
                                   f"but {type(evt_name)} was passed."))
            logger.info('Check documentation for how to create an event_dataset:')
            logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
            logger.info('-' * 10)
            return
        
        for evt_name in evts:    
            # Set input/output directories
            in_dir = self.datapath
            if not out_dir:    
                if not path.exists(self.outpath + '/datasets/'):
                    mkdir(self.outpath + '/datasets/')
                outpath = self.outpath + '/datasets/cluster_fluc'
                if not path.exists(outpath):
                    mkdir(outpath)
                outpath = outpath + f'/{evt_name}'
            else:
                outpath = out_dir
            if not path.exists(outpath):
                mkdir(outpath)
            logger.debug(f'Output being saved to: {outpath}')
            
            xml_dir = select_input_dirs(self.outpath, xml_dir, 'cluster')
            logger.debug(f'Input annotations being read from: {xml_dir}')
            if not path.exists(xml_dir):
                logger.info('')
                logger.critical(f"{xml_dir} doesn't exist. Cluster detection has not "
                                "been run.")
                logger.info('Check documentation for how to run a pipeline:')
                logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
                logger.info('-' * 10)
                return
            
            # Format chan
            if isinstance(chan, str):
                chan = [chan]
            
            fish = FISH(self.rootpath, in_dir, xml_dir, outpath, chan, None, 'eeg', 
                        stage, subs, sessions) 
            
            fish.scalops(chan, evt_name, freq_bands, params, logger)
        
        return
    

    def powerspec_dataset(self, chan, xml_dir = None, out_dir = None, 
                                subs = 'all', sessions = 'all', 
                                stage = ['NREM1','NREM2','NREM3', 'REM'], 
                                concat_stage = False, concat_cycle = True, 
                                cycle_idx = None, grp_name = 'eeg', 
                                rater = None, params = 'all', 
                                general_opts = None, frequency_opts = None, 
                                filter_opts = None, epoch_opts = None, 
                                event_opts = None, outfile=True):
        
        # Set up logging
        if outfile == True:
            today = date.today().strftime("%Y%m%d")
            now = datetime.now().strftime("%H:%M:%S")
            logfile = f'{self.log_dir}/powerspec_subs-{subs}_ses-{sessions}_{today}_{now}_log.txt'
            logger = create_logger_outfile(logfile=logfile, name='Power spectrum dataset')
            logger.info('')
            logger.info(f'-------------- New call of Power spectrum dataset evoked at {now} --------------')
        elif outfile:
            logfile = f'{self.log_dir}/{outfile}'
            logger = create_logger_outfile(logfile=logfile, name='Power spectrum dataset')
        else:
            logger = create_logger('Power spectrum dataset')
        logger = create_logger('Power spectrum dataset')
        
        # Set input/output directories
        in_dir = self.datapath
        log_dir = self.outpath + '/audit/logs/'
        if not path.exists(log_dir):
            mkdir(log_dir)
        if not out_dir:
            if not path.exists(self.outpath + '/datasets/'):
                mkdir(self.outpath + '/datasets/')
            out_dir = f'{self.outpath}/datasets/powerspectrum'
            if not path.exists(out_dir):
                mkdir(out_dir)
        logger.debug(f'Output being saved to: {out_dir}')
        
        xml_dir = select_input_dirs(self.outpath, xml_dir, evt_name = 'powerspectrum')
        logger.debug(f'Input annotations being read from: {xml_dir}')
        if not path.exists(xml_dir):
            logger.info('')
            logger.critical(f"{xml_dir} doesn't exist. Event detection has not "
                            "been run or an incorrect event type has been selected.")
            logger.info('Check documentation for how to run a pipeline:')
            logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
            logger.info('-' * 10)
            return
        
        # Set default parameters
        if not general_opts:
            general_opts = default_general_opts()
        if not frequency_opts:
            frequency_opts = default_frequency_opts()
        if not epoch_opts:
            epoch_opts = default_epoch_opts()  
        if not event_opts:
            event_opts = default_event_opts()
        
        if not filter_opts:
            filter_opts = default_filter_opts()    
        frequency_opts['frequency'] = (filter_opts['highpass'], filter_opts['lowpass'])
        
        # Set suffix for output filename
        if not general_opts['suffix']:
            general_opts['suffix'] = f"{frequency_opts['frequency'][0]}-{frequency_opts['frequency'][1]}Hz"
        
        # Format chan
        if isinstance(chan,str):
            chan = [chan]
        
        # Format concatenation
        cat = (int(concat_cycle),int(concat_stage),1,1)
                            
        spectrum = Spectrum(in_dir, xml_dir, out_dir, chan, None, grp_name, 
                            stage,cat, rater, cycle_idx, subs, sessions) 
        spectrum.powerspec_summary(chan, general_opts, frequency_opts, filter_opts, 
                                   epoch_opts, event_opts, logger)
        
        return

    def bandpower_timecourse(self, band, xml_dir = None, out_dir = None,
                                   subs = 'all', sessions = 'all',
                                   filetype = '.edf', chan = None,
                                   ref_chan = None, grp_name = 'eeg',
                                   rater = None,
                                   stage = ['NREM1','NREM2','NREM3', 'REM'],
                                   cycle_idx = None, concat_cycle = True,
                                   concat_stage = True,
                                   general_opts = None, filter_opts = None,
                                   epoch_opts = None, event_opts = None,
                                   bandpower_opts = None, outfile = True):
        
        """Convenience wrapper that exports band-limited power vs. time."""

        # Local import avoids touching module-level imports per request
        from seapipe.spectrum.bandpower import (
            BandPowerTimecourse,
            default_bandpower_opts,
        )

        # Logging
        if outfile is True:
            subs_str, ses_str = out_names(subs, sessions)
            today = date.today().strftime("%Y%m%d")
            now = datetime.now().strftime("%H%M%S")
            logfile = (f'{self.log_dir}/bandpower_timecourse_'
                       f'subs-{subs_str}_ses-{ses_str}_{today}_{now}_log.txt')
            logger = create_logger_outfile(logfile=logfile,
                                           name='Bandpower timecourse')
            logger.info('')
            logger.info("-------------- New call of 'Bandpower timecourse' "
                        f"evoked at {now} --------------")
        elif outfile:
            logfile = f'{self.log_dir}/{outfile}'
            logger = create_logger_outfile(logfile=logfile,
                                           name='Bandpower timecourse')
        else:
            logger = create_logger('Bandpower timecourse')
        logger.info('')

        # Directories
        in_dir = self.datapath
        xml_dir = select_input_dirs(self.outpath, xml_dir, evt_name = 'staging')
        logger.debug(f'Input annotations being read from: {xml_dir}')
        if not path.exists(xml_dir):
            logger.info('')
            logger.critical(f"{xml_dir} doesn't exist. Sleep staging has not "
                            "been run or hasn't been converted correctly.")
            logger.info('Check documentation for how to set up staging data:')
            logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
            logger.info('-' * 10)
            return

        if not out_dir:
            out_dir = f'{self.outpath}/bandpower_timecourse'
        if not path.exists(out_dir):
            mkdir(out_dir)
        logger.debug(f'Output being saved to: {out_dir}')

        # Channels
        chan, ref_chan = check_chans(self.rootpath, chan, ref_chan, logger)
        if isinstance(chan, str):
            return

        # Defaults
        if general_opts is None:
            general_opts = default_general_opts()
        if filter_opts is None:
            filter_opts = default_filter_opts()
        if epoch_opts is None:
            epoch_opts = default_epoch_opts()
        if event_opts is None:
            event_opts = default_event_opts()
        if bandpower_opts is None:
            bandpower_opts = default_bandpower_opts()

        # Concatenation mask for Spectrum
        cat = (int(concat_cycle),
               int(concat_stage),
               1,
               int(event_opts['concat_events']))

        tracker = BandPowerTimecourse(in_dir, xml_dir, out_dir, chan, ref_chan,
                                      grp_name, stage, cat, rater, cycle_idx,
                                      subs, sessions, self.tracking)
        tracker.bandpower_timecourse(
            band=band,
            general_opts=general_opts,
            filter_opts=filter_opts,
            epoch_opts=epoch_opts,
            bandpower_opts=bandpower_opts,
            event_opts=event_opts,
            filetype=filetype,
            split_by_stage=not concat_stage,
            logger=logger,
        )

        return


def out_names(subs, sessions):
    if isinstance(subs, list):
        subs_str = "_".join(subs).replace('\n', '').replace('\r', '').replace('sub','')
    else:
        subs_str = subs
    if isinstance(sessions, list):
        ses_str = "_".join(sessions).replace('\n', '').replace('\r', '').replace('ses','')
    else:
        ses_str = sessions     
        
    return subs_str, ses_str


def _safe_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _load_cluster_events(csv_path, target_column="Peak power frequency (Hz)", cluster_column="Peak cluster"):
    header = None
    freq_idx = -1
    cluster_idx = None
    events = []
    count_meta = None
    density_meta = None

    with open(csv_path, newline="") as f:
        reader = csv.reader(f)
        for row in reader:
            if not row:
                continue
            if header is None:
                if row[0].strip() == "Count" and len(row) > 1:
                    count_meta = _safe_float(row[1])
                    continue
                if row[0].strip() == "Density" and len(row) > 1:
                    density_meta = _safe_float(row[1])
                    continue
                if target_column in row:
                    header = row
                    freq_idx = row.index(target_column)
                    if cluster_column in row:
                        cluster_idx = row.index(cluster_column)
                continue

            if len(row) <= freq_idx:
                continue
            if not row[0].strip().isdigit():
                continue
            events.append(row)

    return header, cluster_idx, events, count_meta, density_meta


def _summarize_cluster_events(header, cluster_idx, events, cluster_label, count_meta, density_meta):
    if header is None:
        raise ValueError("Missing header row in CSV.")
    if cluster_idx is None:
        raise ValueError("Missing Peak cluster column; run cluster_peaks first.")

    mapping = {name: idx for idx, name in enumerate(header)}
    filtered = [r for r in events if len(r) > cluster_idx and r[cluster_idx] == cluster_label]

    count = len(filtered)
    total_time = None
    if count_meta is not None and density_meta not in (None, 0):
        total_time = count_meta / density_meta
    density = count / total_time if total_time else np.nan

    def values_for(col):
        idx = mapping.get(col)
        vals = []
        if idx is None:
            return vals
        for row in filtered:
            val = _safe_float(row[idx]) if len(row) > idx else None
            if val is not None:
                vals.append(val)
        return vals

    def mean_std(vals):
        if not vals:
            return (np.nan, np.nan)
        arr = np.asarray(vals, dtype=float)
        mean = float(np.mean(arr))
        std = float(np.std(arr, ddof=1)) if len(arr) > 1 else np.nan
        return (mean, std)

    metrics = {}
    for col in [
        "Duration (s)",
        "Min. amplitude (uV)",
        "Max. amplitude (uV)",
        "Peak-to-peak amplitude (uV)",
        "Power (uV^2)",
        "Peak power frequency (Hz)",
    ]:
        metrics[col] = mean_std(values_for(col))

    summary = {
        "count": count,
        "density": density,
        "metrics": metrics,
    }
    return summary, filtered


def _write_cluster_summary_csv(out_path, header, cluster_idx, cluster_label, summary, filtered):
    cluster_column = "Peak cluster"
    if cluster_idx is None:
        header = header + [cluster_column]
        cluster_idx = len(header) - 1

    out_path.parent.mkdir(parents=True, exist_ok=True)
    writer_rows = []
    writer_rows.append(["Count", summary["count"]])
    writer_rows.append(["Density", summary["density"]])
    writer_rows.append(header)

    mean_row = ["Mean"] + [""] * (len(header) - 1)
    sd_row = ["SD"] + [""] * (len(header) - 1)
    for col, (mean, std) in summary["metrics"].items():
        idx = header.index(col) if col in header else None
        if idx is not None:
            mean_row[idx] = mean
            sd_row[idx] = std
    writer_rows.append(mean_row)
    writer_rows.append(sd_row)

    for row in filtered:
        row = list(row)
        if len(row) <= cluster_idx:
            row = row + [""] * (cluster_idx - len(row) + 1)
        row[cluster_idx] = cluster_label
        if len(row) < len(header):
            row = row + [""] * (len(header) - len(row))
        writer_rows.append(row)

    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(writer_rows)


def _build_cluster_param_tree(src_root, dest_root, evt_name, cluster_label, logger):
    src_root = Path(src_root)
    dest_root = Path(dest_root)
    for csv_path in src_root.rglob("*.csv"):
        if not csv_path.stem.endswith(evt_name):
            continue
        header, cluster_idx, events, count_meta, density_meta = _load_cluster_events(csv_path)
        if header is None:
            logger.warning(f"Skipping {csv_path}: no header found.")
            continue
        if not events:
            logger.warning(f"Skipping {csv_path}: no event rows found.")
            continue

        try:
            summary, filtered = _summarize_cluster_events(header, cluster_idx, events, cluster_label, count_meta, density_meta)
        except Exception as exc:  # noqa: BLE001
            logger.warning(f"Skipping {csv_path}: {exc}")
            continue

        rel_parent = csv_path.relative_to(src_root).parent
        if csv_path.stem.endswith(evt_name):
            base = csv_path.stem[: -len(evt_name)]
        else:
            base = f"{csv_path.stem}_"
        new_stem = f"{base}{evt_name}_{cluster_label}"
        dest_file = dest_root / rel_parent / f"{new_stem}{csv_path.suffix}"
        _write_cluster_summary_csv(dest_file, header, cluster_idx, cluster_label, summary, filtered)
        
