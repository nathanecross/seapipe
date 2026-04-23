#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 12:07:36 2023

@author: nathancross
"""
from datetime import datetime, date
from os import listdir, mkdir, path, remove, walk, rename
import shutil
from pathlib import Path
import numpy as np
import csv
from pandas import DataFrame



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
        
        from seapipe.utils.audit import (check_dataset)
        
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
        
        from seapipe.utils.audit import (check_dataset)
        from seapipe.utils.logs import create_logger, create_logger_outfile
        
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
        
        from seapipe.utils.logs import create_logger, create_logger_outfile

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
        
        from seapipe.utils.audit import track_processing
        from seapipe.utils.load import read_tracking_sheet
        
        ## Set up logging
        logger = setup_logging(self.log_dir, 'Tracking', outfile=True)
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
            except Exception:
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
            except Exception:
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
        except Exception:
            self.tracking = {**self.tracking, **tracking}
        
        if show:
            logger.info('')
            logger.info(df)
        if outfile:
            df.to_csv(f'{self.outpath}/audit/{outfile}')

        return   

    def make_bids(self, subs = 'all', origin = 'SCN', filetype = '.edf',
                  outfile = True):
        
        from seapipe.utils.audit import make_bids
        
        # Set up logging
        logger = setup_logging(self.log_dir, 'Make bids', outfile)
        logger.info('')
        logger.debug('Formatting dataset into BIDS.')
        
        make_bids(self.datapath, subs, origin, filetype)
        
    def extract_channels(self, exclude = None):
        
        from seapipe.utils.audit import extract_channels
        
        extract_channels(self.datapath, exclude)
        
    def QC_channels(self, subs = 'all', sessions = 'all', filetype = '.edf',
                    filt = None, chantype = None, 
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
        
        from seapipe.utils.squid import SQUID
        
        # Set up logging
        logger = setup_logging(self.log_dir, 'Channel QC', outfile)
        logger.info('')
        
        # Check chantypes
        chantype = chantype if chantype is not None else ['eeg', 'eog', 'emg', 'ecg']
        valid_chantypes = {'eeg', 'eog', 'emg', 'ecg'}
        if not set(chantype).issubset(valid_chantypes):
            raise ValueError(f"Invalid chantype(s) specified. Allowed: {valid_chantypes}")
        
        # Run detection
        if subs == 'all':
            subs = None
        if sessions == 'all':
            sessions = None
        squid = SQUID(self.rootpath, filetype = filetype, subjects = subs, 
                     sessions = sessions, logger = logger)
        squid.process_all(filt, chantype)
        
    def QC_summary(self, qc_dir = None, chantype = None):
        
        from seapipe.utils.squid import gather_qc_reports
        
        chantype = chantype if chantype is not None else ['eeg', 'eog', 'emg', 'ecg']
        
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
        
        from seapipe.utils.load import load_stages
        
        # Set up logging
        logger = setup_logging(self.log_dir, 'Load sleep stages', outfile)
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
                             stage = None, 
                             cycle_idx = None, concat_cycle = True, 
                             concat_stage = False, general_opts = None, 
                             frequency_opts = None, filter_opts = None, 
                             epoch_opts = None, event_opts = None, 
                             norm = None, norm_opts = None, 
                             outfile = True):
        
        from seapipe.spectrum.psa import (Spectrum, default_epoch_opts, 
                                          default_event_opts,
                                          default_filter_opts, 
                                          default_frequency_opts, 
                                          default_general_opts, 
                                          default_norm_opts)
        
        from seapipe.utils.load import select_input_dirs, check_chans
        
        # Set up logging
        logger = setup_logging(self.log_dir, 'Power spectrum', outfile)
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
        
        if not stage:
            stage = ['NREM1','NREM2','NREM3', 'REM']
        
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
        except Exception:
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
        from seapipe.utils.load import check_chans, read_tracking_sheet

        # Set up logging
        logger = setup_logging(self.log_dir, 'Detect sleep stages', outfile)
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
        except Exception:
            self.tracking = {**self.tracking, **stages.tracking}
        return
    
    
    def detect_artefacts(self, xml_dir = None, out_dir = None, 
                               subs = 'all', sessions = 'all', filetype = '.edf', 
                               method = 'seapipe', win_size = 5,
                               chan = None, ref_chan = None, 
                               label = 'Artefact', allchans_marker = False,
                               rater = None, grp_name = 'eeg', 
                               stage = None,
                               outfile = True):
        
        from seapipe.events.sand import SAND
        from seapipe.utils.load import (check_chans, read_tracking_sheet, 
                                        select_input_dirs)

        # Set up logging
        logger = setup_logging(self.log_dir, 'Detect artefacts', outfile)
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
        
        # Set stage defaults
        if not stage:
            stage = ['NREM1', 'NREM2', 'NREM3', 'REM']
            
        # Check annotations directory exists, run detection
        artefacts = SAND(in_dir, xml_dir, out_dir, chan, ref_chan, 
                         rater, grp_name, subs, sessions, self.tracking) 
        artefacts.detect_artefacts(method, label,  win_size,  filetype, 
                                   allchans_marker, stage, logger)
                                   
    
        try:
            self.tracking = self.tracking | artefacts.tracking
        except Exception:
            self.tracking = {**self.tracking, **artefacts.tracking}
        return
        
        
    
    def detect_spectral_peaks(self, xml_dir = None, out_dir = None, 
                                    subs = 'all', sessions = 'all', chan = None, 
                                    ref_chan = None, grp_name = 'eeg', 
                                    rater = None, frequency = (9,16), 
                                    stage = None, cycle_idx = None,
                                    concat_cycle = True, concat_stage = False,
                                    general_opts = None, frequency_opts = None,
                                    filter_opts = None, epoch_opts = None, 
                                    event_opts = None, fooof_opts = None, 
                                    filetype = '.edf', suffix = None, 
                                    outfile = True):
        
        from seapipe.spectrum.psa import (Spectrum, default_epoch_opts, 
                                          default_event_opts,
                                          default_fooof_opts,
                                          default_filter_opts, 
                                          default_frequency_opts, 
                                          default_general_opts)
        from seapipe.utils.load import (check_chans, read_tracking_sheet, 
                                        select_input_dirs)
        
        # Set up logging
        logger = setup_logging(self.log_dir, 'Detect spectral peaks', outfile)
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
        
        # Set stage defaults
        if not stage:
            stage = ['NREM2','NREM3']
            
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
                                       stage = None, cycle_idx=None, 
                                       duration=(0.2, 2), invert = None,
                                       reject_artf = None,
                                       average_channels = False, outfile = True):
        
        from seapipe.events.seasnakes import seasnakes
        from seapipe.utils.load import (check_chans, read_tracking_sheet, 
                                        select_input_dirs)

        # Set up logging
        logger = setup_logging(self.log_dir, 'Detect slow oscillations', outfile)
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
        
        # Set stage defaults
        if not stage:
            stage = ['NREM2','NREM3']
        
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
        except Exception:
            self.tracking = {**self.tracking, **SO.tracking}
        
        return
    
    
    def detect_spindles(self, xml_dir = None, out_dir = None, subs = 'all', 
                              sessions = 'all', filetype = '.edf', 
                              method = ['Moelle2011'], chan = None, 
                              ref_chan = None, rater = None, 
                              stage = None, grp_name = 'eeg', 
                              cycle_idx = None, concat_cycle = True, 
                              frequency = None, adap_bands = 'Fixed', 
                              adap_bw = 4, duration = (0.5, 3),
                              reject_artf = None, 
                              outfile = True):
        
        from seapipe.events.whales import whales
        from seapipe.utils.load import (check_chans, read_tracking_sheet, 
                                        select_input_dirs, select_output_dirs)
        from seapipe.utils.misc import adap_bands_setup
        
        
        # Set up logging
        logger = setup_logging(self.log_dir, 'Detect spindles', outfile)
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
        
        # Set stage defaults
        if not stage:
            stage = ['NREM2','NREM3']
        
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
                                           sessions, chan, ref_chan, stage, False, 
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
        except Exception:
            self.tracking = {**self.tracking, **spindle.tracking}
            
        return
    
    
    def whales(self, xml_dir = None, out_dir = None, subs = 'all', 
                     sessions = 'all', filetype = '.edf', 
                     method = None, evt_out = 'spindle',
                     merge_type = 'consensus', weights = None,
                     chan = None, ref_chan = None, rater = None, 
                     stage = None, grp_name = 'eeg', 
                     cycle_idx = None, s_freq = None, keyword = None, 
                     duration =(0.5, 3),
                     reject_artf = None, 
                     outfile = True):
        
        from seapipe.events.whales import whales
        from seapipe.utils.load import (check_chans, read_tracking_sheet, 
                                        select_input_dirs)
        
        # Set up logging
        logger = setup_logging(self.log_dir, 'Detect spindles (WHALES)', outfile)
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
        
        # Set detection methods
        if not method:
            method = ['Moelle2011', 'Ray2015']
        
        # Define consensus threshold based on method
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
        
        # Set stage defaults
        if not stage:
            stage = ['NREM2','NREM3']
        
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
                          method = 'YASA', chan = None,
                          ref_chan = None, rater = None, grp_name = 'eeg', 
                          stage = None, cycle_idx = None, 
                          amplitude = (50, 325), duration = (0.3, 1.5),
                          reject_artf = None,
                          average_channels = False, outfile = True):
        
        from seapipe.events.remora import remora
        from seapipe.utils.load import (check_chans, read_tracking_sheet, 
                                        select_input_dirs)

        # Set up logging
        logger = setup_logging(self.log_dir, 'Detect eye movements (REMS)', outfile)
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
        
        # Set stage defaults
        if not stage:
            stage = ['REM']
    
        # Run detection
        REMS = remora(in_dir, xml_dir, out_dir, chan, ref_chan, 
                      rater, grp_name, reject_artf, subs, sessions)
        
        REMS.detect_rems(method, amplitude, duration, stage, cycle_idx, 
                         filetype, logger)
        try:
            self.tracking = self.tracking | REMS.tracking
        except Exception:
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
        from seapipe.spectrum.psa import default_filter_opts
        from seapipe.utils.load import (check_chans, read_tracking_sheet, 
                                        select_input_dirs)

        # Set up logging
        logger = setup_logging(self.log_dir, 'Event spectrogram', outfile)
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
    
    
    def plot_spectrogram(self, xml_dir=None, out_dir=None, subs='all',
                         sessions='all', filetype='.edf', chan=None,
                         ref_chan=None, rater=None, stage=None,
                         grp_name='eeg', cycle_idx=None, 
                         concat_stage=False, concat_cycle=True, 
                         freq_limits = (0, 25),
                         event_opts=None, win_sec = 30,
                         method='yasa', file_extension='.svg',
                         progress=True, outfile=False):

        from seapipe.utils.spectrogram import SONAR
        from seapipe.spectrum.psa import default_event_opts
        from seapipe.utils.load import (check_chans, read_tracking_sheet, 
                                        select_input_dirs)

        # Set up logging
        logger = setup_logging(self.log_dir, 'Spectrogram', outfile)
        logger.info('')
        logger.debug('Creating spectrogram plots.')
        logger.info('')

        # Set input/output directories
        in_dir = self.datapath

        if not xml_dir:
            xml_dir = select_input_dirs(self.outpath, xml_dir, event_opts['evt_type'])

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
            out_dir = f'{self.outpath}/spectrogram'
        if not path.exists(out_dir):
            mkdir(out_dir)
        logger.debug(f'Output figures being saved to: {out_dir}')

        # Check subs
        if not subs:
            tracking = read_tracking_sheet(self.rootpath, logger)
            subs = [x for x in list(set(tracking['sub']))]
            subs.sort()
        if not sessions:
            sessions = read_tracking_sheet(self.rootpath, logger)

        # Set channels
        chan, ref_chan = check_chans(self.rootpath, chan, ref_chan, logger)

        # Event options
        if not event_opts:
            event_opts = default_event_opts()
        if event_opts['evt_type']:
            win_sec = 1

        sonar = SONAR(self.rootpath, in_dir, xml_dir, out_dir, chan, ref_chan,
                      grp_name, stage, rater, subs, sessions,
                      event_opts=event_opts, method=method,
                      file_extension=file_extension, win_sec = win_sec,
                      fmin=freq_limits[0], fmax=freq_limits[1],
                      tracking=self.tracking)

        sonar.spectrogram(cycle_idx=cycle_idx, concat_stage=concat_stage,
                          concat_cycle=concat_cycle,
                          filetype=filetype, progress=progress,
                          logger=logger)

        return
    #--------------------------------------------------------------------------
    '''
    PHASE AMPLITUDE COUPLING.
    
    mean_amps -> runs Phase Amplitude Coupling analyses on sleep EEG data. 
    
    
    '''    
    def pac(self, xml_dir = None, out_dir = None, subs = 'all', sessions = 'all', 
                  filetype = '.edf', chan = None, ref_chan = None, rater = None, 
                  grp_name = 'eeg', stage = None, concat_stage = True, 
                  cycle_idx = None, concat_cycle = True,  
                  method = 'MI', surrogate = 'Time lag', correction = 'Z-score',
                  adap_bands_phase = 'Fixed', frequency_phase = (0.5, 1.25), 
                  adap_bands_amplitude = 'Fixed', frequency_amplitude = (11, 16),
                  adap_bw = 4, min_dur = 1, nbins = 18, invert = None,
                  frequency_opts = None, filter_opts = None, epoch_opts = None, 
                  evt_name = None, event_opts = None, 
                  reject_artf = None, 
                  progress = True, outfile = True):
        
        from seapipe.pac.octopus import octopus, pac_method
        from seapipe.pac.pacats import pacats
        from seapipe.spectrum.psa import (default_epoch_opts, 
                                          default_event_opts,
                                          default_filter_opts, 
                                          default_frequency_opts)
        from seapipe.utils.load import (check_chans, read_tracking_sheet, 
                                        select_input_dirs)
        from seapipe.utils.misc import adap_bands_setup
        

        # Set up logging
        if outfile == True:
            subs_str, ses_str = out_names(subs, sessions)
            pha = f'{frequency_phase[0]}-{frequency_phase[1]}'
            amp = f'{frequency_amplitude[0]}-{frequency_amplitude[1]}'
            today = date.today().strftime("%Y%m%d")
            now = datetime.now().strftime("%H%M%S")
            outfile = f'pac_{pha}_{amp}_subs-{subs_str}_ses-{ses_str}_{today}_{now}_log.txt'
        logger = setup_logging(self.log_dir, 'Phase-amplitude coupling', outfile)
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
        
        # Set stage defaults
        if not stage:
            stage = ['NREM2', 'NREM3']
        
        # Check for adapted bands 
        cat = (int(concat_cycle),int(concat_stage),0,0)
        
        frequency_phase, flag1 = adap_bands_setup(self, adap_bands_phase, frequency_phase, 
                                     subs, sessions, chan, ref_chan, stage, concat_stage, 
                                     concat_cycle, cycle_idx, logger)
        
        frequency_amplitude, flag2 = adap_bands_setup(self, adap_bands_amplitude, 
                                               frequency_amplitude, subs, sessions, 
                                               chan, ref_chan, stage, concat_stage, 
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

    #--------------------------------------------------------------------------
    '''
    EVENT SYNCHRONY

    event_synchrony -> splits a target event based on co-occurrence with a probe.
    '''
    def event_synchrony(self, xml_dir = None, out_dir = None, subs = 'all',
                              sessions = 'all', stage = None,
                              chan = None, ref_chan = None,
                              grp_name = 'eeg', rater = None,
                              evttype_target = None, evttype_probe = None,
                              evttype_tp_target = None, evttype_fn = None,
                              iu_thresh = 0.5, concat_stage = True,
                              concat_cycle = True,
                              reject_artf = None,
                              filetype = ('.edf', '.rec', '.eeg'),
                              outfile = True):

        from seapipe.pac.coral import CORAL
        from seapipe.utils.load import (check_chans, read_tracking_sheet, 
                                        select_output_dirs, select_input_dirs)

        # Set up logging
        logger = setup_logging(self.log_dir, 'Event synchrony', outfile)
        logger.info('')

        if not evttype_target or not evttype_probe or not evttype_tp_target:
            logger.error("Event synchrony requires 'evttype_target', "
                         "'evttype_probe' and 'evttype_tp_target'.")
            return

        # Set input/output directories
        in_dir = self.datapath
        xml_dir = select_input_dirs(self.outpath, xml_dir, evttype_target)

        if not path.exists(xml_dir):
            logger.info('')
            logger.critical(f"{xml_dir} doesn't exist. Events weren't found.")
            logger.info('Check documentation for how to set up staging data:')
            logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
            logger.info('-' * 10)
            return
        else:
            logger.debug(f'Input annotations being read from: {xml_dir}')

        if not out_dir:
            out_dir = select_output_dirs(self.outpath, out_dir, 'sync')
        if not path.exists(out_dir):
            mkdir(out_dir)
        logger.debug(f'Output annotations being saved to: {out_dir}')
        
        
        # Set channels
        chan, ref_chan = check_chans(self.rootpath, chan, ref_chan, logger)
        if not isinstance(chan, DataFrame) and not isinstance(chan, list):
            return
        elif isinstance(ref_chan, str):
            return
        
        # Check subs
        if not subs:
            tracking = read_tracking_sheet(self.rootpath, logger)
            subs = [x for x in list(set(tracking['sub']))]
            subs.sort()
        if not sessions:
            sessions = read_tracking_sheet(self.rootpath, logger)

        # Concatenation
        cat = (int(concat_cycle), int(concat_stage), 0, 0)

        coral = CORAL(self.rootpath, in_dir, xml_dir, out_dir,
                      chan = chan, ref_chan = ref_chan, stage = stage,
                      grp_name = grp_name, rater = rater, subs = subs,
                      sessions = sessions, reject_artf = reject_artf,
                      filetype = filetype, tracking = self.tracking)

        coral.event_sync(evttype_target, evttype_probe, iu_thresh,
                         evttype_tp_target, evttype_fn, cat = cat,
                         logger = logger)

        return

    def event_synchrony_dataset(self, xml_dir = None, out_dir = None, subs = 'all',
                                      sessions = 'all', chan = None, stage = None,
                                      grp_name = 'eeg', rater = None,
                                      evttype_target = None, evttype_probe = None,
                                      evttype_tp_target = None, evttype_fn = None,
                                      iu_thresh = 0.5, concat_stage = True,
                                      concat_cycle = True, outfile_suffix = None,
                                      reject_artf = None,
                                      filetype = ('.edf', '.rec', '.eeg'),
                                      outfile = True):

        from seapipe.pac.synchrony import CORAL
        from seapipe.utils.load import (read_tracking_sheet, 
                                        select_input_dirs)

        # Set up logging
        logger = setup_logging(self.log_dir, 'Event synchrony dataset (CORAL)', outfile)
        logger.info('')

        if not evttype_target or not evttype_probe or not evttype_tp_target:
            logger.error("Event synchrony dataset requires 'evttype_target', "
                         "'evttype_probe' and 'evttype_tp_target'.")
            return

        # Set input/output directories
        in_dir = self.datapath
        xml_dir = select_input_dirs(self.outpath, xml_dir, evttype_target)

        if not path.exists(xml_dir):
            logger.info('')
            logger.critical(f"{xml_dir} doesn't exist. Events weren't found.")
            logger.info('Check documentation for how to set up staging data:')
            logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
            logger.info('-' * 10)
            return
        else:
            logger.debug(f'Input annotations being read from: {xml_dir}')

        if not out_dir:
            if not path.exists(self.outpath + '/datasets/'):
                mkdir(self.outpath + '/datasets/')
            out_dir = f'{self.outpath}/datasets/sync'
        if not path.exists(out_dir):
            mkdir(out_dir)
        logger.debug(f'Output dataset being saved to: {out_dir}')

        if not outfile_suffix:
            outfile_suffix = f'{evttype_target}_x_{evttype_probe}_sync_stats.csv'

        # Check subs
        if not subs:
            tracking = read_tracking_sheet(self.rootpath, logger)
            subs = [x for x in list(set(tracking['sub']))]
            subs.sort()
        if not sessions:
            sessions = read_tracking_sheet(self.rootpath, logger)

        # Concatenation
        cat = (int(concat_cycle), int(concat_stage), 0, 0)

        coral = CORAL(self.rootpath, in_dir, xml_dir, out_dir,
                      chan = chan, stage = stage, grp_name = grp_name,
                      rater = rater, subs = subs, sessions = sessions,
                      reject_artf = reject_artf, filetype = filetype,
                      tracking = self.tracking)

        coral.event_sync_dataset(outfile_suffix, evttype_target, evttype_probe,
                                 iu_thresh, evttype_tp_target, evttype_fn,
                                 cat = cat, logger = logger)

        return

    '''
    DYNAMICS
    
    spindle_clustering -> analysis of spindle event temporal clustering 
    
    sigma_fluctuations -> 
    '''
    
    def cluster_flucs(self, evt_name, xml_dir = None, out_dir = None,
                            freq_bands = None,
                            subs = 'all', sessions = 'all', 
                            filetype = '.edf',
                            chan = None, ref_chan = None, grp_name = 'eeg', 
                            stage = None, concat_stage = False,
                            spectral_method = 'welch',
                            min_bout_length = 300,
                            allowable_interruptions = 1,
                            rejoin_artefact = None,
                            min_total_nrem_sec = None,
                            min_bouts_psd = None,
                            low_snr_percentile = None,
                            plot_fit = False,
                            rater = None, 
                            outfile = True):
        
        from seapipe.events.clam import clam
        from seapipe.utils.load import (check_chans, 
                                        select_input_dirs, select_output_dirs)

        # Force evt_name into list, and loop through events    
        if isinstance(evt_name, str):
            evts = [evt_name]
        elif isinstance(evt_name, list):
            evts = evt_name
        else:
            raise TypeError(f"'evt_name' can only be a str or a list, but {type(evt_name)} was passed.")
        
        # Define default frequency bands
        if not freq_bands:
            freq_bands = {'SWA': (0.5, 4), 'Sigma': (10, 15)}
        
        # Set up logging
        if outfile == True:
            subs_str, ses_str = out_names(subs, sessions)
            evt_out = '_'.join(evt_name)
            today = date.today().strftime("%Y%m%d")
            now = datetime.now().strftime("%H:%M:%S")
            outfile = f'event_clustering_{evt_out}_subs-{subs_str}_ses-{ses_str}_{today}_{now}_log.txt'
        logger = setup_logging(self.log_dir, 'Event clustering', outfile)
        logger.info('')
        
        # Set channels
        chan, ref_chan = check_chans(self.rootpath, chan, ref_chan, logger)
        
        # Set stage defaults
        if not stage:
            stage = ['NREM2', 'NREM3']
        
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
            
            if not isinstance(chan, DataFrame) and chan == 'error':
                return
            
            self.track(subs, sessions, step = 'spindle', show = False, 
                       log = False)
            self.tracking['cluster'] = self.tracking['spindle']
            
            # Run analysis
            CLAM = clam(self.rootpath, self.datapath, xml_dir, out_dir, 
                        chan, grp_name, stage, 
                        rater, subs, sessions, self.tracking)
            
            CLAM.clustering(evt_name, freq_bands, filetype, grp_name,
                            concat_stage, spectral_method,
                            min_bout_length, allowable_interruptions, rejoin_artefact,
                            min_total_nrem_sec, min_bouts_psd, low_snr_percentile,
                            plot_fit, logger)
        
        return


    def plot_infraslow(self, xml_dir = None, out_dir = None,
                             subs = 'all', sessions = 'all',
                             chan = None, ref_chan = None, grp_name = 'eeg',
                             stage = None, rater = None, seed = 1,
                             outfile = True):

        """
        Create infraslow diagnostic panels (A-C) for each subject/session.

        Uses staging annotations to determine sleep onset, NREM bouts, and the
        normalization window; plots are saved under derivatives/clusterfluc.
        """

        from seapipe.events.clam import clam
        from seapipe.utils.load import (check_chans, select_output_dirs, 
                                        select_input_dirs)

        # Set up logging
        logger = setup_logging(self.log_dir, 'Plot infraslow', outfile)
        logger.info('')

        # Resolve input/output directories
        xml_dir = select_input_dirs(self.outpath, xml_dir, 'staging')
        logger.debug(f'Input annotations being read from: {xml_dir}')

        out_dir = select_output_dirs(self.outpath, out_dir, 'cluster')
        logger.debug(f'Output being saved to: {out_dir}')

        # Set channels (tracking sheet)
        chan, ref_chan = check_chans(self.rootpath, chan, ref_chan, logger)
        if not isinstance(chan, DataFrame) and chan == 'error':
            return
        
        # Set stage defaults
        if not stage:
            stage = ['NREM2', 'NREM3']

        # Instantiate CLAM
        CLAM = clam(self.rootpath, self.datapath, xml_dir, out_dir,
                    chan, grp_name, stage,
                    rater, subs, sessions, self.tracking)

        # Run plotting (looping handled inside CLAM)
        CLAM.plot_infraslow_panels(subs=subs, sessions=sessions,
                                   stage=stage,
                                   xml_dir=xml_dir, out_dir=out_dir,
                                   seed=seed, logger=logger)

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
                                 arousal_name = None, 
                                 outfile = True):
        
        from seapipe.stats import sleepstats
        from seapipe.utils.load import (check_chans,
                                        select_output_dirs, 
                                        select_input_dirs)
        
        # Set up logging
        logger = setup_logging(self.log_dir, 'Export macro stats', outfile)
        logger.info('')

        
        # Set input/output directories
        xml_dir = select_input_dirs(self.outpath, xml_dir, 'macro')
        logger.debug(f'Input annotations being read from: {xml_dir}')
        
        out_dir = select_output_dirs(self.outpath, out_dir, 'macro')
        logger.debug(f'Output being saved to: {out_dir}')
        
        # Set channels
        times, ref_chan = check_chans(self.rootpath, None, True, logger)
        
        # Set Arousal name
        if not arousal_name:
            arousal_name = ['Arousal', 'Arou']
        
        self.track(subs = subs, ses = sessions, step = ['staging'], show = False, 
                   log = True)
        
        sleepstats.export_sleepstats(xml_dir, out_dir, subs, sessions, 
                                     rater, times, arousal_name, logger)
        return
    
    def macro_dataset(self, xml_dir = None, out_dir = None, 
                      subs = 'all', sessions = 'all', outfile = True):
        
        from seapipe.stats import sleepstats
        from seapipe.utils.load import (select_output_dirs, 
                                        select_input_dirs)
        
        # Set up logging
        logger = setup_logging(self.log_dir, 'Export macro datast', outfile)
        logger.info('')

        # Set input/output directories
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

    def tide(self, xml_dir = None, out_dir = None, subject_out_dir = None,
             subs = 'all', sessions = 'all', stage = None, rater = None,
             resolution = 'complete', analyses = 'all', keyword = None,
             outfile = True):

        from seapipe.stats.tide import tide as TIDE
        from seapipe.utils.load import (read_tracking_sheet,
                                        select_input_dirs)

        # Set up logging
        logger = setup_logging(self.log_dir, 'TIDE hypnogram dynamics', outfile)
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
            if not path.exists(self.outpath + '/datasets/'):
                mkdir(self.outpath + '/datasets/')
            out_dir = f'{self.outpath}/datasets/hypnogram'
        if not path.exists(out_dir):
            mkdir(out_dir)
        logger.debug(f'Output group dataset being saved to: {out_dir}')

        if not subject_out_dir:
            subject_out_dir = f'{self.outpath}/hypnogram'
        if not path.exists(subject_out_dir):
            mkdir(subject_out_dir)
        logger.debug(f'Output subject-level files being saved to: {subject_out_dir}')

        # Check subs
        if not subs:
            tracking = read_tracking_sheet(self.rootpath, logger)
            subs = [x for x in list(set(tracking['sub']))]
            subs.sort()
        if not sessions:
            tracking = read_tracking_sheet(self.rootpath, logger)
            sessions = [x for x in list(set(tracking['ses']))]
            sessions.sort()

        # Set stage defaults
        if not stage:
            stage = ['Wake', 'NREM1', 'NREM2', 'NREM3', 'REM']

        if analyses == 'all':
            analyses = ['transition_matrix',
                        'stage_duration_distributions',
                        'hypnogram_similarity']
        elif isinstance(analyses, str):
            analyses = [analyses]

        dynamics = TIDE(xml_dir, out_dir, stage, rater, subs, sessions, keyword,
                        subject_out_dir)

        if 'transition_matrix' in analyses:
            dynamics.transition_matrix(stage = stage, resolution = resolution,
                                       logger = logger)
        if 'stage_duration_distributions' in analyses:
            dynamics.stage_duration_distributions(stage = stage,
                                                  logger = logger)
        if 'hypnogram_similarity' in analyses:
            dynamics.hypnogram_similarity(stage = stage, logger = logger)

        return
    
    def export_eventparams(self, evt_name, frequency = None,
                                 xml_dir = None, out_dir = None, subs = 'all', 
                                 sessions = 'all', filetype = '.edf',
                                 chan = None, ref_chan = None, 
                                 stage = None, grp_name = 'eeg', 
                                 rater = None, cycle_idx = None, 
                                 concat_cycle = True, concat_stage = False, 
                                 keyword = None, segs = None,  
                                 adap_bands = 'Fixed',  
                                 adap_bw = 4, params = 'all', epoch_dur = 30, 
                                 average_channels = False, outfile = True):
        
        from seapipe.events.fish import FISH
        from seapipe.utils.misc import adap_bands_setup
        from seapipe.utils.load import (check_chans,
                                        select_output_dirs, 
                                        select_input_dirs)
        
        # Force evt_name into list, and loop through events    
        if isinstance(evt_name, str):
            evts = [evt_name]
        elif isinstance(evt_name, list):
            evts = evt_name
        else:
            raise TypeError(f"'evt_name' can only be a str or a list, but {type(evt_name)} was passed.")
        
        # Set up logging
        subs_str, ses_str = out_names(subs, sessions)
        evt_out = '_'.join(evt_name)
        today = date.today().strftime("%Y%m%d")
        now = datetime.now().strftime("%H:%M:%S")
        logfile = f'{self.log_dir}/export_params_{evt_out}_subs-{subs_str}_ses-{ses_str}_{today}_{now}_log.txt'
        logger = setup_logging(self.log_dir, 'Export event params', logfile)
        logger.info('')
        
        # Set channels
        chan, ref_chan = check_chans(self.rootpath, chan, ref_chan, logger)
        if not isinstance(chan, DataFrame) and not isinstance(chan, list):
            return
        elif isinstance(ref_chan, str):
            return
        
        # Set stage defaults
        if not stage:
            stage = ['NREM2', 'NREM3']
            
        # Set frequency bands
        frequency, flag = adap_bands_setup(self, adap_bands, frequency, 
                                                 subs, sessions, 
                                                 chan, ref_chan, stage, concat_stage, 
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
                            stage = None, concat_stage = False, 
                            concat_cycle = True, cycle_idx = None, 
                            grp_name = 'eeg', adap_bands = 'Fixed',  
                            params = 'all', outfile = True):
        
        from seapipe.events.fish import FISH
        from seapipe.utils.misc import adap_bands_setup
        from seapipe.utils.load import (select_input_dirs)
        
        # Set up logging
        today = date.today().strftime("%Y%m%d")
        now = datetime.now().strftime("%H:%M:%S")
        logfile = f'{self.log_dir}/event_dataset_{evt_name}_subs-{subs}_ses-{sessions}_{today}_{now}_log.txt'
        logger = setup_logging(self.log_dir, 'Export event dataset', logfile)
        logger.info('')
        
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
        
        # Set stage defaults
        if not stage:
            stage = ['NREM2', 'NREM3']
        
        # Set frequency bands
        frequency, flag = adap_bands_setup(self, adap_bands, None, None, 
                                           None, None, None, stage, concat_stage, 
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



    def split_events_freq(self, chan, evt_name, 
                                xml_dir = None, out_dir = None,
                                subs = 'all', sessions = 'all',
                                stage = None, concat_stage = False,
                                concat_cycle = True, cycle_idx = None,
                                grp_name = 'eeg', adap_bands = 'Fixed',
                                params = 'all', 
                                bins = 60, smooth_sigma = 1.5,
                                preferred_split = 13.0, 
                                prefer_window = (11.0, 15.0),
                                prefer_ratio = 0.6,
                                target_column="Peak power frequency (Hz)",
                                cluster_labels = ('low_peak', 'high_peak'),
                                export_only = False,
                                outfile = True):
        """
            Cluster events in exported event parameter CSVs file based on a
            specific parameter column, with a bias toward a preferred split 
            value within a window. Can be skipped with export_only=True.
            Also exports event datasets per cluster.
        """

        from seapipe.events.fish import FISH
        from seapipe.events.whales import cluster_peaks_preferred
        from seapipe.utils.load import select_input_dirs

        # Set up logging
        logger = setup_logging(self.log_dir, 'Cluster peak dataset (preferred) {event_name}', outfile)
        logger.info('')

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
        if cat[0] + cat[1] == 2:
            model = 'whole_night'
        elif cat[0] + cat[1] == 0:
            model = 'stage*cycle'
        elif cat[0] == 0:
            model = 'per_cycle'
        elif cat[1] == 0:
            model = 'per_stage'
        else:
            model = 'whole_night'

        # Format channels
        if isinstance(chan, str):
            chan = [chan]
            
        # Set stage defaults
        if not stage:
            stage = ['NREM2', 'NREM3']

        for evt in evts:
            # Determine input/output directories
            xml_dir_evt = select_input_dirs(self.outpath, xml_dir, evt)
            if not path.exists(xml_dir_evt):
                logger.warning(f"{xml_dir_evt} not found for {evt}; skipping.")
                continue

            if not out_dir:
                if not path.exists(self.outpath + '/datasets/'):
                    mkdir(self.outpath + '/datasets/')
                out_dir = self.outpath + '/datasets/spindle_types'
            if not path.exists(out_dir):
                mkdir(out_dir)
            final_out_dir = path.join(out_dir, f'spindle_{model}')
            if not path.exists(final_out_dir):
                mkdir(final_out_dir)
            
            if not export_only:
                logger.debug(f'Clustering peak frequencies in: {xml_dir_evt}')
                cluster_peaks_preferred(
                    xml_dir_evt, evt, chan,
                    bins=bins, smooth_sigma=smooth_sigma,
                    preferred_split=preferred_split,
                    prefer_window=prefer_window,
                    prefer_ratio=prefer_ratio,
                )

            cluster_root = Path(self.outpath) / "cluster_peaks" / evt
            cluster_root.mkdir(parents=True, exist_ok=True)
            for ch in chan:
                for cluster_label in cluster_labels:
                    logger.debug(f'Building cluster-specific files for {evt} ({cluster_label})')
                    cluster_dir = cluster_root / cluster_label
                    _build_cluster_param_tree(xml_dir_evt, cluster_dir, evt, ch, cluster_label, logger)
                    fish = FISH(self.rootpath, self.datapath, 
                                str(cluster_dir), out_dir,
                                chan, None, grp_name, 
                                stage, subs=subs, sessions=sessions)
                    fish.net(chan, evt, adap_bands, params, cat, cycle_idx, logger)
                    temp_dir = path.join(out_dir, f"{evt}_{model}")
                    if path.exists(temp_dir):
                        for fname in listdir(temp_dir):
                            if not fname.endswith(".csv"):
                                continue
                            new_name = fname.replace(f"{evt}_", f"{evt}_{cluster_label}_", 1)
                            rename(path.join(temp_dir, fname), path.join(final_out_dir, new_name))
                        shutil.rmtree(temp_dir, ignore_errors=True)

        return
    
    
    def pac_dataset(self, chan, evt_name = None, subs = 'all', sessions = 'all',
                          xml_dir = None, out_dir = None,  stage = None, 
                          concat_stage = False, concat_cycle = True, 
                          cycle_idx = None, grp_name = 'eeg', 
                          adap_bands_phase = 'Fixed', frequency_phase = (0.5, 1.25), 
                          adap_bands_amplitude = 'Fixed', frequency_amplitude = (11, 16),  
                          params = 'all', outfile=True):
        
        from seapipe.events.fish import FISH
        from seapipe.utils.load import select_input_dirs
        
        # Set up logging
        logger = setup_logging(self.log_dir, 'PAC dataset {event_name}', outfile)
        logger.info('')
        
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
        
        # Set stage defaults
        if not stage:
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
                                stage = None, freq_bands = ('SWA', 'Sigma'),
                                params = 'all', outfile = True):
        
        from seapipe.events.fish import FISH
        from seapipe.utils.load import select_input_dirs
        
        # Set up logging
        logger = setup_logging(self.log_dir, 'Cluster dataset {event_name}', outfile)
        logger.info('')
        
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
        
        # Set stage defaults
        if not stage:
            stage = ['NREM2','NREM3']
        
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
                                stage = None, 
                                concat_stage = False, concat_cycle = True, 
                                cycle_idx = None, grp_name = 'eeg', 
                                rater = None, params = 'all', 
                                general_opts = None, frequency_opts = None, 
                                filter_opts = None, epoch_opts = None, 
                                event_opts = None, outfile=True):
        
        from seapipe.spectrum.psa import (Spectrum, 
                                          default_epoch_opts, 
                                          default_event_opts,
                                          default_filter_opts, 
                                          default_frequency_opts, 
                                          default_general_opts)
        from seapipe.utils.load import select_input_dirs
        
        # Set up logging
        logger = setup_logging(self.log_dir, 'Power spectrum dataset', outfile)
        logger.info('')
        
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
        
        # Set stage defaults
        if not stage:
            stage = ['NREM1','NREM2','NREM3', 'REM']
        
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
                                   stage = None,
                                   cycle_idx = None, concat_cycle = True,
                                   concat_stage = True,
                                   general_opts = None, filter_opts = None,
                                   epoch_opts = None, event_opts = None,
                                   bandpower_opts = None, outfile = True):
        
        """Convenience wrapper that exports band-limited power vs. time."""

        # Local import avoids touching module-level imports per request
        from seapipe.spectrum.bandpower import (BandPowerTimecourse, 
                                                default_bandpower_opts)
        from seapipe.spectrum.psa import (default_epoch_opts, 
                                          default_event_opts, 
                                          default_filter_opts, 
                                          default_general_opts)
        from seapipe.utils.load import (check_chans,
                                        select_input_dirs)

        # Logging
        logger = setup_logging(self.log_dir, 'Bandpower timecourse', outfile)
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
        
        # Set stage defaults
        if not stage:
            stage = ['NREM1','NREM2','NREM3', 'REM']

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
    wonambi_row = None

    with open(csv_path, newline="") as f:
        reader = csv.reader(f)
        for row in reader:
            if not row:
                continue
            if header is None:
                if wonambi_row is None and row[0].strip().startswith("Wonambi"):
                    wonambi_row = row
                    continue
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

    return header, cluster_idx, events, count_meta, density_meta, wonambi_row


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


def _write_cluster_summary_csv(out_path, header, cluster_idx, cluster_label, summary, filtered, wonambi_row=None):
    cluster_column = "Peak cluster"
    if cluster_idx is None:
        header = header + [cluster_column]
        cluster_idx = len(header) - 1

    out_path.parent.mkdir(parents=True, exist_ok=True)
    writer_rows = []
    if wonambi_row is not None:
        writer_rows.append(wonambi_row)
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


def _build_cluster_param_tree(src_root, dest_root, evt_name, chan, cluster_label, logger):
    src_root = Path(src_root)
    dest_root = Path(dest_root)
    csv_paths = [
        p for p in src_root.rglob("*.csv")
        if p.is_file()
        and evt_name in str(p)
        and chan in str(p)
    ]
    for csv_path in csv_paths:
        if not csv_path.stem.endswith(evt_name):
            continue
        header, cluster_idx, events, count_meta, density_meta, wonambi_row = _load_cluster_events(csv_path)
        if header is None:
            logger.warning(f"Skipping {csv_path}: no header found.")
            continue
        if not events:
            logger.warning(f"Skipping {csv_path}: no event rows found.")
            continue

        try:
            summary, filtered = _summarize_cluster_events(header, cluster_idx, 
                                                          events, cluster_label, 
                                                          count_meta, density_meta)
        except Exception as exc:  # noqa: BLE001
            logger.warning(f"Skipping {csv_path}: {exc}")
            continue

        rel_parent = csv_path.relative_to(src_root).parent
        if csv_path.stem.endswith(evt_name):
            base = csv_path.stem[: -len(evt_name)]
        else:
            base = f"{csv_path.stem}_"
        new_stem = f"{base}{evt_name}"
        dest_file = dest_root / rel_parent / f"{new_stem}{csv_path.suffix}"
        _write_cluster_summary_csv(dest_file, header, cluster_idx, cluster_label, summary, filtered, wonambi_row)
        


def setup_logging(log_dir, logger_name, outfile):
    
    from seapipe.utils.logs import create_logger, create_logger_outfile
    
    # Set up logging
    if outfile == True:
        today = date.today().strftime("%Y%m%d")
        now = datetime.now().strftime("%H%M%S")
        logfile_prefix = '_'.join(logger_name.split(' '))
        logfile = f'{log_dir}/{logfile_prefix}_{today}_{now}_log.txt'
        logger = create_logger_outfile(logfile=logfile, name=f'{logger_name}')
        logger.info('')
        logger.info(f"-------------- New call of '{logger_name}' evoked at {now} --------------")
    elif outfile:
        logfile = f'{log_dir}/{outfile}'
        logger = create_logger_outfile(logfile=logfile, name=f'{logger_name}')
    else:
        logger = create_logger(f'{logger_name}')
        
    return logger
