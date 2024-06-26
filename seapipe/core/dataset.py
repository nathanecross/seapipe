#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 12:07:36 2023

@author: nathancross
"""
from os import listdir, mkdir, path, remove, walk
from pandas import DataFrame, read_csv
from seapipe.events.fish import FISH
from seapipe.events.whales import whales
from seapipe.events.seasnakes import seasnakes
from seapipe.spectrum.psa import (Spectrum, default_epoch_opts, default_event_opts,
                     default_fooof_opts, default_filter_opts, default_frequency_opts, 
                     default_general_opts)
from seapipe.spectrum.spectrogram import event_spectrogram, event_spectrogram_grouplevel
from seapipe.stats import sleepstats
from seapipe.utils.audit import (check_dataset, check_fooof, extract_channels, make_bids,
                        track_processing)
from seapipe.utils.logs import create_logger, create_logger_outfile
from seapipe.utils.load import (check_chans, check_adap_bands, select_input_dirs,
                        select_ouput_dirs)


## TO DO:
#   - add selection of subs to be readable from 'tracking.tsv'
#   - add logging to save to output file (not implemented for all functions)
#   - update adapted bands in tracking.tsv
#   - create catch for errors in tracking sheet around ',' (chans, adap_bands etc.)
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

    Attributes
    ----------
    rootpath : str
        name of the root level directory
    datapath : str
        name of the directory containing the raw data (recordings and annotations)
    outpath : str
        name of the directory containing the output (analysed) data

    """
        
    def __init__(self, indir, outfile=False, filetype='.edf'):
        
        self.rootpath = indir
        self.datapath = indir + '/DATA'
        self.outpath = indir + '/OUT'
        if not path.exists(self.outpath):
            mkdir(self.outpath)
        self.outfile = outfile
        if not path.exists(f'{self.outpath}/audit'):
            mkdir(f'{self.outpath}/audit')
        self.audit_init = check_dataset(self.rootpath, self.outfile, filetype)
        
        self.tracking = {}
        self.track(subs='all', ses='all', 
                   step=['staging','spindle','slowwave','pac','sync','psa'],
                   show=False, log=False)
    
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
            self.audit_update = check_dataset(self.rootpath, False, filetype, 
                                              tracking, logger)
        else:
            if not outfile:
                outfile = self.outfile
            out = f'{out_dir}/{outfile}'
            if path.exists(out):
                remove(out)
            logger = create_logger_outfile(outfile, name = 'Audit')
            logger.propagate = False
            self.audit_update = check_dataset(self.rootpath, out, filetype, 
                                              tracking, logger)
            
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
            subs = [x for x in listdir(self.datapath) if '.' not in x]
        else:
            subs = read_csv(f'{self.rootpath}/{subs}', sep='\t')
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

    def make_bids(self, origin = 'SCN'):
        make_bids(self.datapath, origin=origin)
        
    def extract_channels(self, exclude = None):
        extract_channels(self.datapath, exclude=exclude)
        
    
    
    #--------------------------------------------------------------------------
    '''
    SLEEP EVENTS DETECTIONS
    
    detect_spectral_peaks,
    detect_slow_oscillations,
    detect_spindles,
    
    
    '''
    
    def detect_spectral_peaks(self, xml_dir = None, out_dir = None, 
                              subs = 'all', sessions = 'all', chan = None, 
                              ref_chan = None, grp_name = 'eeg', rater = None, 
                              frequency = (9,16), stage = None, 
                              concat_cycle = True, concat_stage = False,
                              general_opts = None, frequency_opts = None,
                              filter_opts = None, epoch_opts = None, 
                              event_opts = None, fooof_opts = None, 
                              filetype = '.edf', suffix = None):
        
        # Set up logging
        logger = create_logger('Detect spectral peaks')
        
        # Set input/output directories
        in_dir = self.datapath
        log_dir = self.outpath + '/audit/logs/'
        if not path.exists(log_dir):
            mkdir(log_dir)
        if not xml_dir:
            xml_dir = f'{self.outpath}/staging'   
        if not out_dir:
            out_dir = f'{self.outpath}/fooof' 
        if not path.exists(out_dir):
            mkdir(out_dir)
            
        # Set suffix for output filename
        if not suffix:
            suffix = f'{frequency[0]}-{frequency[1]}Hz'
            
        # Set channels
        chan, ref_chan = check_chans(self.rootpath, chan, ref_chan, logger)
        
        # Format concatenation
        cat = (int(concat_cycle),int(concat_stage),1,1)
        
        # Default stage
        if stage == None:
            stage = ['NREM2','NREM3']
        
        # Check annotations directory exists, run detection
        if not path.exists(xml_dir):
            logger.info('')
            logger.critical(f"{xml_dir} doesn't exist. Sleep staging has not been run or hasn't been converted correctly.")
            logger.info('Check documentation for how to set up staging data:')
            logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
            logger.info('-' * 10)
        else:   
            spectrum = Spectrum(in_dir, xml_dir, out_dir, log_dir, chan, ref_chan, 
                             grp_name, stage, frequency, rater, subs, sessions)
            
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
            
            spectrum.fooof_it(general_opts, frequency_opts, filter_opts, 
                                      epoch_opts, event_opts, fooof_opts, 
                                      rater=None, grp_name = 'eeg', 
                                      cat = cat, cycle_idx=None,
                                      filetype = filetype, suffix = suffix)  
        return 
    
    
    def detect_slow_oscillations(self, xml_dir=None, out_dir=None, subs='all', 
                        sessions='all', filetype='.edf', method = None, 
                        chan=None, ref_chan=None, rater=None, grp_name='eeg',
                        stage = None, cycle_idx=None, 
                        duration=(0.2, 2), average_channels = False, 
                        invert = None, outfile=True):
        
        # Set up logging
        logger = create_logger('Detect slow oscillations')
        logger.info('')
        logger.debug("Commencing SO detection pipeline.")
        logger.info('')
        
        # Set input/output directories
        in_dir = self.datapath
        log_dir = self.outpath + '/audit/logs/'
        if not path.exists(log_dir):
            mkdir(log_dir)
        if not xml_dir:
            xml_dir = f'{self.outpath}/staging'   
        if not out_dir:
            out_dir = f'{self.outpath}/slowwave'    
        if not path.exists(out_dir):
            mkdir(out_dir)
        
        # Set channels
        chan, ref_chan = check_chans(self.rootpath, chan, ref_chan, logger)
        
        # Check inversion
        if invert == None:
            invert = check_chans(self.rootpath, None, False, logger)
        elif type(invert) != bool:
            logger.critical(f"The argument 'invert' must be set to either: 'True', 'False' or 'None'; but it was set as {invert}.")
            logger.info('Check documentation for how to set up staging data:')
            logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
            logger.info('-' * 10)
            return
            
        # Format concatenation
        cat = (1,1,1,1)
        
        # Default stage
        if stage == None:
            stage = ['NREM2','NREM3']
            
        # Default method
        if method == None:
            method = ['Staresina2015']
        
        # Check annotations directory exists, run detection
        if not path.exists(xml_dir):
            logger.info('')
            logger.critical(f"{xml_dir} doesn't exist. Sleep staging has not been run or hasn't been converted correctly.")
            logger.info('Check documentation for how to set up staging data:')
            logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
            logger.info('-' * 10)
        else:   
           SO = seasnakes(in_dir, xml_dir, out_dir, log_dir, chan, ref_chan, 
                            grp_name, stage, rater, subs, sessions,
                            self.tracking) 
           SO.detect_slowosc(method, cat, cycle_idx, duration, 
                                  average_channels, invert, filetype, outfile)
           try:
               self.tracking = self.tracking | SO.tracking
           except:
               self.tracking = {**self.tracking, **SO.tracking}
        return
    
    
    def detect_spindles(self, xml_dir = None, out_dir = None, subs = 'all', 
                        sessions = 'all', filetype = '.edf', method = None, 
                        chan = None, ref_chan = None, rater = None, stage = None, 
                        grp_name = 'eeg', cycle_idx = None, frequency = (11,16), 
                        adap_bands = 'Fixed', adap_bw = 4, peaks = None, 
                        reject_artf = ['Artefact', 'Arou', 'Arousal'], 
                        duration =( 0.5, 3), outfile = True):
        
        # Set up logging
        logger = create_logger('Detect spindles')
        logger.info('')
        logger.debug("Commencing spindle detection pipeline.")
        logger.info('')
        
        # Set input/output directories
        in_dir = self.datapath
        log_dir = self.outpath + '/audit/logs/'
        if not path.exists(log_dir):
            mkdir(log_dir)
        if not xml_dir:
            xml_dir = f'{self.outpath}/staging'   
        if not out_dir:
            out_dir = f'{self.outpath}/spindle'    
        if not path.exists(out_dir):
            mkdir(out_dir)
        
        # Set channels
        chan, ref_chan = check_chans(self.rootpath, chan, ref_chan, logger)
        if not isinstance(chan, DataFrame) and not isinstance(chan, list):
            return
        elif isinstance(ref_chan, str):
            return
        
        # Format concatenation
        cat = (1,0,1,1)
        
        # Default stage
        if stage == None:
            stage = ['NREM2','NREM3']
        
        # Default method
        if method == None:
            method = ['Moelle2011']
            
        # Check for adapted bands
        if adap_bands == 'Fixed':
            logger.debug('Detection using FIXED frequency bands has been selected (adap_bands = Fixed)')
        elif adap_bands == 'Manual':
            logger.debug('Detection using ADAPTED (user-provided) frequency bands has been selected (adap_bands = Manual)')
            logger.debug(f"Checking for spectral peaks in {self.rootpath}/'tracking.tsv' ")
            flag = check_adap_bands(self, subs, sessions, chan, logger)
            if flag == 'error':
                return
            elif flag == 'review':
                logger.info('')
                logger.warning(f"Some spectral peak entries in 'tracking.tsv' are inconsistent or missing. In these cases, detection will revert to fixed bands: {frequency[0]}-{frequency[1]}Hz")
                logger.info('')
            peaks = check_chans(self.rootpath, None, False, logger)
        elif adap_bands == 'Auto':   
            frequency = (9,16)             
            logger.debug('Detection using ADAPTED (automatic) frequency bands has been selected (adap_bands = Auto)')
            self.track(step='fooof', show = False, log = False)
            if not type(chan) == type(DataFrame()):
                logger.critical("For adap_bands = Auto, the argument 'chan' must be 'None' and specfied in 'tracking.csv'")
                return
            else:
                flag, pk_chan, pk_sub, pk_ses = check_fooof(self, frequency, 
                                                            chan, ref_chan, 
                                                            stage, cat, 
                                                            cycle_idx, logger)
            if flag == 'error':
                logger.critical('Error in reading channel names, check tracking sheet.')
                logger.info("Check documentation for how to set up channel names in tracking.tsv':")
                logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
                logger.info('-' * 10)
                return
            elif flag == 'review':
                logger.debug('Spectral peaks have not been found for all subs, analysing the spectral parameters prior to spindle detection..')
                if 'NREM2' in stage and 'NREM3' in stage:
                    concat_cycle = True
                    concat_stage = True
                else:
                    concat_cycle = True 
                    concat_stage = False
                self.detect_spectral_peaks(subs = pk_sub, 
                                           sessions = pk_ses, 
                                           chan = pk_chan, 
                                           concat_cycle=concat_cycle, 
                                           concat_stage=concat_stage)
  
        # Check annotations directory exists, run detection
        if not path.exists(xml_dir):
            logger.info('')
            logger.critical(f"{xml_dir} doesn't exist. Sleep staging has not been run or hasn't been converted correctly.")
            logger.info('Check documentation for how to set up staging data:')
            logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
            logger.info('-' * 10)
        else:   
           spindle = whales(in_dir, xml_dir, out_dir, log_dir, chan, ref_chan, 
                            grp_name, stage, frequency, rater, subs, sessions, 
                            reject_artf, self.tracking) 
           spindle.whale_it(method, cat, cycle_idx, adap_bands, peaks, adap_bw, 
                            duration, filetype, outfile)
           try:
               self.tracking = self.tracking | spindle.tracking
           except:
               self.tracking = {**self.tracking, **spindle.tracking}
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
        
        # Set up logging
        logger = create_logger('Event spectrogram')
        logger.info('')
        logger.debug("Creating spectrogram of events.")
        logger.info('')
        
        # Set input/output directories
        in_dir = self.datapath
        log_dir = self.outpath + '/audit/logs/'
        if not path.exists(log_dir):
            mkdir(log_dir)
        if not xml_dir:
            xml_dir = f'{self.outpath}/staging'   
        if not out_dir:
            out_dir = f'{self.outpath}/spindle'    
        if not path.exists(out_dir):
            mkdir(out_dir)
        
        # Format concatenation
        cat = (int(concat_cycle),int(concat_stage),1,1)
        
        # Check inversion
        if invert == None:
            invert = check_chans(self.rootpath, chan, False, logger)
        elif type(invert) != bool:
            logger.critical(f"The argument 'invert' must be set to either: 'True', 'False' or 'None'; but it was set as {invert}.")
            logger.info('Check documentation for how to set up staging data:')
            logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
            logger.info('-' * 10)
            return
        
        if not filter_opts:
            filter_opts = default_filter_opts()
        
        
        if not evt_type:
            logger.warning('No event type (evt_type) has been specified. Spectrogram will be run on epochs instead. This may take some time...')
        
        event_spectrogram(self, in_dir, xml_dir, out_dir, subs, sessions, stage, 
                              cycle_idx, chan, ref_chan, rater, grp_name, 
                              evt_type, buffer, invert, cat, filter_opts,  
                              outfile, filetype, progress, self.tracking)
        
        
        
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
                                 times = None, rater = None, outfile = True):
        
        # Set up logging
        logger = create_logger('Export macro stats')
        
        # Set input/output directories
        log_dir = self.outpath + '/audit/logs/'
        if not path.exists(log_dir):
            mkdir(log_dir)
        xml_dir = select_input_dirs(self, xml_dir, 'macro')
        out_dir = select_ouput_dirs(self, out_dir, 'macro')
        
        # Set channels
        times, ref_chan = check_chans(self.rootpath, None, True, logger)
        
        self.track(subs='all', ses='all', 
                   step=['staging'],
                   show=False, log=True)
        
        sleepstats.export_sleepstats(xml_dir, out_dir, subs, sessions, 
                                     rater, times, log_dir, outfile)
        return
    
    def macro_dataset(self, xml_dir = None, out_dir = None, rater = None, 
                      subs = 'all', sessions = 'all', outfile = True):
         
         # Set input/output directories
         log_dir = self.outpath + '/audit/logs/'
         if not path.exists(log_dir):
             mkdir(log_dir)
         if not path.exists(self.outpath + '/datasets/'):
             mkdir(self.outpath + '/datasets/')
         out_dir = self.outpath + '/datasets/macro/'
         
         xml_dir = select_input_dirs(self, xml_dir, 'macro')
         out_dir = select_ouput_dirs(self, out_dir, 'macro')
         
         sleepstats.sleepstats_from_csvs(xml_dir, out_dir, rater,   
                                 subs, sessions, log_dir, outfile)
         return
    
    def export_eventparams(self, xml_dir = None, out_dir = None, subs = 'all', 
                           sessions = 'all', chan = None, ref_chan = None, 
                           rater=None, stage = None, grp_name = 'eeg', 
                           concat_cycle = True, concat_stage = False, 
                           cycle_idx = None, keyword = None, evt_name = 'spindle', 
                           frequency = (11,16), segs = None, adap_bands = False, 
                           param_keys = 'all',  
                           epoch_dur = 30, n_fft_sec = 4, Ngo = {'run':False},
                           outfile = True):
        
        # Set up logging
        logger = create_logger('Export params')
        
        # Set input/output directories
        in_dir = self.datapath
        log_dir = self.outpath + '/audit/logs/'
        if not path.exists(log_dir):
            mkdir(log_dir)
        xml_dir = select_input_dirs(self, xml_dir, evt_name)
        
        # Check annotations directory exists
        if not path.exists(xml_dir):
            logger.info('')
            logger.critical(f"{xml_dir} doesn't exist. Event detection has not been run or an incorrect event type has been selected.")
            logger.info('Check documentation for how to run a pipeline:')
            logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
            logger.info('-' * 10)
            return
        
        out_dir = select_ouput_dirs(self, out_dir, evt_name)
       
        if adap_bands:
            xml_dir = f'{xml_dir}_adapted'
            out_dir = f'{out_dir}_adapted'
        
        # Set channels
        chan, ref_chan = check_chans(self.rootpath, chan, ref_chan, logger)
        
        # Format concatenation
        cat = (int(concat_cycle),int(concat_stage),1,1)
        
        # Default stage
        if stage == None:
            stage = ['NREM2','NREM3']
            
        # Run line
        fish = FISH(in_dir, xml_dir, out_dir, log_dir, chan, ref_chan, grp_name, 
                          stage, frequency, rater, subs, sessions) 
        fish.line(keyword, evt_name, cat, segs, cycle_idx, adap_bands, 
                            param_keys, epoch_dur, 
                            n_fft_sec, Ngo, outfile)
        return
    
    
    def event_dataset(self, xml_dir = None, out_dir = None, subs = 'all', 
                            sessions = 'all', chan = None, ref_chan = None, 
                            rater = None, stage = None, concat_stage = False, 
                            concat_cycle = True, cycle_idx = None, grp_name = 'eeg', 
                            frequency = (11,16), adap_bands = False, 
                            evt_name = 'spindle', params = 'all', outfile=True):
        
        # Set up logging
        logger = create_logger('Create dataset')
        
        if type(evt_name) is not str:
            logger.critical("'evt_name' MUST be a string only (i.e.) only 1 event at a time.")
            return
        
        # Set input/output directories
        in_dir = self.datapath
        log_dir = self.outpath + '/audit/logs/'
        if not path.exists(log_dir):
            mkdir(log_dir)
        if not path.exists(self.outpath + '/datasets/'):
            mkdir(self.outpath + '/datasets/')
        out_dir = self.outpath + f'/datasets/{evt_name}'
        xml_dir = select_input_dirs(self, xml_dir, evt_name)
        
        if not path.exists(xml_dir):
            logger.info('')
            logger.critical(f"{xml_dir} doesn't exist. Event detection has not been run or an incorrect event type has been selected.")
            logger.info('Check documentation for how to run a pipeline:')
            logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
            logger.info('-' * 10)
            return
        
        out_dir = select_ouput_dirs(self, out_dir, evt_name)
        
        if adap_bands:
            xml_dir = f'{xml_dir}_adapted'
            out_dir = f'{out_dir}_adapted'
            
        # Set channels
        chan, ref_chan = check_chans(self.rootpath, chan, ref_chan, logger)
        
        # Format concatenation
        cat = (int(concat_cycle),int(concat_stage),1,1)
        
        # Default stage
        if stage == None:
            stage = ['NREM2','NREM3']
            
        fish = FISH(in_dir, xml_dir, out_dir, log_dir, chan, ref_chan, grp_name, 
                          stage, frequency, rater, subs, sessions) 
        fish.net(chan, evt_name, params, cat, cycle_idx, outfile)
        
        return
    

    
    
    
    
        
        