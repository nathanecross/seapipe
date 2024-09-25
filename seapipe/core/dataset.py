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
from seapipe.events.seabass import seabass
from seapipe.pac.octopus import octopus, pac_method
from seapipe.spectrum.psa import (Spectrum, default_epoch_opts, default_event_opts,
                     default_fooof_opts, default_filter_opts, default_frequency_opts, 
                     default_general_opts,default_norm_opts)
from seapipe.spectrum.spectrogram import event_spectrogram, event_spectrogram_grouplevel
from seapipe.stats import sleepstats
from seapipe.utils.audit import (check_dataset, check_fooof, extract_channels, make_bids,
                        track_processing)
from seapipe.utils.logs import create_logger, create_logger_outfile
from seapipe.utils.load import (check_chans, check_adap_bands, read_tracking_sheet, 
                                select_input_dirs, select_output_dirs,)


## TO DO:
#   - add frequency option in export_eventparams to be default for SO & spindle
#   - add in log for detection whether auto, fixed or adapted bands was run
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
        elif not isinstance(subs, list):
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
        
    
    
    def power_spectrum(self, xml_dir = None, out_dir = None, 
                             subs = 'all', sessions = 'all', chan = None, 
                             ref_chan = None, grp_name = 'eeg', rater = None, 
                             stage = ['NREM1','NREM2','NREM3', 'REM'], 
                             cycle_idx = None, concat_cycle = True, 
                             concat_stage = False, general_opts = None, 
                             frequency_opts = None, filter_opts = None, 
                             epoch_opts = None, event_opts = None, 
                             norm = None, norm_opts = None, 
                             filetype = '.edf'):
        
        # Set up logging
        logger = create_logger('Power spectrum')
        
        # Set input/output directories
        in_dir = self.datapath
        log_dir = self.outpath + '/audit/logs/'
        if not path.exists(log_dir):
            mkdir(log_dir)
        if not xml_dir:
            xml_dir = f'{self.outpath}/staging'   
        if not out_dir:
            out_dir = f'{self.outpath}/spectrum' 
        if not path.exists(out_dir):
            mkdir(out_dir)
            
        # Set channels
        chan, ref_chan = check_chans(self.rootpath, chan, ref_chan, logger)
        
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
               int(epoch_opts['concat_signal']),
               int(event_opts['concat_events']),
               )
        
        # Set suffix for output filename
        if not general_opts['suffix']:
            general_opts['suffix'] = f"{frequency_opts['frequency'][0]}-{frequency_opts['frequency'][1]}Hz"
        
        # Check annotations directory exists, run detection
        if not path.exists(xml_dir):
            logger.info('')
            logger.critical(f"{xml_dir} doesn't exist. Sleep staging has not been run or hasn't been converted correctly.")
            logger.info('Check documentation for how to set up staging data:')
            logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
            logger.info('-' * 10)
        else:   
            spectrum = Spectrum(in_dir, xml_dir, out_dir, log_dir, chan, ref_chan, 
                                grp_name, stage, cat, rater, cycle_idx, subs, 
                                sessions, self.tracking)
            
                                  
            spectrum.powerspec_it(general_opts, frequency_opts, filter_opts, 
                                  epoch_opts, event_opts, norm, norm_opts, 
                                  filetype) 
            
            try:
                self.tracking = self.tracking | spectrum.tracking
            except:
                self.tracking = {**self.tracking, **spectrum.tracking}
            
        return 
        

    #--------------------------------------------------------------------------
    '''
    SLEEP EVENTS DETECTIONS
    
    sleep_staging
    detect_spectral_peaks,
    detect_slow_oscillations,
    detect_spindles,
    
    
    '''
    def detect_sleep_stages(self, eeg_chan = None, xml_dir = None, out_dir = None, 
                                  subs = 'all', sessions = 'all', 
                                  method = 'Vallat2021', qual_thresh = 0.5,
                                  filetype = '.edf', ref_chan = None, 
                                  eog_chan = None, emg_chan = None, 
                                  rater = None, invert = False, outfile = True):
        
        # Set up logging
        logger = create_logger('Detect sleep stages')
        logger.info('')
        logger.debug("Commencing sleep stage detection pipeline.")
        logger.info('')
        
        # Set input/output directories
        in_dir = self.datapath
        log_dir = self.outpath + '/audit/logs/'
        if not path.exists(log_dir):
            mkdir(log_dir)
        if not xml_dir:
            xml_dir = f'{self.outpath}/staging'   
        if not out_dir:
            out_dir = f'{self.outpath}/staging'    
        if not path.exists(out_dir):
            mkdir(out_dir)
        
        # Check subs
        if not subs:
            tracking = read_tracking_sheet(self.rootpath, logger)
            subs = [x for x in list(set(tracking['sub']))]
            subs.sort()
        if not sessions:
            sessions = read_tracking_sheet(self.rootpath, logger)
        
        # Set channels
        eeg_chan, ref_chan = check_chans(self.rootpath, eeg_chan, ref_chan, logger)
        
        # Check inversion
        if invert == None:
            invert = check_chans(self.rootpath, None, False, logger)
        elif type(invert) != bool:
            logger.critical(f"The argument 'invert' must be set to either: 'True', 'False' or 'None'; but it was set as {invert}.")
            logger.info('')
            logger.info('Check documentation for how to set up staging data: https://seapipe.readthedocs.io/en/latest/index.html')
            logger.info('-' * 10)
            logger.critical('Sleep stage detection finished with ERRORS. See log for details.')
            return
    
        # Check annotations directory exists, run detection
        if not path.exists(xml_dir):
            logger.info('')
            logger.critical(f"{xml_dir} doesn't exist. Sleep staging has not been run or hasn't been converted correctly.")
            logger.info('')
            logger.info('Check documentation for how to set up staging data: https://seapipe.readthedocs.io/en/latest/index.html')
            logger.info('-' * 10)
            logger.critical('Sleep stage detection finished with ERRORS. See log for details.')
        else:   
           stages = seabass(in_dir, xml_dir, out_dir, log_dir, eeg_chan, 
                            ref_chan, eog_chan, emg_chan, rater, subs, sessions, 
                            self.tracking) 
           stages.detect_stages(method, qual_thresh, invert, filetype, 
                                outfile)
           try:
               self.tracking = self.tracking | stages.tracking
           except:
               self.tracking = {**self.tracking, **stages.tracking}
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
        
        # Check annotations directory exists, run detection
        if not path.exists(xml_dir):
            logger.info('')
            logger.critical(f"{xml_dir} doesn't exist. Sleep staging has not been run or hasn't been converted correctly.")
            logger.info('Check documentation for how to set up staging data:')
            logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
            logger.info('-' * 10)
        else:   
            spectrum = Spectrum(in_dir, xml_dir, out_dir, log_dir, chan, 
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
                              filetype = filetype)  
                
        return 
    
    
    def detect_slow_oscillations(self, xml_dir=None, out_dir=None, subs='all', 
                                       sessions='all', filetype='.edf', 
                                       method = ['Staresina2015'], chan=None,
                                       ref_chan=None, rater=None, grp_name='eeg', 
                                       stage = ['NREM2','NREM3'], cycle_idx=None, 
                                       duration=(0.2, 2), invert = None,
                                       reject_artf = ['Artefact', 'Arou', 'Arousal'],
                                       average_channels = False, outfile=True):
        
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
            logger.critical(f"The argument 'invert' must be set to either: 'True', 'False' or 'None'; but it was set as {invert}.")
            logger.info('Check documentation for how to set up staging data:')
            logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
            logger.info('-' * 10)
            logger.critical('SO detection finished with ERRORS. See log for details.')
            return
            
        # Format concatenation
        cat = (1,1,1,1)
    
        # Check annotations directory exists, run detection
        if not path.exists(xml_dir):
            logger.info('')
            logger.critical(f"{xml_dir} doesn't exist. Sleep staging has not been run or hasn't been converted correctly.")
            logger.info('Check documentation for how to set up staging data:')
            logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
            logger.info('-' * 10)
            logger.critical('SO detection finished with ERRORS. See log for details.')
        else:   
           SO = seasnakes(in_dir, xml_dir, out_dir, log_dir, chan, ref_chan, 
                            grp_name, stage, rater, subs, sessions, 
                            self.tracking, reject_artf) 
           SO.detect_slowosc(method, cat, cycle_idx, duration, 
                                  average_channels, invert, filetype, outfile)
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
                              adap_bw = 4, duration =( 0.5, 3),
                              reject_artf = ['Artefact', 'Arou', 'Arousal'], 
                              outfile = True):
        
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
            for met in method:
                out_dir = select_output_dirs(self.outpath, out_dir, met)  
        if not path.exists(out_dir):
            mkdir(out_dir)
        
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
        
        # Format concatenation
        if concat_cycle == True:
            cat = (1,0,1,1)
        else:
            cat = (0,0,1,1)
        
        # Check for adapted bands
        if adap_bands == 'Fixed':
            logger.debug('Detection using FIXED frequency bands has been selected (adap_bands = Fixed)')
            if not frequency:
                frequency = (11,16)
        elif adap_bands == 'Manual':
            logger.debug('Detection using ADAPTED (user-provided) frequency bands has been selected (adap_bands = Manual)')
            logger.debug(f"Checking for spectral peaks in {self.rootpath}/'tracking.tsv' ")
            flag = check_adap_bands(self.rootpath, subs, sessions, chan, logger)
            if flag == 'error':
                logger.critical('Spindle detection finished with ERRORS. See log for details.')
                return
            elif flag == 'review':
                logger.info('')
                logger.warning(f"Some spectral peak entries in 'tracking.tsv' are inconsistent or missing. In these cases, detection will revert to fixed bands: {frequency[0]}-{frequency[1]}Hz")
                logger.info('')
        elif adap_bands == 'Auto': 
            if not frequency:
                frequency = (9,16)           
            logger.debug('Detection using ADAPTED (automatic) frequency bands has been selected (adap_bands = Auto)')
            self.track(subs, sessions, step = 'fooof', show = False, log = False)
            if not type(chan) == type(DataFrame()):
                logger.critical("For adap_bands = Auto, the argument 'chan' must be 'None' and specfied in 'tracking.csv'")
                logger.critical('Spindle detection finished with ERRORS. See log for details.')
                return
            else:
                flag, pk_chan, pk_sub, pk_ses = check_fooof(self, frequency, 
                                                                  chan, ref_chan, 
                                                                  stage, 
                                                                  cat,
                                                                  cycle_idx, 
                                                                  logger)
            if flag == 'error':
                logger.critical('Error in reading channel names, check tracking sheet.')
                logger.info("Check documentation for how to set up channel names in tracking.tsv':")
                logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
                logger.info('-' * 10)
                logger.critical('Spindle detection finished with ERRORS. See log for details.')
                return
            elif flag == 'review':
                logger.debug('Spectral peaks have not been found for all subs, analysing the spectral parameters prior to spindle detection..')
                for (sub,ses) in zip(pk_sub,pk_ses):
                    self.detect_spectral_peaks(subs = [sub], 
                                           sessions = [ses], 
                                           chan = pk_chan, 
                                           frequency = frequency,
                                           stage = stage, cycle_idx = cycle_idx,
                                           concat_cycle=concat_cycle, 
                                           concat_stage=True)
        # Check annotations directory exists, run detection
        self.track(step='fooof', show = False, log = False)
        if not path.exists(xml_dir):
            logger.info('')
            logger.critical(f"{xml_dir} doesn't exist. Sleep staging has not been run or hasn't been converted correctly.")
            logger.info('Check documentation for how to set up staging data:')
            logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
            logger.info('-' * 10)
            logger.critical('Spindle detection finished with ERRORS. See log for details.')
        else:   
           spindle = whales(self.rootpath, in_dir, xml_dir, out_dir, log_dir, 
                            chan, ref_chan, grp_name, stage, frequency, rater, 
                            subs, sessions, reject_artf, self.tracking) 
           spindle.whale_it(method, cat, cycle_idx, adap_bands, adap_bw, 
                            duration, filetype, outfile)
           try:
               self.tracking = self.tracking | spindle.tracking
           except:
               self.tracking = {**self.tracking, **spindle.tracking}
        return
    
    
    def whales(self, xml_dir = None, out_dir = None, subs = 'all', 
                     sessions = 'all', filetype = '.edf', 
                     method = ['Moelle2011', 'Ray2015'], 
                     merge_type = 'consensus', chan = None, 
                     ref_chan = None, rater = None, 
                     stage = ['NREM2','NREM3'], grp_name = 'eeg', 
                     cycle_idx = None,
                     s_freq = None, keyword = None, min_duration = 0.3,
                     frequency = None, adap_bands = 'Fixed', 
                     adap_bw = 4, peaks = None, duration =( 0.5, 3),
                     reject_artf = ['Artefact', 'Arou', 'Arousal'], 
                     outfile = True):
        
        # Set up logging
        logger = create_logger('Detect spindles (WHALES)')
        logger.info('')
        logger.debug("Commencing spindle optimisation pipeline.")
        logger.info('')
        
        # Set input/output directories
        in_dir = self.datapath
        log_dir = self.outpath + '/audit/logs/'
        if not path.exists(log_dir):
            mkdir(log_dir)
            
        xml_dir = select_input_dirs(self.outpath, xml_dir, 'spindle')  
        if not out_dir:
            out_dir = xml_dir  
        if not path.exists(out_dir):
            mkdir(out_dir)
        
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
            return
        elif isinstance(ref_chan, str):
            return
        
        spindle = whales(self.rootpath, in_dir, xml_dir, out_dir, log_dir, chan, ref_chan, 
                         grp_name, stage, frequency, rater, subs, sessions, 
                         reject_artf, self.tracking) 
        spindle.whales(method, merge_type, chan, rater, stage, ref_chan, grp_name, 
                       keyword, cs_thresh, min_duration, s_freq, 
                       frequency = (11, 16), duration= (0.5, 3), 
                       evt_out = 'spindle', weights = None,
                       outfile = 'export_params_log.txt', filetype = '.edf')
    
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
    PHASE AMPLITUDE COUPLING.
    
    mean_amps -> runs Phase Amplitude Coupling analyses on sleep EEG data. 
    
    
    '''    
    def pac(self, xml_dir = None, out_dir = None, subs = 'all', sessions = 'all', 
                  chan = None, ref_chan = None, rater = None, grp_name = 'eeg', 
                  stage = ['NREM2','NREM3'], concat_stage = True, 
                  cycle_idx = None, concat_cycle = True, filetype = '.edf', 
                  method = 'MI', surrogate = 'Time lag', correction = 'Z-score',
                  adap_bands_phase = 'Fixed', frequency_phase = (0.5, 1.25), 
                  adap_bands_amplitude = 'Fixed', frequency_amplitude = (11, 16),
                  evt_name = 'Staresina2015', min_dur = 1, nbins = 18, invert = None,
                  adap_bw = 4, peaks = None, frequency_opts = None, 
                  filter_opts = None, epoch_opts = None, event_opts = None, 
                  reject_artf = ['Artefact', 'Arou', 'Arousal'], 
                  progress = True, outfile = True):
        
        # Set up logging
        logger = create_logger('Phase-amplitude coupling')
        logger.info('')
        
        # Set input/output directories
        in_dir = self.datapath
        log_dir = self.outpath + '/audit/logs/'
        if not path.exists(log_dir):
            mkdir(log_dir) 
        if not xml_dir:
            xml_dir = select_input_dirs(self.outpath, xml_dir, evt_name) 
        if not out_dir:
            out_dir = select_output_dirs(self.outpath, out_dir, 'pac')  
        if not path.exists(out_dir):
            mkdir(out_dir)
        
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
        
        # Format concatenation
        cat = (int(concat_cycle),int(concat_stage),0,0)
        
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
            logger.critical(f"The argument 'invert' must be set to either: 'True', 'False' or 'None'; but it was set as {invert}.")
            logger.info('Check documentation for how to set up staging data:')
            logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
            logger.info('-' * 10)
            logger.critical('Phase amplitude coupling finished with ERRORS. See log for details.')
            return

        Octopus = octopus(self.rootpath, in_dir, xml_dir, out_dir, log_dir, chan, ref_chan, 
                          grp_name, stage, rater, subs, sessions, reject_artf,
                          self.tracking)
        
        Octopus.pac_it(cycle_idx, cat, nbins, filter_opts, epoch_opts, 
                       frequency_opts, event_opts, filetype, idpac, evt_name, 
                       min_dur, adap_bands_phase, frequency_phase, 
                       adap_bands_amplitude, frequency_amplitude, 
                       peaks, adap_bw, invert, progress, outfile)
    
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
        xml_dir = select_input_dirs(self.outpath, xml_dir, 'macro')
        out_dir = select_output_dirs(self.outpath, out_dir, 'macro')
        
        # Set channels
        times, ref_chan = check_chans(self.rootpath, None, True, logger)
        
        self.track(subs='all', ses='all', 
                   step=['staging'],
                   show=False, log=True)
        
        sleepstats.export_sleepstats(xml_dir, out_dir, subs, sessions, 
                                     rater, times, log_dir, outfile)
        return
    
    def macro_dataset(self, xml_dir = None, out_dir = None, 
                      subs = 'all', sessions = 'all', outfile = True):
         
         # Set input/output directories
         log_dir = self.outpath + '/audit/logs/'
         if not path.exists(log_dir):
             mkdir(log_dir)
         if not path.exists(self.outpath + '/datasets/'):
             mkdir(self.outpath + '/datasets/')
         out_dir = self.outpath + '/datasets/macro/'
         
         xml_dir = select_input_dirs(self.outpath, xml_dir, 'macro')
         out_dir = select_output_dirs(self.outpath, out_dir, 'macro')
         
         sleepstats.sleepstats_from_csvs(xml_dir, out_dir,   
                                 subs, sessions, log_dir, outfile)
         return
    
    def export_eventparams(self, evt_name, frequency = None,
                                 xml_dir = None, out_dir = None, subs = 'all', 
                                 sessions = 'all', chan = None, ref_chan = None, 
                                 stage = ['NREM2','NREM3'], grp_name = 'eeg', 
                                 rater = None, cycle_idx = None, 
                                 concat_cycle = True, concat_stage = False, 
                                 keyword = None, segs = None,  
                                 adap_bands = 'Fixed',  
                                 adap_bw = 4, params = 'all', epoch_dur = 30, 
                                 average_channels = False, outfile = True):
        
        # Set up logging
        logger = create_logger('Export params')
        
        # Set input/output directories
        in_dir = self.datapath
        log_dir = self.outpath + '/audit/logs/'
        if not path.exists(log_dir):
            mkdir(log_dir)
        
        # Force evt_name into list, and loop through events    
        if isinstance(evt_name, str):
            evts = [evt_name]
        elif isinstance(evt_name, list):
            evts = evt_name
        else:
            logger.error(TypeError(f"'evt_name' can only be a str or a list, but {type(evt_name)} was passed."))
            return
        for evt_name in evts:
            out_dir = select_output_dirs(self.outpath, out_dir, evt_name)
            xml_dir = select_input_dirs(self.outpath, xml_dir, evt_name)
            
            
            # Check annotations directory exists
            if not path.exists(xml_dir):
                logger.info('')
                logger.critical(f"{xml_dir} doesn't exist. Event detection has not been run or an incorrect event type has been selected.")
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
            fish = FISH(self.rootpath, in_dir, xml_dir, out_dir, log_dir, chan, ref_chan, grp_name, 
                              stage, rater, subs, sessions, self.tracking) 
            fish.line(keyword, evt_name, cat, segs, cycle_idx, frequency, adap_bands, 
                      peaks, adap_bw, params, epoch_dur, Ngo, outfile)
        return
    
    
    def event_dataset(self, chan, evt_name, 
                            xml_dir = None, out_dir = None, subs = 'all', 
                            sessions = 'all', stage = None, 
                            concat_stage = False, concat_cycle = True, 
                            cycle_idx = None, grp_name = 'eeg', 
                            adap_bands = 'Fixed',  params = 'all', outfile=True):
        
        # Set up logging
        logger = create_logger('Event dataset')
        
        # Force evt_name into list, and loop through events    
        if isinstance(evt_name, str):
            evts = [evt_name]
        elif isinstance(evt_name, list):
            evts = evt_name
        else:
            logger.error(TypeError(f"'evt_name' can only be a str or a list, but {type(evt_name)} was passed."))
            return
        for evt_name in evts:
            # Set input/output directories
            in_dir = self.datapath
            log_dir = self.outpath + '/audit/logs/'
            if not path.exists(log_dir):
                mkdir(log_dir)
            if not path.exists(self.outpath + '/datasets/'):
                mkdir(self.outpath + '/datasets/')
            
            if adap_bands in ['Auto','Manual']:
                evt_name = f'{evt_name}_adap'
                self.track(step='fooof', show = False, log = False)
            out_dir = self.outpath + f'/datasets/{evt_name}'
            
            xml_dir = select_input_dirs(self.outpath, xml_dir, evt_name)
            if not path.exists(xml_dir):
                logger.info('')
                logger.critical(f"{xml_dir} doesn't exist. Event detection has not been run or an incorrect event type has been selected.")
                logger.info('Check documentation for how to run a pipeline:')
                logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
                logger.info('-' * 10)
                return
            
            # Format concatenation
            cat = (int(concat_cycle),int(concat_stage),1,1)
            
            # Format chan
            if isinstance(chan, str):
                chan = [chan]
            
            # Default stage
            if stage == None:
                stage = ['NREM2','NREM3']
                
            fish = FISH(self.rootpath, in_dir, xml_dir, out_dir, log_dir, chan, None, grp_name, 
                              stage, subs = subs, sessions = sessions) 
            fish.net(chan, evt_name, adap_bands, params,  cat, cycle_idx, outfile)
        
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
        logger = create_logger('Power spectrum dataset')
        
        # Set input/output directories
        in_dir = self.datapath
        log_dir = self.outpath + '/audit/logs/'
        if not path.exists(log_dir):
            mkdir(log_dir)
        if not path.exists(self.outpath + '/datasets/'):
            mkdir(self.outpath + '/datasets/')
        
        if not out_dir:
            out_dir = select_output_dirs(self.outpath, out_dir, 'powerspec')
        
        xml_dir = select_input_dirs(self.outpath, xml_dir, evt_name='spectrum')
        if not path.exists(xml_dir):
            logger.info('')
            logger.critical(f"{xml_dir} doesn't exist. Event detection has not been run or an incorrect event type has been selected.")
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

        spectrum = Spectrum(in_dir, xml_dir, out_dir, log_dir, chan, 
                            ref_chan = None, grp_name = grp_name, stage = stage, 
                            cat = cat, rater = rater, cycle_idx = cycle_idx, 
                            subs = subs, sessions = sessions) 
        spectrum.powerspec_summary(chan, general_opts, frequency_opts, filter_opts, 
                                   epoch_opts, event_opts, logger)
        
        return
    
        
        