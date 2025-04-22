#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 18 16:29:30 2025

@author: ncro8394
"""


from numpy import (argsort, asarray, concatenate, cumsum, diff, exp, isnan, 
                   logical_and, mean, nan, nanmean, nanstd, sort, var, zeros)
from os import listdir, mkdir, path
from pandas import concat, DataFrame, read_csv
from pathlib import Path
from wonambi import Dataset
from wonambi.detect.spindle import transform_signal
from scipy.integrate import trapz
from scipy.signal import find_peaks, periodogram
from scipy.stats import variation
from sleepecg import detect_heartbeats
from .logs import create_logger
from ..events.sand import detect_ecg_artifacts





class SQUID:
    """
    SQUID: Signal Quality Inference and Detection

    This class orchestrates signal quality assessment over a BIDS-formatted dataset,
    using modular QC pipelines for EEG, EMG, and ECG recordings.

    Parameters
    ----------
    bids_root : str or Path
        Path to the root of the BIDS-formatted dataset.

    filetype : str, optional
        Modality type to analyze ('eeg', 'ecg', 'emg'), by default 'eeg'.

    subjects : list of str, optional
        List of subject IDs to include (e.g., ['01', '02']). If None, all subjects are processed.

    sessions : list of str, optional
        List of session IDs. If None, sessions are inferred automatically.

    logger : logging.Logger, optional
        Custom logger instance. If None, a default logger is created.
    """

    def __init__(self, bids_root, filetype = '.edf', subjects = None, 
                 sessions = None, logger = None):
        self.bids_root = Path(bids_root)
        
        self.datapath = Path(f'{self.bids_root}/sourcedata/')
        if not path.exists(self.datapath):
            self.datapath = Path(f'{self.bids_root}/DATA/')

        if not path.exists(f'{bids_root}/derivatives'):
            mkdir(f'{bids_root}/derivatives')
        if not path.exists(f'{bids_root}/derivatives/QC'):
            mkdir(f'{bids_root}/derivatives/QC') 
        self.outpath = Path(f'{bids_root}/derivatives/QC')
            
        self.filetype = filetype
        self.subjects = subjects
        self.sessions = sessions
        self.logger = create_logger('SQUID')
        self.results = []

    
    def process_all(self, chantype = ['eeg', 'eog', 'emg', 'ecg']):
        """
        Iterates over all subjects and sessions to run quality control.
        """
        subjects = self.subjects or self._get_all_subjects()

        for subj in subjects:
            sessions = self.sessions or self._get_sessions(subj)

            for sess in sessions:
                try:
                    self.logger.debug(f"Processing {subj}, {sess}")
                    self.process_subject(subj, chantype, sess)
                except Exception as e:
                    self.logger.warning(f"Failed to process {subj}, {sess}: {e}")
                    continue

    def process_subject(self, subject_id, chantype, session_id=None):
        """
        Process a single subject/session: infer, choose best, and compute QC.

        Parameters
        ----------
        subject_id : str
        session_id : str or None
        chantype : list
        """
        
        dset = self._load_dataset(subject_id, session_id)
        outdir = self._set_outdir(subject_id, session_id)
    
        # Define function mappings
        candidate_funcs = {
            'eeg': get_candidate_eeg,
            'eog': get_candidate_eog,
            'emg': get_candidate_emg
        }
    
        squid_funcs = {
            'eeg': squid_eeg,
            'eog': squid_eeg,  # assuming same function as eeg
            'emg': squid_emg
        }
    
        for c in chantype:
            if c not in candidate_funcs or c not in squid_funcs:
                self.logger.warning(f"Unsupported chantype: {c}")
                continue
    
            chans = candidate_funcs[c](dset, logger=self.logger)
            if not chans:
                continue
    
            results = []
            for chan in chans:
                try:
                    signal = dset.read_data(chan=[chan]).data[0][0]
                    s_freq = dset.header['s_freq']
                    qc_row = squid_funcs[c](signal, s_freq, name=chan)
                    qc_row['subject'] = subject_id
                    qc_row['session'] = session_id
                    results.append(qc_row)
                except Exception as e:
                    self.logger.warning(f"Failed to process channel {chan}: {e}")
    
            if results:
                savename = f'{subject_id}_{session_id}' if session_id else subject_id
                df = concat(results, ignore_index=True)
                cols_first = ['subject', 'session', 'channel', 
                              'quality_score', 'quality_label']
                cols_rest = [col for col in df.columns if col not in cols_first]
                df = df[cols_first + cols_rest]
                filepath = f'{outdir}/{savename}_qc_report_{c}.csv'
                df.to_csv(filepath, index=False)


        # dset = self._load_dataset(subject_id, session_id)
        # outdir = self._set_outdir(subject_id, session_id)
        
        # for c in chantype:
        #     if c in ('eeg', 'eog'):
        #         if c == 'eeg':
        #             chans = get_candidate_eeg(dset, logger = self.logger)
        #         else:
        #             chans = get_candidate_eog(dset, logger = self.logger)
        #         if chans:
        #             results = []
        #             for chan in chans:
        #                 try:
        #                     signal = dset.read_data(chan=[chan]).data[0][0]  # Adapt as needed
        #                     s_freq = dset.header['s_freq']
        #                     qc_row['subject'] = subject_id
        #                     qc_row['session'] = session_id
        #                     qc_row = squid_eeg(signal, s_freq, name=chan)
        #                     results.append(qc_row)
        #                 except Exception as e:
        #                     self.logger.warning(f"Failed to process channel {chan}: {e}")
                        
        #             filepath = f'{outdir}/qc_report_{c}.csv'
        #             df = concat(results, ignore_index=True)
        #             df.to_csv(filepath, index=False)


        #     if c == 'emg':
        #         chans = get_candidate_emg(dset, logger = self.logger)
        #         if chans:
        #             results = []
        #             for chan in chans:
        #                 try:
        #                     signal = dset.read_data(chan=[chan]).data[0][0]  # Adapt as needed
        #                     s_freq = dset.header['s_freq']
        #                     qc_row['subject'] = subject_id
        #                     qc_row['session'] = session_id
        #                     qc_row = squid_emg(signal, s_freq, name=chan)
        #                     results.append(qc_row)
        #                 except Exception as e:
        #                     self.logger.warning(f"Failed to process channel {chan}: {e}")
                            
        #         filepath = f'{outdir}/qc_report_{c}.csv'
        #         df = concat(results, ignore_index=True)
        #         df.to_csv(filepath, index=False)
            
            # TO DO:
            # if c == 'ecg':
            #     chans = infer_ecg(dset)
            #     if chans:
            #         best = choose_best_emg(chans, dset)
            #         for chan in best:
            #             signal = dset.read_data(chan=chan).data[0][0]  # Adapt as needed
            #             s_freq = dset.get_sampling_rate(chan=chan)
            #             qc_row = squid_emg(signal, s_freq, name=chan)
            #             qc_row['subject'] = subject_id
            #             qc_row['session'] = session_id
            #             self.results.append(qc_row)
                
            #     filepath = f'{outdir}/qc_report_{c}.csv'
            #     self.save_results(filepath)
            

    def _load_dataset(self, subject_id, session_id):
        """
        Load a subject/session dataset. 
        """
        
        dpath = f'{self.datapath}/{subject_id}/{session_id}/eeg'
        if not path.exists(dpath):
            self.logger.warning(f'Directory not found for {subject_id}, {session_id}.')
            return
        
        filename = [x for x in listdir(dpath) if self.filetype in x]
        
        if len(filename) == 0:
            self.logger.warning(f"No {self.filetype} files found for "
                                f"{subject_id}, {session_id}")

        dset = Dataset(f'{dpath}/{filename[0]}')
        
        return dset
    
    def _set_outdir(self, subject_id, session_id):
        """
        Sets the output directory to save the audit. 
        """
        outdir = f'{self.outpath}/{subject_id}/'
        if not path.exists(f'{outdir}'):
            mkdir(f'{outdir}')
        outdir = f'{self.outpath}/{subject_id}/{session_id}/'
        if not path.exists(f'{outdir}'):
            mkdir(f'{outdir}')
            
        return outdir

    def _get_all_subjects(self):
        """
        Infer all subject IDs from the BIDS root.
        """
        subj_dirs = sorted(self.datapath.glob("sub-*"))
        return [d.name.split("/")[-1] for d in subj_dirs]

    def _get_sessions(self, subject_id):
        """
        Infer session IDs from a subject folder, or return [''] if no sessions.
        """
        subj_path = self.datapath / f"{subject_id}"
        ses_dirs = sorted(subj_path.glob("ses-*"))
        if ses_dirs:
            return [d.name.split("/")[-1] for d in ses_dirs]
        else:
            return ['']



def get_candidate_eeg(dset, filt=None, logger=create_logger('Get Candidate EEG')):
    """Return a list of plausible EEG channel names based on known naming conventions.

    Parameters
    ----------
    dset : object
        The dataset object containing signal metadata.
    filt : list of str, optional
        List of string patterns to filter EEG channels explicitly.
    logger : Logger object, optional
        Logger for status and error messages.

    Returns
    -------
    list of str
        Candidate EEG channel names (may be empty if no matches found).
    """
    h = dset.header

    if filt:
        if not isinstance(filt, list):
            raise ValueError("`filt` must be a list")
        eeg = [x for x in h["chan_name"] if any(l in x for l in filt)]

    else:
        eeg = [x for x in h["chan_name"] if any(l in x for l in ['C1', 'C2', 
                                                                 'C3', 'C4', 
                                                                 'C5', 'C6', 
                                                                 'Cz', 'CZ',
                                                                 'F1', 'F2', 
                                                                 'F3', 'F4', 
                                                                 'F5', 'F6', 
                                                                 'Fz', 'FZ',
                                                                 'Fpz'])]

        if len(eeg) == 0:
            eeg = [x for x in h["chan_name"] if 'E' in x and not 
                   any(l in x for l in ['EMG', 'ECG', 'EOG', 'Chin'])]

    return eeg


def get_candidate_eog(dset, logger=create_logger('Get Candidate EOG')):
    """
    Identify potential EOG channels from the dataset header based on naming heuristics.

    Parameters
    ----------
    dset : object
        The dataset object with header information.
    logger : logging.Logger, optional
        Logger for debug messages.

    Returns
    -------
    list
        List of candidate EOG channel names (may be empty if no candidates found).
    """
    h = dset.header

    # First and most obvious EOG label
    eyes = [x for x in h['chan_name'] if 'EOG' in x]

    # Fallback: channels with 'E' but excluding non-EOG types
    if len(eyes) == 0:
        eyes = [x for x in h['chan_name'] if 'E' in x and 
                all(excl not in x for excl in ['EMG', 'ECG'])]

    return eyes


def get_candidate_emg(dset, logger=create_logger('Get Candidate EMG')):
    """
    Identify potential EMG channels from the dataset header based on naming heuristics.

    Parameters
    ----------
    dset : object
        The dataset object with header information.
    logger : logging.Logger, optional
        Logger for debug messages.

    Returns
    -------
    list
        List of candidate EMG channel names.
    """
    h = dset.header

    # Most common EMG indicators
    emg = [x for x in h['chan_name'] if any(l in x for l in ['EMG', 'Chin', 'CHIN'])]

    return emg


def infer_eeg(dset, filt=None, num=1, logger=create_logger('Infer EEG')):
    
    """Return top-N EEG channels based on quality scores using SQUID QC.

    Parameters
    ----------
    dset : object
        Dataset object with header and signal data.
    filt : list of str, optional
        Optional patterns to filter candidate EEGs.
    num : int, optional
        Number of top-quality EEG channels to return.
    logger : Logger, optional
        Logger for messages.

    Returns
    -------
    list of str
        EEG channel names ranked by quality (or None if no suitable channels).
    """
    eeg = get_candidate_eeg(dset, filt=filt, logger=logger)

    if len(eeg) == 0 or len(eeg) > 10:
        logger.error("Unable to determine EEG channels from recording. "
                     "EEG names must be explicitly specified.")
        logger.info("Check documentation for how to specify channels:")
        logger.info("https://seapipe.readthedocs.io/en/latest/index.html")
        logger.info("-" * 10)
        return None

    # Try to find ECG to pass into choose_best_eeg if needed
    ecg = [x for x in dset.header["chan_name"] if any(l in x for l in ["ECG", "EKG"])]
    if len(ecg) > 1:
        ecg = choose_best_ecg(ecg, dset, num=1, logger=logger)

    return choose_best_eeg(eeg, ecg, dset, num=num, logger=logger)
    

def infer_eog(dset, logger=create_logger('Infer EOG')):
    """
    Use candidate EOG channels and quality-check to return best pair.

    Parameters
    ----------
    dset : object
        The dataset object containing signal data and metadata.
    logger : logging.Logger, optional
        Logger for messages.

    Returns
    -------
    list
        List of selected EOG channel names (ideally left + right).
    """
    eyes = get_candidate_eog(dset, logger=logger)

    if len(eyes) == 0 or len(eyes) > 6:
        logger.error('Unable to determine EOG channels from recording. EOG names must be explicitly specified.')
        logger.info('Check documentation for how to specify channels:')
        logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
        logger.info('-' * 10)
        return None

    if len(eyes) > 2:
        loc = [x for x in eyes if any(l in x for l in ['l', 'L'])]
        loc = choose_best_eog(loc, dset, 1, logger) if loc else []

        roc = [x for x in eyes if any(l in x for l in ['r', 'R'])]
        roc = choose_best_eog(roc, dset, 1, logger) if roc else []

        eyes = roc + loc

    return eyes


def infer_emg(dset, logger=create_logger('Infer EMG')):
    """
    Use candidate EMG channels and apply quality checks to select the best one.

    Parameters
    ----------
    dset : object
        The dataset object containing signal data and metadata.
    logger : logging.Logger, optional
        Logger for messages.

    Returns
    -------
    list
        List containing the best EMG channel (or None).
    """
    emg = get_candidate_emg(dset, logger=logger)

    if len(emg) == 0 or len(emg) > 6:
        logger.error('Unable to determine EMG channels from recording. EMG names must be explicitly specified.')
        logger.info('Check documentation for how to specify channels:')
        logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
        logger.info('-' * 10)
        return None

    if len(emg) > 1:
        emg = choose_best_emg(emg, dset, num=1, logger=logger)

    return emg

        
# def infer_eeg(self, dset, filt = None, logger = create_logger('Infer EEG')):
    
#     ''' This function reads a recording file (e.g. edf) and searches for likely
#         EEG channel names. IF the EEG channels are not easly identifiable from
#         the edf, then nothing is returned. IF there are multiple (>2) options,
#         then the channels will be QC'd and the channels with the best quality 
#         will be returned. #NOTE: Central channels will be given highest priority.
#     '''
    
#     h = dset.header
    
#     if filt:
#         if not isinstance(filt, list):
#             raise ValueError('`filt` must be a list')
    
#         eeg = [x for x in h['chan_name'] if any(l in x for l in filt)] 
    
#     else:
#         # 1. Look for central channels first
#         eeg = [x for x in h['chan_name'] if any(l in x for l in ['C1', 'C2', 'C3',
#                                                                  'C4', 'C5', 'C6',
#                                                                  'Cz', 'CZ'])]  
        
#         # 2. If no central channels, look for frontal:
#         if len(eeg) == 0:
#             eeg = [x for x in h['chan_name'] if any(l in x for l in ['F1', 'F2', 'F3',
#                                                                      'F4', 'F5', 'F6',
#                                                                      'Fz', 'FZ'])]          
#         # 3. Look for any channels names beginning with E
#         if len(eeg) == 0:
#             eeg = [x for x in h['chan_name'] if 'E' in x if not 
#                                                     any(l in x for l in ['EMG', 
#                                                                          'ECG',
#                                                                          'EOG',
#                                                                          'Chin',
#                                                                          ])]
               
#     # Check for ECG channel
#     ecg = [x for x in h['chan_name'] if any(l in x for l in ['ECG', 'EKG'])]  
#     if len(ecg) > 1:
#         ecg = self.choose_best_ecg(ecg, dset, num = 1, logger = logger)

    
#     # If still nothing or when all electrode names contain 'E' (e.g. EGI nets) 
#     if len(eeg) == 0 or len(eeg) > 10:
#         logger.error('Unable to determine EEG channels from recording. EEG names must be explicitly specified.')
#         logger.info('Check documentation for how to specify channels:')
#         logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
#         logger.info('-' * 10)
#         return None
#     else:
#         eeg = self.choose_best_eeg(eeg, ecg, dset, 1, logger)
            
#     return eeg 
    
# def infer_eog(self, dset, logger = create_logger('Infer EOG')):
    
#     ''' This function reads a recording file (e.g. edf) and searches for likely
#         EOG channel names. IF the EOG channels are not easly identifiable from
#         the edf, then nothing is returned. IF there are multiple (>2) options,
#         then the channels will be QC'd and the channels with the best quality 
#         will be returned.
#     '''
    
#     h = dset.header
    
#     # First and most obvious EOG name
#     eyes = [x for x in h['chan_name'] if 'EOG' in x]  
    
#     # If EOG are not labelled as 'EOG'
#     if len(eyes) == 0:
#         eyes = [x for x in h['chan_name'] if 'E' in x if not 'EMG' in x 
#                                                       if not 'ECG' in x]    
    
#     # If still nothing or when all electrode names contain 'E' (e.g. EGI nets) 
#     if len(eyes) == 0 or len(eyes) > 6:
#         logger.error('Unable to determine EOG channels from recording. EOG names must be explicitly specified.')
#         logger.info('Check documentation for how to specify channels:')
#         logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
#         logger.info('-' * 10)
#         return None

#     elif len(eyes) > 2:
        
#         loc = [x for x in eyes if any(l in x for l in ['l', 'L'])]
#         loc = self.choose_best_eog(loc, dset, 1, logger)
        
#         roc = [x for x in eyes if any(l in x for l in ['r', 'R'])]
#         roc = self.choose_best_eog(roc, dset, 1, logger)
        
#         eyes = roc + loc
            
#     return eyes    

# def infer_emg(self, dset, logger = create_logger('Infer EMG')):
    
#     ''' This function reads a recording file (e.g. edf) and searches for likely
#         EOG channel names. IF the EOG channels are not easly identifiable from
#         the edf, then nothing is returned. IF there are multiple (>2) options,
#         then the channels will be QC'd and the channels with the best quality 
#         will be returned.
#     '''
    
#     h = dset.header
    
#     # First and most obvious EMG name
#     emg = [x for x in h['chan_name'] if any(l in x for l in ['EMG', 
#                                                              'Chin', 
#                                                              'CHIN'])]  
     
    
#     # If still nothing or when all electrode names contain 'E' (e.g. EGI nets) 
#     if len(emg) == 0 or len(emg) > 6:
#         logger.error('Unable to determine EMG channels from recording. EOG names must be explicitly specified.')
#         logger.info('Check documentation for how to specify channels:')
#         logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
#         logger.info('-' * 10)
#         return None

#     elif len(emg) > 1:
#         emg = self.choose_best_emg(emg, dset, 1, logger)
            
#     return emg  


def choose_best_eeg(chans, ecg, dset, num = 2, logger = create_logger('QC EEG')):
    """
    Choose the best EEG channel(s) based on quality metrics computed using `squid`.

    Parameters
    ----------
    chans : list of str
        List of EEG channel names to evaluate.
    
    dset : object
        The dataset object containing signal data and associated metadata.
    
    num : int, optional, default=1
        The number of top channels to return based on their quality scores.
    
    logger : logging.Logger, optional
        Logger instance for logging any warnings or debug information.
    
    Returns
    -------
    eeg_names : list of str
        List of top `num` EEG channel names based on their quality score.
    """
    ch_dat = dset.read_data(chan = chans)
  
    all_quals = []

    for i, name in enumerate(chans):
        logger.debug(f'Checking quality of {name}')
        s_freq = dset.header['orig']['n_samples_per_record'][chans.index(name)]
        
        quals = SQUID.squid_eeg(ch_dat.data[0][i], s_freq, name)
        quals['channel'] = name  # optionally tag with channel name
        
        if quals.loc['time_not_flatline'] > 50:
            logger.warning(f'Channel {name} is mostly flatlined... Review.')
    
        if quals.loc['time_below_200'] < 0.75:
            logger.warning(f"{quals.loc[i,'time_below_200']*100}% of the channel "
                           f"{name} is > 200 µV. This is implausible. Check if "
                           "the signal contains a lot of artefact or has proper "
                           "scaling.")
        
        all_quals.append(quals)
    
    # Combine all into a single DataFrame
    all_quals_df = concat(all_quals, ignore_index = True)  
    
    # Mask only numeric and relevant columns
    numeric_cols = all_quals_df.select_dtypes(include = 'number').columns
    
    # Apply ranking only to those columns
    winners = all_quals_df[numeric_cols].rank(ascending=False)
    
    # Put back the channel names for clarity
    winners['channel'] = all_quals_df['channel']
    
    #winners = DataFrame([all_quals_df[x].rank() for x in all_quals_df.columns]).T
    idx = winners.sum(axis = 1).nlargest(num).index.to_list()
    eeg_names = [chans[x] for x in idx]
    
    logger.debug(f'Best EEG channels based on auto-QC are: {eeg_names}')
    
    return eeg_names


def choose_best_eog(chans, dset, num = 2, logger = create_logger('QC EOG')):
    """
    Choose the best EOG channel(s) based on quality metrics computed using `squid`.

    Parameters
    ----------
    chans : list of str
        List of EOG channel names to evaluate.
    
    dset : object
        The dataset object containing signal data and associated metadata.
    
    num : int, optional, default=1
        The number of top channels to return based on their quality scores.
    
    logger : logging.Logger, optional
        Logger instance for logging any warnings or debug information.
    
    Returns
    -------
    eog_names : list of str
        List of top `num` EOG channel names based on their quality score.
    """
    ch_dat = dset.read_data(chan = chans)
  
    all_quals = []

    for i, name in enumerate(chans):
        logger.debug(f'Checking quality of {name}')
        s_freq = dset.header['orig']['n_samples_per_record'][chans.index(name)]
        
        quals = SQUID.squid_eeg(ch_dat.data[0][i], s_freq, name)
        quals['channel'] = name  # optionally tag with channel name
        all_quals.append(quals)
    
    # Combine all into a single DataFrame
    all_quals_df = concat(all_quals, ignore_index = True)  
    
    # Mask only numeric and relevant columns
    numeric_cols = all_quals_df.select_dtypes(include = 'number').columns
    
    # Apply ranking only to those columns
    winners = all_quals_df[numeric_cols].rank(ascending=False)
    
    # Put back the channel names for clarity
    winners['channel'] = all_quals_df['channel']
    
    #winners = DataFrame([all_quals_df[x].rank() for x in all_quals_df.columns]).T
    idx = winners.sum(axis = 1).nlargest(num).index.to_list()
    eog_names = [chans[x] for x in idx]
    
    logger.debug(f'Best EOG channels based on auto-QC are: {eog_names}')
    
    return eog_names


def choose_best_emg(chans, dset, num=1, logger=create_logger('QC EMG')):
    """
    Chooses the best EMG channels based on quality metrics computed using `squid_emg`.
    
    Parameters
    ----------
    chans : list of str
        List of channel names to evaluate.
    
    dset : object
        The dataset object containing signal data and associated metadata.
        Assumes `dset.read_data(chan=chans)` can retrieve the EMG signal for each channel.
    
    num : int, optional, default=1
        The number of top channels to return based on their quality scores.
    
    logger : logging.Logger, optional
        Logger instance for logging any warnings or debug information.
    
    Returns
    -------
    best_channels : list of str
        List of `num` top EMG channels, ranked by quality score.
    """
    
    quals = DataFrame(columns=['time_not_flatline', 
                               'stdev', 
                               'power_low', 
                               'power_high', 
                               'broadband_ratio',
                               'snr_proxy', 
                               'cv', 
                               'var_ratio', 
                               'quality_score', 
                               'reasons'
                               ])
    
    emg_dat = dset.read_data(chan=chans)
    
    for i, name in enumerate(chans):
        logger.debug(f'Checking quality of {name}')
        s_freq = dset.header['orig']['n_samples_per_record'][i]
        signal = emg_dat.data[0][i]
        
        # Run the squid_emg function to get quality metrics
        quals_row = SQUID.squid_emg(signal, s_freq, name, logger)
        
        # Append the results to the dataframe
        quals = concat([quals, quals_row], ignore_index=True)
    
    # Rank the channels based on quality score (higher is better)
    winners = quals[['quality_score']].rank(axis=0, ascending=False)
    idx = winners.sum(axis=1).nlargest(num).index.to_list()
    emg_name = [chans[x] for x in idx]
    
    logger.debug(f'Best EMG channels based on auto-QC are: {emg_name}')
    
    return emg_name
        

def choose_best_ecg(chans, dset, num=1, logger=None):
    """
    Choose the best ECG channel(s) based on quality metrics computed using `squid_ecg`.

    Parameters
    ----------
    chans : list of str
        List of ECG channel names to evaluate.
    
    dset : object
        The dataset object containing signal data and associated metadata.
    
    num : int, optional, default=1
        The number of top channels to return based on their quality scores.
    
    logger : logging.Logger, optional
        Logger instance for logging any warnings or debug information.
    
    Returns
    -------
    ecg_names : list of str
        List of top `num` ECG channel names based on their quality score.
    """
    if logger is None:
        logger = create_logger('QC_ECG')

    ecg_dat = dset.read_data(chan=chans)
    results = []

    for i, name in enumerate(chans):
        logger.debug(f'Checking quality of {name}')
        s_freq = dset.header['orig']['n_samples_per_record'][i]
        e = ecg_dat.data[0][i]

        # Compute quality metrics for this channel
        quals = SQUID.squid_ecg(e, s_freq, name, logger)
        results.append(quals)

    # Combine Results into DataFrame
    df_quals = concat(results, ignore_index=True)

    # Rank numeric columns only
    numeric_cols = df_quals.select_dtypes(include='number').columns
    df_ranks = df_quals[numeric_cols].rank(ascending=False)

    # Add overall quality score and rank
    df_quals['quality_score'] = df_ranks.sum(axis=1)
    df_quals['quality_rank'] = df_quals['quality_score'].rank(ascending=False)

    # Select best ECG channel(s)
    top_channels = df_quals.sort_values('quality_score', ascending=False).head(num)
    ecg_names = top_channels['channel'].tolist()

    logger.debug(f'Best ECG channel(s) based on auto-QC: {ecg_names}')

    return ecg_names

    
    
def squid_eeg(data, s_freq, name = None, weights = None, 
          logger = create_logger('SQUID')):
    
    """
        S.Q.U.I.D. – EEG Edition:

        Computes EEG-specific QC metrics and returns a DataFrame row with flags.

        Parameters
        ----------
        data : array-like
            Input timeseries data. Expected to be a nested structure accessible 
            via `ch_dat.data[0][i]`, where `i` indexes the channel. Assumes that 
            associated metadata such as sampling frequency is available globally 
            or via a `dset` object.
            
        Returns
        -------
        quals : pandas.DataFrame
            A dataframe containing the following quality metrics per channel:
                - time_not_flatline : float
                    Proportion of time the signal is not flat (based on moving average).
                - time_above_10 : float
                    Proportion of time the signal amplitude is above 10 µV.
                - time_below_200 : float
                    Proportion of time the signal stays below 200 µV.
                - gini_coeff : float
                    Gini coefficient of signal amplitude distribution.
                - inverse_power_ratio : float
                    Inverse of the ratio of high-frequency to low-frequency power.
                - ecg_artefact : float
                    Measure of ECG contamination based on average peak correlation.
                - ecg_artefact_perc : float
                    Proportion of ECG artefact time with high (r>0.5) correlation.
    
    """

    quals = DataFrame(dtype=float, index=[0])
    reasons = []

    # --- Flatline ---
    e = moving_average_time(abs(data), s_freq, 5)
    flat_ratio = sum(diff(e) == 0) / len(diff(e))
    flat_score = 1 - flat_ratio
    quals.loc[0, 'time_not_flatline'] = flat_score
    quals.loc[0, 'is_flatlined'] = flat_ratio > 0.5
    if flat_ratio > 0.5:
        reasons.append("Flatline > 50%")

    # --- Amplitude constraints ---
    above_10 = logical_and(e >= 0, e <= 10)
    amp_score = round(1 - sum(e > 200) / len(e), 2)
    quals.loc[0, 'time_above_10'] = 1 - sum(above_10) / len(e)
    quals.loc[0, 'time_below_200'] = amp_score
    quals.loc[0, 'has_high_amplitude'] = amp_score < 0.75
    if amp_score < 0.75:
        reasons.append("Signal > 200µV > 25%")

    # --- Gini coefficient ---
    g = gini(e)
    quals.loc[0, 'gini_coeff'] = g

    # --- Power ratio ---
    ratio = transform_signal(data, s_freq, 'moving_power_ratio', 
                             method_opt={'dur': 5,
                                         'freq_narrow': (30, 50),
                                         'freq_broad': (0.5, 5),
                                         'fft_dur': 30,
                                         'step': 15
                                         })
    inv_ratio = 1 / (abs(nanmean(ratio)) + 1e-6)
    quals.loc[0, 'inverse_power_ratio'] = inv_ratio

    # --- ECG contamination ---
    (ecg_artifacts, 
     ecg_template, 
     ecg_corr) = detect_ecg_artifacts(data, sfreq = s_freq)
    (avg_peak, 
     prop_high_corr, 
     total_ecg_area) = compute_ecg_impact_metrics(ecg_artifacts, ecg_corr)
    ecg_score = round(prop_high_corr, 2)
    quals.loc[0, 'ecg_artefact'] = abs(avg_peak)
    quals.loc[0, 'ecg_artefact_perc'] = ecg_score
    quals.loc[0, 'ecg_contaminated'] = prop_high_corr > 0.25
    if prop_high_corr > 0.25:
        reasons.append("ECG contamination > 25%")

    # --- Flags & reasons ---
    quals.loc[0, 'is_bad_channel'] = quals.loc[0, ['is_flatlined', 
                                                   'has_high_amplitude', 
                                                   'ecg_contaminated']].any()
    quals.loc[0, 'reason'] = '; '.join(reasons) if reasons else "OK"
    if name:
        quals.loc[0, 'channel'] = name

    # --- Quality scoring ---
    default_weights = {'time_not_flatline': 0.25,
                       'time_below_200': 0.2,
                       'ecg_clean': 0.2,
                       'low_inv_ratio': 0.2,
                       'gini_balance': 0.15
                       }
    w = weights or default_weights
    total_weight = sum(w.values())

    # Nonlinear transforms
    inv_ratio_scaled = 1 / (1 + exp(-5 * (inv_ratio - 1)))  # sigmoid
    gini_scaled = exp(-10 * (g - 0.5) ** 2)          # bell curve

    quality_score = (w['time_not_flatline'] * flat_score +
                     w['time_below_200'] * amp_score +
                     w['ecg_clean'] * (1 - ecg_score) +
                     w['low_inv_ratio'] * inv_ratio_scaled +
                     w['gini_balance'] * gini_scaled
                     ) / total_weight

    quals.loc[0, 'quality_score'] = round(quality_score, 3)

    # --- Quality label ---
    if quality_score < 0.5:
        label = 'bad'
    elif quality_score < 0.7:
        label = 'medium'
    else:
        label = 'good'
    quals.loc[0, 'quality_label'] = label

    return quals


def squid_emg(signal, s_freq, name=None, weights = None, 
              logger = create_logger('SQUID')):
    """
    S.Q.U.I.D. – EMG Edition: 

    Computes EMG-specific QC metrics and returns a DataFrame row with flags.

    Parameters
    ----------
    signal : array-like
        The raw EMG signal for a given channel, expected to be a 1D array of voltage 
        measurements (in µV).
    
    s_freq : float
        The sampling frequency (Hz) of the signal. Typically derived from the 
        dataset's header or metadata.
    
    name : str, optional
        The name of the channel. Used for logging and reporting purposes.
    
    logger : logging.Logger, optional
        Logger instance for recording warnings or debug information.
    
    Returns
    -------
    quals : pandas.DataFrame
        A dataframe containing the following quality metrics per channel:
            - time_not_flatline : float
                Proportion of time the signal is not flat, based on the difference
                of the moving average (smoothed absolute signal).
            - stdev : float
                Standard deviation of the raw signal, reflecting the signal variability.
            - power_low : float
                Total low-frequency power (20-50 Hz), used as a proxy for signal activity.
            - power_high : float
                Total high-frequency power (70-150 Hz), used to assess the presence of
                muscle bursts or noise.
            - broadband_ratio : float
                Ratio of high-frequency to low-frequency power. High values indicate
                bursts or sharp activity in the signal.
            - snr_proxy : float
                Signal-to-noise ratio, approximated as the ratio of high-frequency power
                to total power.
            - cv : float
                Coefficient of variation of the signal, a measure of burstiness.
            - var_ratio : float
                Variance ratio between maximum windowed variance and mean variance,
                indicating burst characteristics in the signal.
            - quality_score : float
                An overall quality score, calculated as the mean of several key metrics.
            - reasons : str
                Comma-separated string of reasons for signal quality issues (e.g., "flatline",
                "low_stdev", "low_cv").
    
    """

    quals = DataFrame(index=[0], dtype=float)
    reasons = []

    # --- Flatline detection using smoothed rectified signal ---
    e = moving_average_time(signal, s_freq, 5)
    flatness = 1 - sum(diff(e) == 0) / len(e)
    quals.loc[0, 'time_not_flatline'] = flatness

    if flatness < 0.95:
        reasons.append("flatline")
        
    # --- Standard deviation ---
    stdev = nanstd(signal)
    quals.loc[0, 'stdev'] = stdev

    if stdev < 1:
        reasons.append("low_stdev")

    # --- Power spectral density ---
    f, Pxx = periodogram(signal, s_freq)
    low_band = logical_and(f >= 20, f <= 50)
    high_band = logical_and(f >= 70, f <= 150)

    power_low = trapz(Pxx[low_band], f[low_band])
    power_high = trapz(Pxx[high_band], f[high_band])
    total_power = trapz(Pxx, f)

    quals.loc[0, 'power_low'] = power_low
    quals.loc[0, 'power_high'] = power_high
    quals.loc[0, 'broadband_ratio'] = power_high / (power_low + 1e-6)
    quals.loc[0, 'snr_proxy'] = power_high / (total_power + 1e-6)

    if power_high < 1e-2:
        reasons.append("low_high_freq_power")

    # --- Coefficient of variation (burstiness) ---
    cv = variation(signal)
    quals.loc[0, 'cv'] = cv

    if cv < 0.1:
        reasons.append("low_cv")

    # --- Variance ratio: burst detection proxy ---
    win_size = int(s_freq * 2)  # 2-second windows
    n_windows = len(signal) // win_size
    windowed_vars = [var(signal[i*win_size:(i+1)*win_size]) for i in range(n_windows)]

    if windowed_vars:
        var_ratio = max(windowed_vars) / (mean(windowed_vars) + 1e-6)
        quals.loc[0, 'var_ratio'] = var_ratio
        if var_ratio < 2:
            reasons.append("low_burstiness")
    else:
        quals.loc[0, 'var_ratio'] = nan
        reasons.append("short_signal")
        if logger: logger.warning(f"{name}: Too short for variance-based burstiness")
 
    # --- Quality scoring ---
    default_weights = {'time_not_flatline': 0.25,
                       'stdev_scaled': 0.2,
                       'broadband_scaled': 0.15,
                       'cv_scaled': 0.15,
                       'var_ratio_scaled': 0.15,
                       'snr_scaled': 0.1
                       }
    w = weights or default_weights
    
    # --- Nonlinear transforms / scaling ---
    # These help map all features to ~[0, 1], where higher = better
    quals['stdev_scaled'] = 1 / (1 + exp(-1 * (quals['stdev'] - 2)))  # sigmoid centered ~2µV
    quals['broadband_scaled'] = exp(-0.1 * quals['broadband_ratio'])  # penalize very high ratios
    quals['cv_scaled'] = 1 / (1 + exp(-10 * (quals['cv'] - 0.2)))     # sigmoid centered around healthy CV
    quals['var_ratio_scaled'] = 1 / (1 + exp(-2 * (quals['var_ratio'] - 2)))  # ideal burstiness > 2
    quals['snr_scaled'] = 1 / (1 + exp(-10 * (quals['snr_proxy'] - 0.05)))   # sigmoid tuned to low SNR
    
    # Ensure no NaNs interfere
    score_fields = list(default_weights.keys())
    valid_scores = [w[k] * quals.loc[0, k] for k in score_fields if not isnan(quals.loc[0, k])]
    valid_weights = [w[k] for k in score_fields if not isnan(quals.loc[0, k])]
    
    # --- Final summary score (mean of key features) ---
    quals.loc[0, 'quality_score'] = (sum(valid_scores) / sum(valid_weights) 
                                     if valid_scores else nan)

    if name:
        quals.loc[0, 'channel'] = name
        
    # --- Quality label ---
    if quals.loc[0, 'quality_score'] < 0.5:
        label = 'bad'
    elif quals.loc[0, 'quality_score'] < 0.7:
        label = 'medium'
    else:
        label = 'good'
    quals.loc[0, 'quality_label'] = label

    # --- Add flags as string list ---
    quals.loc[0, 'reasons'] = ", ".join(reasons) if reasons else ""

    return quals
    
def squid_ecg(data, s_freq, name=None, logger=None):
    """
    S.Q.U.I.D. – ECG Edition: 
    
    Computes a series of quality metrics for ECG signals, including:
    - Heart rate (BPM)
    - Standard deviation
    - Flatline detection
    - Power spectral density (PSD) analysis

    Parameters
    ----------
    data : array-like
        ECG signal data for a given channel.
    
    s_freq : float
        Sampling frequency (Hz) of the ECG signal.
    
    name : str, optional
        The name of the channel (for logging purposes).
    
    logger : logging.Logger, optional
        Logger instance for logging any warnings or debug information.

    Returns
    -------
    quals : pandas.DataFrame
        A dataframe containing the following quality metrics per channel:
            - hr_qual : int
                Heart rate quality (0, 1, 2).
            - stdev : float
                Standard deviation of the absolute ECG signal.
            - time_not_flatline : float
                Proportion of time the signal is not flat.
            - PSD_qual : int
                Quality of the power spectral density (0, 1, 2).
    """
    quals = DataFrame(data=None, columns=['hr_qual', 'stdev', 'time_not_flatline', 'PSD_qual'])

    # Heart Rate Check
    beats = detect_heartbeats(data, s_freq)
    hr = len(beats) / (len(data) / (s_freq * 60))

    if 40 < hr < 120:
        quals['hr_qual'] = 2
        if logger:
            logger.debug(f'Heart rate ({round(hr)} bpm) is plausible.')
    elif 120 <= hr < 200:
        quals['hr_qual'] = 1
        if logger:
            logger.warning(f'Heart rate ({round(hr)} bpm) is plausible but high.')
    else:
        quals['hr_qual'] = 0
        if logger:
            logger.warning(f'Heart rate ({round(hr)} bpm) is NOT plausible.')

    # Standard Deviation
    quals['stdev'] = nanstd(abs(data))

    # Flatline Check
    quals['time_not_flatline'] = 1 - (sum(diff(data) == 0) / len(diff(data)))

    # PSD Quality
    f_p, Pxx_spec = periodogram(data, s_freq)
    Pxx_spec_sm = concatenate((
        zeros(2000),
        moving_average(Pxx_spec, 4000),
        zeros(1999)
    ))
    top = int((len(data) / s_freq) * 20)  # up to ~20 Hz
    peaks = find_peaks(Pxx_spec_sm[:top], prominence=800)[0]

    if len(peaks) > 4:
        quals['PSD_qual'] = 2
    elif len(peaks) > 2:
        quals['PSD_qual'] = 1
    else:
        quals['PSD_qual'] = 0
        if logger:
            logger.warning(f'PSD looks unusual for channel: {name}')
    
    return quals


def gather_qc_reports(qc_root, modality='eeg'):
    """
    Gather all QC reports for a given modality into a single DataFrame.
    
    Parameters
    ----------
    qc_root : str or Path
        Root path to the derivatives/QC folder.
    modality : str
        Modality to search for (e.g., 'eeg', 'emg').
        
    Returns
    -------
    pd.DataFrame
        Combined DataFrame of all qc_report_<modality>.csv files.
    """
    qc_root = Path(qc_root)
    all_reports = []
    
    # Look for files like sub-*/ses-*/qc_report_<modality>.csv
    for csv_file in qc_root.glob("sub-*/ses-*/sub-*_ses-*qc_report_{}.csv".format(modality)):
        try:
            df = read_csv(csv_file, sep=None, engine='python')  # Handles tabs or commas
            df['file_path'] = str(csv_file)
            all_reports.append(df)
        except Exception as e:
            print(f"Could not read {csv_file}: {e}")
    
    if all_reports:
        combined_df = concat(all_reports, ignore_index=True)
        return combined_df
    else:
        return DataFrame()  # Return empty dataframe if none found
    
       
        
#%% Quality Control Checks 

def gini(x, w = None):
    
    ''' Calculates the Gini coefficient:
                    a measure of statistical dispersion calculated by comparing 
                    the actual Lorenz curve to the diagonal line of equality.
    '''

    x = asarray(x)
    if w is not None:
        w = asarray(w)
        sorted_indices = argsort(x)
        sorted_x = x[sorted_indices]
        sorted_w = w[sorted_indices]
        # Force float dtype to avoid overflows
        cumw = cumsum(sorted_w, dtype=float)
        cumxw = cumsum(sorted_x * sorted_w, dtype=float)
        return (sum(cumxw[1:] * cumw[:-1] - cumxw[:-1] * cumw[1:]) / 
                (cumxw[-1] * cumw[-1]))
    else:
        sorted_x = sort(x)
        n = len(x)
        cumx = cumsum(sorted_x, dtype=float)
        # The above formula, with all weights equal to 1 simplifies to:
        return (n + 1 - 2 * sum(cumx) / cumx[-1]) / n


def compute_ecg_impact_metrics(ecg_artifacts, ecg_corr, percentile_threshold=99):
    """
    Compute 3 summary metrics for ECG artifact contamination.

    Parameters:
    - ecg_artifacts: array of sample indices with detected ECG events
    - ecg_corr: full correlation time series from template matching
    - percentile_threshold: for "high correlation" cutoff (default=99th percentile)

    Returns:
    - avg_peak_corr: mean ECG correlation at detected artifact peaks (Metric 1)
    - high_corr_ratio: % of samples with ECG-like correlation above threshold (Metric 2)
    - ecg_area: total correlation "area" across detected ECG events (Metric 3)
    """

    # Metric 1: Avg correlation at detected artifact peaks
    avg_peak_corr = mean(ecg_corr[ecg_artifacts])

    # Metric 2: Ratio of signal strongly affected by ECG
    threshold = 0.5
    high_corr_ratio = mean(ecg_corr > threshold)

    # Metric 3: Area under the correlation curve at ECG peaks
    ecg_area = sum(ecg_corr[ecg_artifacts])

    return avg_peak_corr, high_corr_ratio, ecg_area    



## Other functions for SQUID:   
def moving_average_time(a, s_freq, win=5):
    n = int(win*s_freq*60)
    ret = cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

def moving_average(a, win=100):
    n = int(win)
    ret = cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n    
    
