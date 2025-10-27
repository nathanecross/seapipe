#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  7 15:47:20 2025

@author: ncro8394
"""

from numpy import (angle, arange, array, average, concatenate, degrees, diff, exp,
                   histogram,  inf, isnan, linspace, max, min, nan, nanmean,
                   nanmedian, nanstd, pi, random, sort, vstack, where, zeros)
from pandas import DataFrame, Series
from collections import Counter, defaultdict
from copy import deepcopy
from os import listdir, mkdir, path
import pywt
import shutil
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.signal import (butter, firwin, filtfilt, find_peaks, detrend, 
                          hilbert, morlet2, periodogram, resample, welch) 
from sklearn.cluster import DBSCAN
from wonambi import Dataset
from wonambi.attr import Annotations
from wonambi.detect import consensus, DetectSpindle
from wonambi.trans import fetch
from ..utils.logs import create_logger, create_logger_outfile
from ..utils.load import (load_channels, load_adap_bands, load_sessions, 
                          rename_channels, read_manual_peaks)
from ..utils.misc import remove_duplicate_evts, merge_epochs


class clam:

        def __init__(self, rootpath, rec_dir, xml_dir, out_dir, chan, 
                     grp_name, stage, rater = None, subs = 'all', 
                     sessions = 'all', tracking = None):
            
            self.rootpath = rootpath
            self.rec_dir = rec_dir
            self.xml_dir = xml_dir
            self.out_dir = out_dir
            self.chan = chan
            self.grp_name = grp_name
            self.stage = stage
            self.rater = rater
            
            self.subs = subs
            self.sessions = sessions
            
            if tracking == None:
                tracking = {'cluster':{}}
            self.tracking = tracking

        def clustering(self, evt_type = 'spindle', 
                             freq_bands = {'SWA': (0.5, 4), 'Sigma': (10, 15)},
                             filetype = '.edf', grp_name = 'eeg',
                             logger = create_logger('Event clustering')):
            
            
            
            ### 0.a Set up logging
            tracking = self.tracking
            flag = 0
            
            logger.info('')            
            logger.debug(rf"""Commencing clustering and fluctuations pipeline...
       
                                  __.--:__:--.__
                             _.--°      |       °--._
                            |    \      |      /     ;
                            \     \     |     /     /
                             `\    \    |    /    /'
                               `\   \   |   /   /'
                                 `\  \  |  /  /'
                                _.-`\ \ | / /'-._
                               (_____`\\|//'_____)      
                                       `-'
                  
                    Clustering and Low-frequency Activity Modulations 
                    (C.L.A.M)   
                    
                    Event type: {evt_type}
                                                """,)
            
            
            ### 1. First we check the directories
            # a. Check for output folder, if doesn't exist, create
            if path.exists(self.out_dir):
                    logger.debug("Output directory: " + self.out_dir + " exists")
            else:
                mkdir(self.out_dir)
            
            # b. Check input list
            subs = self.subs
            if isinstance(subs, list):
                None
            elif subs == 'all':
                    subs = listdir(self.rec_dir)
                    subs = [p for p in subs if not '.' in p]
            else:
                logger.error("'subs' must either be an array of subject ids or = 'all' ")       
            
            ### 2. Begin loop through dataset
           
            # a. Begin loop through participants
            subs.sort()
            for i, sub in enumerate(subs):
                if not sub in self.tracking['cluster'].keys():
                    self.tracking['cluster'][sub] = {}
                
                # b. Begin loop through sessions
                flag, sessions = load_sessions(sub, self.sessions, self.rec_dir, flag, 
                                         logger, verbose=2)   
                for v, ses in enumerate(sessions):
                    logger.info('')
                    logger.debug(f'Commencing {sub}, {ses}')
                    if not ses in self.tracking['cluster'][sub].keys():
                        self.tracking['cluster'][sub][ses] = {} 
        
                    
                    ## c. Load recording
                    rdir = self.rec_dir + '/' + sub + '/' + ses + '/eeg/'
                    try:
                        edf_file = [x for x in listdir(rdir) if x.endswith(filetype)]
                        dset = Dataset(rdir + edf_file[0])
                    except:
                        logger.warning(f' No input {filetype} file in {rdir}')
                        flag+=1
                        break
                    
                    ## d. Load annotations
                    xdir = self.xml_dir + '/' + sub + '/' + ses + '/'
                    if not path.exists(xdir):
                        logger.warning(f'{evt_type} has not been detected for'
                                       f'{sub}, {ses}. Skipping..')
                        flag += 1
                        continue
                    xml_file = [x for x in listdir(xdir) if x.endswith('.xml')][0]
                    # Copy annotations file before beginning
                    if not path.exists(self.out_dir):
                        mkdir(self.out_dir)
                    if not path.exists(self.out_dir + '/' + sub):
                        mkdir(self.out_dir + '/' + sub)
                    if not path.exists(self.out_dir + '/' + sub + '/' + ses):
                        mkdir(self.out_dir + '/' + sub + '/' + ses)
                    outpath = self.out_dir + '/' + sub + '/' + ses + '/'
                    backup_file = (f'{outpath}{sub}_{ses}.xml')
                    if not path.exists(backup_file):
                        shutil.copy(xdir + xml_file, backup_file)

                    # Read annotations file
                    annot = Annotations(backup_file, rater_name=self.rater)
                    
                    ## e. Channel setup 
                    pflag = deepcopy(flag)
                    flag, chanset = load_channels(sub, ses, self.chan, self.chan,
                                                  flag, logger)
                    if flag - pflag > 0:
                        logger.warning(f'Skipping {sub}, {ses}...')
                        flag += 1
                        continue
                    
                    newchans = rename_channels(sub, ses, self.chan, logger)
                    
                    for c, ch in enumerate(chanset):
                        
                        # 5.b Rename channel for output file (if required)
                        if newchans:
                            fnamechan = newchans[ch]
                        else:
                            fnamechan = ch
                            
                        stagename = '-'.join(self.stage)
                        
                        
            
                        # ---- Analys event clustering ----
                        logger.debug(f'Analysing cluster metrics for {evt_type}, {ch}')
                        stats = cluster_metrics(annot, ch, evt_type, 
                                                   self.stage, grp_name)
                        
                        # Save stats to file
                        if not isinstance(stats, dict):
                            flag += 1
                            continue
                        else:
                            df = DataFrame([stats])
                            df.insert(0, 'sub', sub)  
                            df.insert(1, 'ses', ses)   
                            logger.debug(f'Saving cluster metrics: {evt_type}')
                            df.to_csv(f'{outpath}/{sub}_{ses}_{fnamechan}_{stagename}_{evt_type}_clustering.csv')
                            
                            
                            
                        # ---- Analyse low-frequency fluctuations ----
                        logger.debug('Analysing low frequency fluctuations for '
                                     f"{', '.join([band for band in freq_bands])} "
                                     f'in channel: {fnamechan}: {chanset[ch]}')
                        segs = fetch(dset, annot, cat = (0,1,1,1), 
                                     stage = self.stage,
                                      reject_artf = True)
                        segs.read_data(chan = [ch], ref_chan = chanset[ch])
                        baseline_data = segs[0]['data'].data[0][0]
                        sampling_rate = segs[0]['data'].s_freq
                        baseline_data = baseline_data[0:int(sampling_rate*100)]
                        
                        segs = fetch(dset, annot, cat = (0,0,0,0), 
                                     stage = self.stage,
                                      reject_artf = True)
                        segs.read_data(chan = [ch], ref_chan = chanset[ch])
                        
                        psd = {band:[] for band in freq_bands}
                        dur = {band:[] for band in freq_bands}
                        phases = {band:[] for band in freq_bands}
                        trace = {band:zeros((2,1)) for band in freq_bands}
                        results = None
                        for i, seg in enumerate(segs):
                            signal_data = seg['data'].data[0][0]
                            if len(signal_data)/sampling_rate > 300:
                                results = extract_infraslow_spectral_profile(signal_data, 
                                                                             sampling_rate,
                                                                             baseline_data,
                                                                             band_defs=freq_bands)
                                
                                for band in freq_bands:
                                    psd[band].append(results[band]['psd'])
                                    dur[band].append(results[band]['duration'])

                                    ## ---- PHASE ANGLE OF SPINDLE OCCURENCES ----
                                    # Get spindle events
                                    start_seg = seg['data'].axis['time'][0].min()
                                    end_seg = seg['data'].axis['time'][0].max()
                                    events = [x for x in annot.get_events(evt_type, chan = f'{ch} ({grp_name})') if 
                                              x['start'] > start_seg and x['end'] < end_seg]
                                    evts_mid = [ (x['start'] - start_seg) + 
                                                ((x['end'] - start_seg) - (x['start'] - start_seg))/2
                                                for x in events]
                                    #evts_idx = round(array(evts_mid) * resample_rate).astype(int)
                                    evts_idx = [int(x*2) for x in evts_mid]
                                    
                                    # Calculate phase angle of spindle mid-point
                                    sigma_trace = results[band]['smoothed_power_trace']
                                    sigma_trace = sigma_trace - nanmean(sigma_trace)
                                    
                                    smoothed_signal = Series(sigma_trace).rolling(window=5, 
                                                                                center=True,
                                                                                min_periods=1).mean().to_numpy()
                                    
                                    smoothed_signal = resample(smoothed_signal, 
                                                               len(sigma_trace))
                                    analytic_signal = hilbert(smoothed_signal)
                                    instantaneous_phase = angle(analytic_signal)
                                    event_phases = instantaneous_phase[evts_idx]
                                    phases[band].extend(event_phases.tolist())
                                    
                                    # Save filtered trace
                                    l = linspace(int(start_seg),
                                                 int(end_seg), 
                                                 int(end_seg)- int(start_seg)+1)
                                    
                                    signal_1s = smoothed_signal[::2]
                                    
                                    if len(signal_1s) > len(l):
                                        signal_1s = signal_1s[:len(l)]
                                    elif len(signal_1s) < len(l):
                                        l = l[:len(signal_1s)]
                                    trace[band] = concatenate((trace[band], 
                                                               vstack((l, signal_1s))),
                                                              axis = 1)
                                    
                                # End the analysis (only analyse the first segment >300 sec)    
                                break
                                
      
                                 
                        # Peak of averaged power spectrum
                        if results:
                            stats = {}
                            for band in freq_bands:
                                arrays = array(psd[band])           
                                weights = array(dur[band])
                                weighted_avg = average(arrays, 
                                                       axis=0, 
                                                       weights=weights) if len(weights)>0 else average(arrays, axis=0)
                                mean_psd = weighted_avg / nanmean(weighted_avg)
                                             
                                # Gaussian fit 
                                freqs_low = results[band]['infraslow_freqs']
                                (peak_freq, 
                                 sigma, 
                                 avg_peak_power) = fit_gaussian(
                                                                freqs_low[(freqs_low > 0.01) & (freqs_low < 0.05)],
                                                                mean_psd[(freqs_low > 0.01) & (freqs_low < 0.05)]
                                                                )
                                logger.debug(f"Chan: {fnamechan}: {chanset[ch]} - "
                                             f"{band} infraslow peak: {round(peak_freq,3)}; "
                                             f"Peak power: {round(avg_peak_power,3)}")
                                
                                stats[f'{band}_peak_freq'] = round(peak_freq,3)
                                stats[f'{band}_avg_peak_power'] = round(avg_peak_power,3)
                                
                                # Histogram of events
                                bins = linspace(-pi, pi, 19)
                                hist, bin_edges = histogram(phases[band], bins=bins)
                                bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
                                hist_df = DataFrame(hist).transpose()
                                hist_df.columns = [int(round(degrees(x))) for x in bin_centers]
                                hist_df.to_csv(f'{outpath}/{sub}_{ses}_{fnamechan}_{stagename}_{evt_type}_{band}_coupling.csv')
                                
                                # PSD of LowFreq-Fluctuation
                                psd_df = DataFrame(mean_psd).transpose()
                                psd_df.columns = [f'{round(f,3)} Hz' for f in freqs_low]
                                psd_df.to_csv(f'{outpath}/{sub}_{ses}_{fnamechan}_{stagename}_{band}_fluctuations_psd.csv')
                                
                                # Filtered trace
                                ts_df = DataFrame([trace[band][1]], 
                                                  columns = trace[band][0])
                                ts_df.to_csv(f'{outpath}/{sub}_{ses}_{fnamechan}_{stagename}_{band}_fluctuations_timeseries.csv')
    
                                
                            # Save stats to file
                            df = DataFrame([stats])
                            df.insert(0, 'sub', sub)  
                            df.insert(1, 'ses', ses)   
                            bandnames = '-'.join([band for band in freq_bands])
                            df.to_csv(f'{outpath}/{sub}_{ses}_{fnamechan}_{stagename}_{bandnames}_fluctuation_stats.csv')
                        else:
                            logger.warning(f'No valid segments found for {sub}. {ses}. Skipping...')
                            flag +=1

            ### 3. Check completion status and print
            if flag == 0:
                logger.debug('Event clustering analysis finished without error.')  
            else:
                logger.warning(f'Event clustering analysis finished with {flag} WARNINGS. See log for details.')
            
            #self.tracking = tracking   ## TO UPDATE - FIX TRACKING
            
            return 
            
            
            
def cluster_metrics(annot, chan, evt_name = 'spindle', stage = None, 
                       grp_name = 'eeg', 
                       logger = create_logger('Event clustering')):
    
    
    """
    Compute temporal clustering metrics for detected sleep events (e.g., spindles).

    This function quantifies the degree to which events (such as spindles) are temporally clustered 
    within a recording, including the coefficient of variation (CV), variance-to-mean ratio (VMR), 
    and DBSCAN-based cluster statistics. It also performs permutation-based null testing to evaluate 
    statistical significance of clustering.

    Parameters
    ----------
    annot : wonambi.Annotations
        Wonambi annotations object containing event and staging information.
    chan : str
        Name of the EEG channel to extract events from.
    evt_name : str, optional
        Name of the event type to analyze (default is 'spindle').
    stage : list of str or None, optional
        List of sleep stages (e.g., ['NREM2', 'NREM3']) to restrict event analysis to.
        If None, all stages are used.
    grp_name : str, optional
        Name of the signal group in the annotations (default is 'eeg').
    logger : logging.Logger, optional
        Logger instance for progress and warning messages.

    Returns
    -------
    stats : dict
        Dictionary of summary statistics describing temporal clustering of events:

        - 'cv' : float
            Coefficient of variation of inter-event intervals.
        - 'cv_p_value' : float
            P-value from permutation test for CV.
        - 'VMR' : float
            Variance-to-mean ratio of rolling event counts.
        - 'VMR_p_value' : float
            P-value from Monte Carlo test for VMR.
        - 'num_clusters_n' : int
            Number of DBSCAN-detected clusters.
        - 'average_cluster_size_n' : float
            Average number of events per cluster.
        - 'max_cluster_size_n' : int
            Size of the largest cluster.
        - 'min_cluster_size_n' : int
            Size of the smallest cluster.
        - 'clusters_per_60s' : float
            Number of clusters per minute of analyzed time.
        - 'ave_within_cluster_freq' : float
            Average frequency (Hz) of events within clusters.
        - 'mean_iei_s' : float
            Mean inter-event interval (in seconds) within clusters.
        - 'median_iei_s' : float
            Median inter-event interval (in seconds) within clusters.
        - 'min_interval_s' : float
            Minimum interval (in seconds) between clusters.
        - 'mean_interval_s' : float
            Mean interval (in seconds) between clusters.
        - 'max_interval_s' : float
            Maximum interval (in seconds) between clusters.
        - 'min_cluster_duration_s' : float
            Shortest cluster duration (in seconds).
        - 'avg_cluster_duration_s' : float
            Average cluster duration (in seconds).
        - 'max_cluster_duration_s' : float
            Longest cluster duration (in seconds).

    Notes
    -----
    - Events must be pre-detected and stored in the provided annotations object.
    - Clustering uses DBSCAN with fixed parameters (eps=20, min_samples=3).
    - If fewer than 3 events are found, clustering statistics will not be returned.
    """


    # Stages
    stage_key = {'Wake':0,
                 'NREM1':1,
                 'NREM2':2,
                 'NREM3':3,
                 'REM':5,
                 'Artefact':0,
                 'Undefined':0,
                 'Unknown':0}
    
    epochs = annot.get_epochs()
    leng = epochs[-1]['end']
    
    hypno = zeros(leng)
    
    for epoch in epochs:
        start = int(epoch['start'])
        end = int(epoch['end'])
        
        hypno[start:end] = stage_key[epoch['stage']]

    merged = merge_epochs(epochs)
    if stage:
        merged = [x for x in merged if x['stage'] in stage]
    
    # Get Spindles
    chan_full = chan + ' (' + grp_name + ')'
    events = annot.get_events(evt_name, chan = chan_full, stage = stage)
    if len(events) == 0:
        logger.warning(f'No {evt_name} events found in Annotations file.')
        return 'error'
    evt_array = zeros(leng)
    
    for ev in events:
        start = int(ev['start'])
        end = int(ev['end'])
        evt_array[start:end] = 1
    
    evt_array_trim = array([0])
    for epoch in merged:
            segment = evt_array[epoch['start']:epoch['end']]
            evt_array_trim = concatenate((evt_array_trim,segment))
    
    num_permutations = 1000
    # Calculate temporal Coefficient of Variation
    event_times = where(evt_array_trim == 1)[0]  # Indices of events
    inter_event_times = diff(event_times)
    
    if len(inter_event_times) > 1:
        cv = nanstd(inter_event_times) / nanmean(inter_event_times)
    
        # Permutation test: Shuffle event locations to create a null distribution
        cv_null = []
        event_indices = sort(random.choice(arange(100), 
                                                 size=10, 
                                                 replace=False))  # Random non-adjacent events
        for _ in range(num_permutations):
            shuffled_events = random.choice(arange(100), 
                                               size=len(event_indices),
                                               replace=False)
            shuffled_iet = diff(sort(shuffled_events))
            cv_null.append(nanstd(shuffled_iet) / nanmean(shuffled_iet))
        
        # Compute p-value: fraction of null CVs greater than observed CV
        cv_p_value = sum(array(cv_null) >= cv) / num_permutations
        #print(f"Observed CV of inter-event times: {cv:.2f}, p-value: {cv_p_value:.3f}")
    
    
    # Calculate the observed VMR (Variance-to-Mean Ratio)
    window_size = 30
    rolling_counts = Series(evt_array_trim).rolling(window=window_size).sum()
    rolling_counts = rolling_counts[::5]
    observed_vmr = rolling_counts.var() / rolling_counts.mean()
    
    # Monte Carlo Simulation
    vmr_null = []
    for _ in range(num_permutations):
        shuffled_series = zeros(100)
        shuffled_events = random.choice(arange(100), size=len(event_indices), replace=False)
        shuffled_series[shuffled_events] = 1
        shuffled_counts = Series(shuffled_series).rolling(window=window_size).sum()
        shuffled_counts = shuffled_counts[::5]
        vmr_null.append(shuffled_counts.var() / shuffled_counts.mean())
    
    # Compute p-value
    vmr_p_value = sum(array(vmr_null) >= observed_vmr) / num_permutations
    #print(f"Observed VMR: {observed_vmr:.2f}, p-value: {vmr_p_value:.3f}")
    

    #### Find clusters and get parameters    
    if len(event_times) >= 3:  # DBSCAN needs at least min_samples
        db = DBSCAN(eps=20, min_samples=3, metric="euclidean").fit(event_times.reshape(-1, 1))
        labels = db.labels_
        
        # ---- Get unique clusters (excluding noise: label -1) ----
        valid_cluster_labels = set(labels) - {-1}
        num_clusters = len(valid_cluster_labels)
        
        # ---- Average cluster size ----
        cluster_counts = Counter(labels)
        if -1 in cluster_counts:  # Remove noise label (-1)
            del cluster_counts[-1]
        average_cluster_size = nanmean(list(cluster_counts.values()))

        # ---- Clusters per unit time (e.g., per 30 seconds) ----
        time_range = sum([x['end'] - x['start'] for x in merged])
        clusters_per_60s = num_clusters / (time_range / 60)
        
        
        # ---- Group events by cluster ----
        clusters = defaultdict(list)
        for t, label in zip(event_times.flatten(), labels):
            if label != -1:  # skip noise
                clusters[label].append(t)
                
        # ---- Min/Max cluster size ----
        cluster_sizes = [len(times) for times in clusters.values() if len(times) > 0]
        if len(cluster_sizes) > 0:
            max_cluster_size = max(cluster_sizes)
            min_cluster_size = min(cluster_sizes)   
        else:
            max_cluster_size = nan
            min_cluster_size = nan
        
        
        # ---- Calculate within-cluster frequency ----
        cluster_freqs = {}
        for label, times in clusters.items():
            times = sorted(times)
            duration = times[-1] - times[0]
            n_events = len(times)
        
            if duration > 0:  # avoid division by zero
                freq = n_events / duration  # events per second
                cluster_freqs[label] = freq
            else:
                cluster_freqs[label] = nan  # maybe a single event or nearly simultaneous
       
        # Summary stats 
        valid_freqs = [f for f in cluster_freqs.values() if not isnan(f)]
        average_within_cluster_freq = nanmean(valid_freqs)
        
        
        # ---- Calculate IEIs per cluster ----
        cluster_ieis = {}  # store IEIs per cluster
        all_ieis = []      # collect all IEIs if you want overall stats
        
        for label, times in clusters.items():
            sorted_times = sorted(times)
            
            if len(sorted_times) >= 2:
                ieis = diff(sorted_times)
                cluster_ieis[label] = ieis
                all_ieis.extend(ieis)
            else:
                cluster_ieis[label] = array([])  # cluster with only one event
        
        # Summary stats
        mean_iei = nanmean(all_ieis)
        median_iei = nanmedian(all_ieis)
        
        
        # ---- Between cluster intervals & cluster duration ----
        # Build (start, end) pairs for each cluster
        cluster_bounds = []
        cluster_durations = []
        for label, times in clusters.items():
            if len(times) > 1:
                start = min(times)
                end = max(times)
                cluster_bounds.append((start, end))
                duration = max(times) - min(times)
                cluster_durations.append(duration)

        cluster_bounds.sort() # Sort clusters by start time
        
        # Compute between-cluster intervals
        between_cluster_intervals = []
        for i in range(len(cluster_bounds) - 1):
            current_end = cluster_bounds[i][1]
            next_start = cluster_bounds[i+1][0]
            interval = next_start - current_end
            between_cluster_intervals.append(interval)
        
        # Cluster duration
        if len(cluster_durations) > 0:
            avg_cluster_duration = nanmean(cluster_durations)
            min_cluster_duration = min(cluster_durations)
            max_cluster_duration = max(cluster_durations)
        else:
            avg_cluster_duration = min_cluster_duration = max_cluster_duration = nan
        
        # Between cluster intervals
        if len(between_cluster_intervals) > 0:
            min_interval = min(between_cluster_intervals)
            mean_interval = nanmean(between_cluster_intervals)
            max_interval = max(between_cluster_intervals)
        else:
            min_interval = mean_interval = max_interval = nan

        stats = {'cv': round(cv, 2),
                 'cv_p_value': round(cv_p_value, 3),
                 'VMR': round(observed_vmr,2),
                 'VMR_p_value': round(vmr_p_value, 3),
                 'num_clusters_n': num_clusters,
                 'average_cluster_size_n': round(average_cluster_size,2),
                 'max_cluster_size_n': round(max_cluster_size,2),
                 'min_cluster_size_n': round(min_cluster_size,2),
                 'clusters_per_60s': round(clusters_per_60s,2),
                 'ave_within_cluster_freq': round(average_within_cluster_freq,4),
                 'mean_iei_s': round(mean_iei, 2),
                 'median_iei_s': round(median_iei, 2),
                 'min_interval_s': round(min_interval, 2),
                 'mean_interval_s': round(mean_interval, 2),
                 'max_interval_s': round(max_interval, 2),
                 'min_cluster_duration_s': round(min_cluster_duration, 2),
                 'avg_cluster_duration_s': round(avg_cluster_duration, 2),
                 'max_cluster_duration_s': round(max_cluster_duration, 2),
                }
        
    else:
        logger.warning("Not enough events in segment to compute clustering.")
        stats = 'error'
        
    return stats            
            


def extract_infraslow_spectral_profile(eeg_signal, fs, baseline_data,
                                        wavelet='cmor4.0-1.0',
                                        band_defs={'SWA': (0.5, 4), 'Sigma': (10, 15)},
                                        smoothing_window_sec = 4,
                                        infraslow_freq_range = (0.001, 0.12),
                                        infraslow_step = 0.001,
                                        downsample_step_sec = 0.5):
    
    """
    Extract infraslow power fluctuations in EEG band-specific power using a two-stage wavelet decomposition.

    This function computes smoothed and normalized power time courses for specified frequency bands
    (e.g., SWA, sigma), then applies a second wavelet transform to quantify infraslow fluctuations
    in these band power signals.
    
    This function is based upon the paper: 
        Lecci et al. (2017) - Science Advances;  
        'Coordinated infraslow neural and cardiac oscillations mark fragility 
        and offline periods in mammalian sleep.'

    Parameters
    ----------
    eeg_signal : array_like
        1D array of EEG signal values (time series).
    fs : float
        Sampling frequency of the EEG signal in Hz.
    baseline_data : array_like
        EEG data (typically from a baseline period) used to compute average band power for normalization.
    wavelet : str, optional
        Wavelet name used for CWT (default is 'cmor4.0-1.0').
    band_defs : dict, optional
        Dictionary defining frequency bands to analyze, with names as keys and (fmin, fmax) tuples as values.
    smoothing_window_sec : float, optional
        Duration of the moving average smoothing window in seconds (default is 4).
    infraslow_freq_range : tuple, optional
        Range of infraslow frequencies (in Hz) for the second-stage wavelet transform (default is (0.001, 0.12)).
    infraslow_step : float, optional
        Step size (Hz) for infraslow frequency resolution (default is 0.001).
    downsample_step_sec : float, optional
        Time resolution (seconds) for downsampling the smoothed power traces before infraslow analysis (default is 0.5).

    Returns
    -------
    results : dict
        Dictionary with a key for each band in `band_defs`. Each value is a dictionary with:
        
        - 'infraslow_freqs' : ndarray
            Frequencies used in the second-stage (infraslow) wavelet transform.
        - 'psd' : ndarray
            Average power spectral density across time points for the band-specific infraslow transform.
        - 'smoothed_power_trace' : ndarray
            The downsampled, smoothed, and normalized power time series used for the infraslow transform.
        - 'duration' : int
            Number of original EEG samples (i.e., length of `eeg_signal`).
    """
    
    dt = 1 / fs
    n = len(eeg_signal)

    # Time-frequency decomposition (0.5–24 Hz)
    freqs_high = arange(0.5, 24.1, 0.2)
    scales_high = pywt.central_frequency(wavelet) / (freqs_high * dt)
    coeffs, freqs_used = pywt.cwt(eeg_signal, scales_high, wavelet, sampling_period=dt)
    power = abs(coeffs) ** 2

    results = {}
    
    # Calculate PSD on entire segment, for normalisation.
    f_fft, pxx = welch(baseline_data, fs=fs, nperseg=fs*4)
    

    for band_name, (fmin, fmax) in band_defs.items():
        
        # Band power time series
        band_idx = (freqs_used >= fmin) & (freqs_used <= fmax)
        band_power = power[band_idx].mean(axis=0)

        # Moving average smoothing
        win_samples = int(smoothing_window_sec * fs)
        smooth = Series(band_power).rolling(win_samples, center=True, min_periods=1).mean().to_numpy()

        # Normalize to mean power
        idx_band = (f_fft >= fmin) & (f_fft <= fmax)
        band_baseline = nanmean(pxx[idx_band])
        normalized = smooth / band_baseline

        # Downsample for infraslow resolution
        step = int(downsample_step_sec * fs)
        smooth_ds = normalized[::step]
        dt_ds = downsample_step_sec

        # Infraslow wavelet transform (e.g. 0.001–0.12 Hz)
        freqs_low = arange(infraslow_freq_range[0], infraslow_freq_range[1], infraslow_step)
        scales_low = pywt.central_frequency(wavelet) / (freqs_low * dt_ds)
        coeffs_low, _ = pywt.cwt(smooth_ds, scales_low, wavelet, sampling_period=dt_ds)
        power_low = abs(coeffs_low) ** 2
        mean_power = power_low.mean(axis=1)

        results[band_name] = {
                              'infraslow_freqs': freqs_low,
                              'psd': mean_power,
                              'smoothed_power_trace': smooth_ds,
                              'duration': n
                             }

    return results


def gaussian(x, a, mu, sigma):
    return a * exp(-((x - mu) ** 2) / (2 * sigma ** 2))

def fit_gaussian(freqs_low, psd_vals):
    """
    Fit a Gaussian curve to a 1D power spectrum and compute peak characteristics.

    Parameters
    ----------
    freqs_low : array_like
        1D array of frequency values (e.g., in the infraslow range).
    psd_vals : array_like
        1D array of power spectral density (PSD) values corresponding to `freqs_low`.

    Returns
    -------
    peak_freq : float
        Estimated peak frequency (mu) of the fitted Gaussian [Hz].
    sigma : float
        Standard deviation (spread) of the fitted Gaussian [Hz].
    avg_peak_power : float
        Average PSD value in the frequency window centered at the peak frequency
        and spanning ±0.5 * sigma, computed via linear interpolation.
        Returns NaN for all outputs if the fit fails.
    """
    
    try:
        # Example: constrain mu to [0.01, 0.05], and sigma > 0
        bounds = ([0, 0.01, 1e-4], [inf, 0.05, 0.1])  # [a, mu, sigma]
        
        # Fit Gaussian
        popt, _ = curve_fit(gaussian, freqs_low, psd_vals, 
                            p0 = [1.0, 0.02, 0.01],
                            bounds = bounds)
        a, peak_freq, sigma = popt
        peak_band = (peak_freq - 0.5 * sigma, peak_freq + 0.5 * sigma)
        
        # Clip window to valid frequency range
        peak_band = (
                    max((peak_band[0], freqs_low.min())),
                    min((peak_band[1], freqs_low.max()))
                    )
        
        # Interpolate PSD for smoother power estimate in peak window
        interp_fn = interp1d(freqs_low, psd_vals, 
                             kind='linear', 
                             bounds_error=False, 
                             fill_value='extrapolate')
        fine_freqs = linspace(*peak_band, 100)
        fine_power = interp_fn(fine_freqs)
        
        avg_peak_power = nanmean(fine_power)
        
    except RuntimeError:
        peak_freq, sigma, avg_peak_power = nan, nan, nan

    return peak_freq, sigma, avg_peak_power
    
            