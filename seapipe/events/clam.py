#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  7 15:47:20 2025

@author: ncro8394
"""

from numpy import (angle, arange, array, average, ceil, concatenate, degrees, diff, exp,
                   hamming, hanning, histogram,  inf, isfinite, isnan, linspace, max, min, nan, nanargmax, nan_to_num,
                   nanargmin, nanmax, nanmean, nanmedian, nanmin, nanpercentile, nansum, nanstd, pi, random, sort, vstack, where, zeros, zeros_like)
from numpy.fft import rfft, rfftfreq
from pandas import DataFrame, Series
from collections import Counter, defaultdict
from copy import deepcopy
from os import listdir, mkdir, path
import pywt
import shutil
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter1d
from scipy.signal import (butter, firwin, filtfilt, find_peaks, detrend, 
                          hilbert, morlet2, periodogram, resample, sosfiltfilt, welch) 
import builtins
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
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
                             concat_stage = False,
                             spectral_method = 'welch',
                             min_total_nrem_sec = None,
                             min_bouts_psd = None,
                             low_snr_percentile = None,
                             plot_fit = False,
                             logger = create_logger('Event clustering')):
            
            
            
            ### 0.a Set up logging
            tracking = self.tracking
            flag = 0

            spectral_method = str(spectral_method).lower()
            if spectral_method not in ('welch', 'morlet'):
                logger.warning(f"Unknown spectral_method '{spectral_method}', defaulting to 'welch'.")
                spectral_method = 'welch'
            if low_snr_percentile is not None:
                try:
                    low_snr_percentile = float(low_snr_percentile)
                except Exception:
                    logger.warning('low_snr_percentile is not numeric; disabling SNR filtering.')
                    low_snr_percentile = None
                if low_snr_percentile is not None and not (0 < low_snr_percentile < 100):
                    logger.warning('low_snr_percentile must be between 0 and 100; disabling SNR filtering.')
                    low_snr_percentile = None
            
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
                            
                        if isinstance(self.stage, (list, tuple)):
                            stage_all = list(self.stage)
                        else:
                            stage_all = [self.stage]

                        if concat_stage:
                            stage_runs = [{
                                'stage_list': stage_all,
                                'stagename': '-'.join(stage_all),
                                'cat': (0, 1, 0, 0)
                            }]
                        else:
                            stage_runs = [{
                                'stage_list': [st],
                                'stagename': st,
                                'cat': (0, 0, 0, 0)
                            } for st in stage_all]

                        for stage_cfg in stage_runs:
                            stage_list = stage_cfg['stage_list']
                            stagename = stage_cfg['stagename']
                            cat = stage_cfg['cat']

                            # ---- Analys event clustering ----
                            logger.debug(f'Analysing cluster metrics for {evt_type}, {ch}')
                            stats = cluster_metrics(annot, ch, evt_type, 
                                                    stage_list, grp_name)

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

                            # ---- Analyse low-frequency fluctuations (human-style settings) ----
                            logger.debug('Analysing low frequency fluctuations for '
                                         f"{', '.join([band for band in freq_bands])} "
                                         f'in channel: {fnamechan}: {chanset[ch]}')
                            # Human method: Morlet CWT 0.5–24 Hz, 0.2 Hz step, 4-cycle; 4-s smoothing; 0.5 s sampling
                            window_sec = 5
                            overlap_sec = 4.5
                            step_sec = window_sec - overlap_sec  # 0.5 s steps
                            baseline_limit_sec = 100 * 60  # first 100 min for normalization
                            min_bout_sec = 120  # human method minimum bout length
                            # Human methods: include bouts >=120 s; Welch window capped at 300 s
                            min_psd_bout_sec = 120
                            infraslow_step_hz = 0.001
                            # Morlet method targets 0.001–0.12 Hz; Welch can compute up to 0.5 Hz
                            infraslow_max_hz = 0.12 if spectral_method == 'morlet' else 0.5

                            per_seg_series = []
                            baseline_vals = {band: [] for band in freq_bands}
                            freqs_high = arange(0.5, 24.1, 0.2)

                            epochs = annot.get_epochs()
                            sleep_onset = None
                            n1_set = ['NREM1', 'N1', 'S1']
                            n1n2_set = ['NREM1', 'N1', 'S1', 'NREM2', 'N2', 'S2']
                            # Sleep onset: first NREM1 followed by at least 2 consecutive NREM1/NREM2 epochs.
                            for i in range(len(epochs) - 2):
                                if (epochs[i]['stage'] in n1_set and
                                        epochs[i + 1]['stage'] in n1n2_set and
                                        epochs[i + 2]['stage'] in n1n2_set):
                                    sleep_onset = epochs[i]['start']
                                    break
                            if sleep_onset is None:
                                # Fallback to first NREM2 if NREM1 transition is absent.
                                for ep in epochs:
                                    if ep['stage'] in ['NREM2', 'N2', 'S2']:
                                        sleep_onset = ep['start']
                                        break
                            if sleep_onset is None and epochs:
                                sleep_onset = epochs[0]['start']
                            if sleep_onset is None:
                                logger.warning('Sleep onset not found; skipping...')
                                flag += 1
                                continue

                            analysis_start = sleep_onset
                            analysis_end = sleep_onset + (210 * 60)

                            allowed_stages = set()
                            for st in stage_list:
                                if st in ['NREM1', 'N1', 'S1']:
                                    allowed_stages.update(['NREM1', 'N1', 'S1'])
                                elif st in ['NREM2', 'N2', 'S2']:
                                    allowed_stages.update(['NREM2', 'N2', 'S2'])
                                elif st in ['NREM3', 'N3', 'S3', 'S4']:
                                    allowed_stages.update(['NREM3', 'N3', 'S3', 'S4'])
                                else:
                                    allowed_stages.add(st)

                            epoch_stages = {ep['stage'] for ep in epochs}
                            stage_fetch = [st for st in allowed_stages if st in epoch_stages]
                            if not stage_fetch:
                                stage_fetch = [st for st in stage_list if st in epoch_stages]

                            segs = fetch(dset, annot, cat = cat,
                                         stage = stage_fetch,
                                         reject_artf = True)
                            segs.read_data(chan = [ch], ref_chan = chanset[ch])

                            seg_infos = []
                            for seg in segs:
                                signal_data = seg['data'].data[0][0]
                                sampling_rate = seg['data'].s_freq
                                if sampling_rate is None or sampling_rate <= 0:
                                    logger.warning('Sampling rate not found for segment; skipping.')
                                    continue
                                time_axis = seg['data'].axis['time'][0]
                                seg_infos.append({
                                    'start': time_axis.min(),
                                    'end': time_axis.max(),
                                    'time': time_axis,
                                    'data': signal_data,
                                    'fs': sampling_rate
                                })

                            nrem_epochs = []
                            for ep in epochs:
                                if ep['stage'] not in allowed_stages:
                                    continue
                                if ep['start'] < analysis_start or ep['end'] > analysis_end:
                                    continue
                                if any(ep['start'] >= s['start'] and ep['end'] <= s['end'] for s in seg_infos):
                                    nrem_epochs.append(ep)

                        if len(nrem_epochs) == 0:
                            logger.warning(f'No valid NREM epochs found for {sub}. {ses}. Skipping...')
                            flag += 1
                            continue

                        nrem_epochs.sort(key=lambda x: x['start'])
                        bouts = []
                        bout_start = nrem_epochs[0]['start']
                        bout_end = nrem_epochs[0]['end']
                        for ep in nrem_epochs[1:]:
                            if abs(ep['start'] - bout_end) < 1e-6:
                                bout_end = ep['end']
                            else:
                                bouts.append((bout_start, bout_end))
                                bout_start = ep['start']
                                bout_end = ep['end']
                        bouts.append((bout_start, bout_end))

                        for bout_start, bout_end in bouts:
                            if (bout_end - bout_start) < min_bout_sec:
                                continue
                            seg_match = next((s for s in seg_infos if s['start'] <= bout_start and s['end'] >= bout_end), None)
                            if seg_match is None:
                                continue
                            fs = seg_match['fs']
                            time_axis = seg_match['time']
                            mask = (time_axis >= bout_start) & (time_axis <= bout_end)
                            if not mask.any():
                                continue
                            signal_data = seg_match['data'][mask]

                            dt = 1 / fs
                            scales_high = pywt.central_frequency('cmor4.0-1.0') / (freqs_high * dt)
                            coeffs, freqs_used = pywt.cwt(signal_data, scales_high, 'cmor4.0-1.0', sampling_period=dt)
                            power = abs(coeffs) ** 2

                            band_power_raw = {}
                            for band, (fmin, fmax) in freq_bands.items():
                                band_idx = (freqs_used >= fmin) & (freqs_used <= fmax)
                                band_power = power[band_idx].mean(axis=0)
                                win_samples = int(4 * fs)
                                band_power = Series(band_power).rolling(win_samples, center=True, min_periods=1).mean().to_numpy()
                                step = int(step_sec * fs)
                                band_power = band_power[::step]
                                band_power_raw[band] = band_power
                                if baseline_limit_sec > 0:
                                    n_take = int(builtins.min(len(band_power), builtins.max(0, baseline_limit_sec / step_sec)))
                                    if n_take > 0:
                                        baseline_vals[band].extend(band_power[:n_take].tolist())
                            times = bout_start + arange(len(list(band_power_raw.values())[0])) * step_sec
                            per_seg_series.append({
                                'start_time': times[0],
                                'end_time': times[-1],
                                'times': times,
                                'band_power_raw': band_power_raw,
                                'sampling_rate_band': 1 / step_sec
                            })
                            baseline_limit_sec = builtins.max(0, baseline_limit_sec - len(times) * step_sec)

                        if len(per_seg_series) == 0:
                            logger.warning(f'No valid NREM bouts found for {sub}. {ses}. Skipping...')
                            flag += 1
                            continue

                        # Minimum total NREM coverage (toggle)
                        if min_total_nrem_sec is not None:
                            total_nrem_sec = sum(len(seg['times']) * step_sec for seg in per_seg_series)
                            if total_nrem_sec < min_total_nrem_sec:
                                logger.warning(f'Total NREM duration {total_nrem_sec:.1f}s < '
                                               f'{min_total_nrem_sec}s; skipping {sub}, {ses}, {fnamechan}, {stagename}.')
                                flag += 1
                                continue

                        baseline_band_power = {}
                        for band in freq_bands:
                            vals = array(baseline_vals[band])
                            baseline_band_power[band] = nanmean(vals) if len(vals) > 0 else nan

                        for seg_series in per_seg_series:
                            for band in freq_bands:
                                denom = baseline_band_power.get(band, nan)
                                seg_series[f'{band}_norm'] = seg_series['band_power_raw'][band] / denom

                        trace_segments = {band: [] for band in freq_bands}
                        psd_accum = {band: [] for band in freq_bands}
                        psd_weights = {band: [] for band in freq_bands}

                        # Precompute eligible bouts by duration
                        eligible_idx = [i for i, seg in enumerate(per_seg_series)
                                        if len(seg['times']) * step_sec >= min_psd_bout_sec]

                        # Minimum number of PSD-eligible bouts (toggle)
                        if min_bouts_psd is not None:
                            try:
                                min_bouts_psd_val = int(min_bouts_psd)
                            except Exception:
                                min_bouts_psd_val = None
                            if min_bouts_psd_val is not None and len(eligible_idx) < min_bouts_psd_val:
                                logger.warning(f'Only {len(eligible_idx)} bouts >= {min_psd_bout_sec}s '
                                               f'(min required {min_bouts_psd_val}); skipping {sub}, {ses}, '
                                               f'{fnamechan}, {stagename}.')
                                flag += 1
                                continue

                        # Low-SNR bout filter thresholds per band (toggle)
                        snr_thresh = {}
                        if low_snr_percentile is not None and len(eligible_idx) > 0:
                            for band in freq_bands:
                                vals = []
                                for i in eligible_idx:
                                    series_raw = per_seg_series[i]['band_power_raw'][band]
                                    vals.append(nanmean(series_raw))
                                if len(vals) > 0:
                                    snr_thresh[band] = nanpercentile(array(vals, dtype=float), low_snr_percentile)

                        for seg_series in per_seg_series:
                            duration_sec = len(seg_series['times']) * step_sec
                            if duration_sec < min_psd_bout_sec:
                                continue
                            for band in freq_bands:
                                if band in snr_thresh:
                                    band_mean = nanmean(seg_series['band_power_raw'][band])
                                    if isfinite(band_mean) and band_mean < snr_thresh[band]:
                                        continue
                                series = nan_to_num(array(seg_series[f'{band}_norm'], dtype=float),
                                                     nan=nanmean(seg_series[f'{band}_norm']))
                                fs_band = 1 / step_sec

                                if spectral_method == 'morlet':
                                    # Morlet wavelet on band-power timecourse (0.001–0.12 Hz, 0.001 step)
                                    freqs_isf = arange(0.001, 0.120 + infraslow_step_hz, infraslow_step_hz)
                                    dt = step_sec
                                    scales = pywt.central_frequency('cmor4.0-1.0') / (freqs_isf * dt)
                                    coeffs, freqs_used = pywt.cwt(series, scales, 'cmor4.0-1.0', sampling_period=dt)
                                    power = abs(coeffs) ** 2
                                    power_mean = nanmean(power, axis=1)
                                    freqs_low = arange(0, infraslow_max_hz + infraslow_step_hz, infraslow_step_hz)
                                    interp_fn = interp1d(freqs_used, power_mean, bounds_error=False, fill_value=nan)
                                    power_interp = interp_fn(freqs_low)
                                else:
                                    # High-pass to suppress slow drift (MATLAB-style FIR) or detrend if too short
                                    hp_cut = 0.005
                                    transition_hz = 0.01
                                    forder = int(ceil(3.3 / (transition_hz / fs_band)))
                                    numtaps = forder + 1
                                    if numtaps % 2 == 0:
                                        numtaps += 1
                                    if len(series) <= numtaps * 3:
                                        series = detrend(series)
                                    else:
                                        b = firwin(numtaps, hp_cut, pass_zero=False, fs=fs_band, window='hamming')
                                        series = filtfilt(b, [1.0], series)
                                    centered = series
                                    
                                    # PWELCH on band-power timecourse
                                    nperseg = int(builtins.min(len(centered), int(300 / step_sec)))
                                    if nperseg < 4:
                                        continue
                                    noverlap = int(nperseg * 0.5)
                                    freqs_pwel, power_pwel = welch(centered, fs=fs_band,
                                                                   window='hann',
                                                                   nperseg=nperseg,
                                                                   noverlap=noverlap,
                                                                   nfft=nperseg,
                                                                   detrend=False,
                                                                   scaling='density')
                                    mask = freqs_pwel <= infraslow_max_hz
                                    freqs_pwel = freqs_pwel[mask]
                                    power_pwel = power_pwel[mask]
                                    
                                    # Interpolate to common grid (0..infraslow_max_hz, fixed spacing)
                                    freqs_low = arange(0, infraslow_max_hz + infraslow_step_hz, infraslow_step_hz)
                                    interp_fn = interp1d(freqs_pwel, power_pwel, bounds_error=False, fill_value='extrapolate')
                                    power_interp = interp_fn(freqs_low)

                                psd_accum[band].append(power_interp)
                                psd_weights[band].append(duration_sec)

                        # Bout-length summary (all bouts + PSD-eligible bouts)
                        all_bout_lengths = array([len(seg['times']) * step_sec for seg in per_seg_series], dtype=float)
                        psd_bout_lengths = all_bout_lengths[all_bout_lengths >= min_psd_bout_sec]

                        stats = {
                            'n_bouts_total': int(len(all_bout_lengths)),
                            'mean_bout_len_s': float(nanmean(all_bout_lengths)) if len(all_bout_lengths) > 0 else nan,
                            'std_bout_len_s': float(nanstd(all_bout_lengths)) if len(all_bout_lengths) > 0 else nan,
                            'n_bouts_psd': int(len(psd_bout_lengths)),
                            'mean_bout_len_psd_s': float(nanmean(psd_bout_lengths)) if len(psd_bout_lengths) > 0 else nan,
                            'std_bout_len_psd_s': float(nanstd(psd_bout_lengths)) if len(psd_bout_lengths) > 0 else nan,
                        }
                        for band in freq_bands:
                            arrays = array(psd_accum[band])
                            if len(arrays) == 0:
                                logger.warning(f'No valid bouts >= {min_psd_bout_sec}s for {band}; skipping.')
                                continue
                            weights = array(psd_weights[band], dtype=float)
                            if len(weights) != len(arrays) or nansum(weights) == 0:
                                weighted_avg = average(arrays, axis=0)
                            else:
                                weighted_avg = average(arrays, axis=0, weights=weights)
                            
                            # Normalize by mean power in 0.0075–0.1 Hz (Matlab style)
                            norm_mask = (freqs_low >= 0.0075) & (freqs_low <= 0.1)
                            mean_psd = weighted_avg / nanmean(weighted_avg[norm_mask])
                            
                            # No extra smoothing beyond PWELCH interpolation
                            fit_mask = (freqs_low >= 0.0075) & (freqs_low <= 0.1)
                            if not fit_mask.any():
                                logger.warning(f'No frequencies in fit window for {band}; skipping.')
                                continue
                            peak_freq, sigma_val, avg_peak_power, fit_curve, fit_params = fit_gaussian(
                                freqs_low[fit_mask], mean_psd[fit_mask], num_components=2,
                                full_freqs=freqs_low, return_params=True)

                            logger.debug(f"Chan: {fnamechan}: {chanset[ch]} - "
                                         f"{band} infraslow peak: {round(peak_freq,3)}; "
                                         f"Peak power: {round(avg_peak_power,3)}")

                            stats[f'{band}_peak_freq'] = round(peak_freq, 5)
                            stats[f'{band}_avg_peak_power'] = round(avg_peak_power, 3)
                            stats[f'{band}_sigma'] = round(sigma_val, 5)
                            stats[f'{band}_fit_adjR2'] = round(fit_params.get('adjr2', nan), 3)
                            stats[f'{band}_fit_model'] = fit_params.get('model_used', nan)
                            stats[f'{band}_fit_poor'] = fit_params.get('poor_fit', False)
                            stats[f'{band}_fit_peak_source'] = fit_params.get('peak_source', 'fit')

                            # Phase estimation using band-pass around peak ±1*SD
                            if isnan(peak_freq) or isnan(sigma_val):
                                logger.warning(f'Gaussian fit failed for {band}; skipping phase and trace outputs.')
                                continue

                            passband = (builtins.max(peak_freq - sigma_val, 0.0005), builtins.min(peak_freq + sigma_val, infraslow_max_hz))
                            if passband[0] <= 0 or passband[0] >= passband[1]:
                                logger.warning(f'Invalid passband for {band}: {passband}; skipping.')
                                continue
                            sos = butter(4, [passband[0], passband[1]], btype='bandpass',
                                         fs=1 / step_sec, output='sos')
                            band_phases = []
                            for seg_series in per_seg_series:
                                series = nan_to_num(array(seg_series[f'{band}_norm'], dtype=float),
                                                     nan=nanmean(seg_series[f'{band}_norm']))
                                if len(series) * step_sec < min_bout_sec:
                                    continue
                                filtered = sosfiltfilt(sos, series)
                                analytic_signal = hilbert(filtered)
                                inst_phase = angle(analytic_signal)
                                start_seg = seg_series['start_time']
                                end_seg = seg_series['end_time']
                                events = [x for x in annot.get_events(evt_type, chan=f'{ch} ({grp_name})')
                                          if x['start'] > start_seg and x['end'] < end_seg]
                                evts_mid = [((x['start'] + x['end']) / 2) for x in events]
                                times_vec = seg_series['times']
                                evts_idx = [nanargmin(abs(times_vec - mid)) for mid in evts_mid] if len(evts_mid) > 0 else []
                                if len(evts_idx) > 0:
                                    band_phases.extend(inst_phase[evts_idx].tolist())
                                
                                # Save filtered trace (aligned to original times)
                                trace_segments[band].append(vstack((times_vec, filtered)))
                            
                            # Histogram of events
                            bins = linspace(-pi, pi, 19)
                            hist, bin_edges = histogram(band_phases, bins=bins)
                            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
                            hist_df = DataFrame(hist).transpose()
                            hist_df.columns = [int(round(degrees(x))) for x in bin_centers]
                            hist_df.to_csv(f'{outpath}/{sub}_{ses}_{fnamechan}_{stagename}_{evt_type}_{band}_coupling.csv')

                            # PSD of LowFreq-Fluctuation
                            psd_df = DataFrame({
                                'freq_hz': freqs_low,
                                'welch_psd': mean_psd,
                                'fit_psd': fit_curve
                            })
                            psd_df.to_csv(f'{outpath}/{sub}_{ses}_{fnamechan}_{stagename}_{band}_fluctuations_psd.csv')

                            if plot_fit:
                                plot_mu_bounds = fit_params.get('mu_bounds') if isinstance(fit_params, dict) else (0.0075, 0.04)
                                plot_out = f'{outpath}/{sub}_{ses}_{fnamechan}_{stagename}_{band}_fluctuations_fit.png'
                                plot_psd_fit(freqs_low, mean_psd, fit_curve, fit_params,
                                             plot_mu_bounds, plot_out,
                                             title=f'{sub} {ses} {fnamechan} {band}')

                            # Filtered trace
                            if trace_segments[band]:
                                combined_trace = concatenate(trace_segments[band], axis=1)
                                ts_df = DataFrame([combined_trace[1]],
                                                  columns=combined_trace[0])
                                ts_df.to_csv(f'{outpath}/{sub}_{ses}_{fnamechan}_{stagename}_{band}_fluctuations_timeseries.csv')

                        # Save stats to file
                        if len(stats) == 0:
                            logger.warning(f'No valid infraslow stats computed for {sub}. {ses}. Skipping...')
                            flag += 1
                        else:
                            df = DataFrame([stats])
                            df.insert(0, 'sub', sub)
                            df.insert(1, 'ses', ses)
                            bandnames = '-'.join([band for band in freq_bands])
                            df.to_csv(f'{outpath}/{sub}_{ses}_{fnamechan}_{stagename}_{bandnames}_fluctuation_stats.csv')

            ### 3. Check completion status and print
            if flag == 0:
                logger.debug('Event clustering analysis finished without error.')  
            else:
                logger.warning(f'Event clustering analysis finished with {flag} WARNINGS. See log for details.')
            
            #self.tracking = tracking   ## TO UPDATE - FIX TRACKING
            
            return 


        def plot_infraslow_panels(self, sub = None, ses = None,
                                  subs = None, sessions = None,
                                  stage = None, channel = None,
                                  xml_dir = None, out_dir = None, seed = 1,
                                  target_fs = 100,
                                  logger = create_logger('Infraslow panels')):
            """
            Create diagnostic panels (A-C) for the infraslow workflow on a single subject/session.

            Panel A: Timeline, hypnogram (with 210-min window shaded), NREM bouts (with 100-min
                     normalization window shaded), and bouts >=120s.
            Panel B: EEG signal from one randomly selected bout >=120s.
            Panel C: Normalized sigma band-power time series (raw) and its 4-s smoothed version.
            """

            # ---- Loop over dataset if sub/ses not specified ----
            if sub is None or ses is None:
                return _plot_infraslow_panels_all(self, subs=subs, sessions=sessions,
                                                  stage=stage, xml_dir=xml_dir,
                                                  out_dir=out_dir, seed=seed,
                                                  target_fs=target_fs, logger=logger)

            # ---- Resolve staging XML ----
            xml_root = xml_dir if xml_dir is not None else self.xml_dir
            if not path.isabs(xml_root):
                xml_root = path.join(self.rootpath, xml_root)
            xml_path = path.join(xml_root, sub, ses, f'{sub}_{ses}_eeg.xml')
            if not path.exists(xml_path):
                xml_dir_alt = path.join(xml_root, sub, ses)
                if path.exists(xml_dir_alt):
                    xmls = [f for f in listdir(xml_dir_alt) if f.endswith('.xml')]
                    if len(xmls) > 0:
                        xml_path = path.join(xml_dir_alt, xmls[0])
            if not path.exists(xml_path):
                logger.warning(f'Staging XML not found for {sub} {ses}: {xml_path}')
                return

            # ---- Parse dataset path and duration ----
            eeg_path = None
            last_second = None
            try:
                tree = ET.parse(xml_path)
                root = tree.getroot()
                ds = root.find('dataset')
                if ds is not None:
                    last_txt = ds.findtext('last_second')
                    if last_txt:
                        last_second = float(last_txt)
                    eeg_path = ds.findtext('path') or ds.findtext('filename')
            except Exception as exc:
                logger.warning(f'Failed to parse dataset info from XML: {exc}')

            if eeg_path is None or not path.exists(eeg_path):
                rec_root = self.rec_dir
                if rec_root is not None:
                    if not path.isabs(rec_root):
                        rec_root = path.join(self.rootpath, rec_root)
                    guess = path.join(rec_root, sub, ses, 'eeg', f'{sub}_{ses}_eeg.edf')
                    if path.exists(guess):
                        eeg_path = guess

            if eeg_path is None or not path.exists(eeg_path):
                logger.warning(f'EEG file not found for {sub} {ses}: {eeg_path}')
                return

            # ---- Load dataset and annotations ----
            dset = Dataset(eeg_path)
            annot = Annotations(xml_path)
            epochs = annot.get_epochs()
            if not epochs:
                logger.warning(f'No epochs found in XML for {sub} {ses}')
                return

            # ---- Channel selection ----
            chan_names = dset.header.get('chan_name', [])
            def pick_channel(names, desired):
                if desired:
                    for name in names:
                        if desired == name or desired in name:
                            return name
                for pref in ['C3', 'C4', 'Cz', 'Fz', 'F3', 'F4', 'Pz', 'P3', 'P4']:
                    for name in names:
                        if pref in name:
                            return name
                for name in names:
                    if all(tag not in name.upper() for tag in ['EOG', 'EMG', 'ECG', 'EKG']):
                        return name
                return names[0] if len(names) > 0 else None

            picked_chan = pick_channel(chan_names, channel)
            if picked_chan is None:
                logger.warning('No EEG channel found to plot.')
                return

            fnamechan = channel if channel is not None else picked_chan
            fnamechan = fnamechan.replace(' ', '').replace(':', '-')

            # ---- Sleep onset (first NREM1 followed by >=2 consecutive NREM1/NREM2 epochs) ----
            sleep_onset = None
            n1_set = ['NREM1', 'N1', 'S1']
            n1n2_set = ['NREM1', 'N1', 'S1', 'NREM2', 'N2', 'S2']
            for i in range(len(epochs) - 2):
                if (epochs[i]['stage'] in n1_set and
                        epochs[i + 1]['stage'] in n1n2_set and
                        epochs[i + 2]['stage'] in n1n2_set):
                    sleep_onset = epochs[i]['start']
                    break
            if sleep_onset is None:
                for ep in epochs:
                    if ep['stage'] in ['NREM2', 'N2', 'S2']:
                        sleep_onset = ep['start']
                        break
            if sleep_onset is None and epochs:
                sleep_onset = epochs[0]['start']
            if sleep_onset is None:
                logger.warning('Sleep onset not found; skipping...')
                return

            analysis_start = sleep_onset
            analysis_end = sleep_onset + (210 * 60)

            # ---- Allowed stages ----
            stage_list = stage if stage is not None else self.stage
            if isinstance(stage_list, str):
                stage_list = stage_list.split('-')
            stage_list = stage_list or ['NREM2', 'NREM3']

            allowed_stages = set()
            for st in stage_list:
                if st in ['NREM1', 'N1', 'S1']:
                    allowed_stages.update(['NREM1', 'N1', 'S1'])
                elif st in ['NREM2', 'N2', 'S2']:
                    allowed_stages.update(['NREM2', 'N2', 'S2'])
                elif st in ['NREM3', 'N3', 'S3', 'S4']:
                    allowed_stages.update(['NREM3', 'N3', 'S3', 'S4'])
                else:
                    allowed_stages.add(st)

            epoch_stages = {ep['stage'] for ep in epochs}
            stage_fetch = [st for st in allowed_stages if st in epoch_stages]
            if not stage_fetch:
                stage_fetch = [st for st in stage_list if st in epoch_stages]

            # ---- Fetch artifact-free segments ----
            segs = fetch(dset, annot, cat = (0, 0, 0, 0),
                         stage = stage_fetch, reject_artf = True)
            segs.read_data(chan = [picked_chan])

            seg_infos = []
            for seg in segs:
                signal_data = seg['data'].data[0][0]
                sampling_rate = seg['data'].s_freq
                time_axis = seg['data'].axis['time'][0]
                seg_infos.append({
                    'start': time_axis.min(),
                    'end': time_axis.max(),
                    'time': time_axis,
                    'data': signal_data,
                    'fs': sampling_rate
                })

            # ---- NREM epochs within analysis window and artifact-free ----
            nrem_epochs = []
            for ep in epochs:
                if ep['stage'] not in allowed_stages:
                    continue
                if ep['start'] < analysis_start or ep['end'] > analysis_end:
                    continue
                if any(ep['start'] >= s['start'] and ep['end'] <= s['end'] for s in seg_infos):
                    nrem_epochs.append(ep)

            if len(nrem_epochs) == 0:
                logger.warning(f'No valid NREM epochs found for {sub} {ses}')
                return

            nrem_epochs.sort(key=lambda x: x['start'])

            # ---- Build contiguous NREM bouts ----
            bouts = []
            bout_start = nrem_epochs[0]['start']
            bout_end = nrem_epochs[0]['end']
            for ep in nrem_epochs[1:]:
                if abs(ep['start'] - bout_end) < 1e-6:
                    bout_end = ep['end']
                else:
                    bouts.append((bout_start, bout_end))
                    bout_start = ep['start']
                    bout_end = ep['end']
            bouts.append((bout_start, bout_end))

            bouts_long = [b for b in bouts if (b[1] - b[0]) >= 120]

            # ---- Normalization interval: first 100 min of NREM ----
            norm_intervals = []
            norm_remaining = 100 * 60
            for b_start, b_end in bouts:
                if norm_remaining <= 0:
                    break
                dur = b_end - b_start
                if dur <= norm_remaining:
                    norm_intervals.append((b_start, b_end))
                    norm_remaining -= dur
                else:
                    norm_intervals.append((b_start, b_start + norm_remaining))
                    norm_remaining = 0

            # ---- Panel A: timeline + hypnogram + NREM bouts ----
            if last_second is None:
                last_second = dset.header['n_samples'] / dset.header['s_freq']
            last_min = last_second / 60

            fig_a, axes = plt.subplots(4, 1, figsize=(12, 8), sharex=True,
                                       gridspec_kw={'height_ratios': [0.3, 1.2, 0.6, 0.6]})
            ax0, ax1, ax2, ax3 = axes

            # A1: timeline
            ax0.plot([0, last_min], [0, 0], color='k', linewidth=1)
            ax0.set_yticks([])
            ax0.set_ylabel('Time')

            # A2: hypnogram
            stage_map = {
                'NREM3': 0, 'N3': 0, 'S3': 0, 'S4': 0,
                'NREM2': 1, 'N2': 1, 'S2': 1,
                'NREM1': 2, 'N1': 2, 'S1': 2,
                'REM': 3, 'R': 3,
                'Wake': 4, 'W': 4
            }
            t_step = [epochs[0]['start'] / 60]
            v_step = []
            for ep in epochs:
                t_step.append(ep['end'] / 60)
                v_step.append(stage_map.get(ep['stage'], nan))
            if len(v_step) > 0:
                v_step.append(v_step[-1])
                ax1.step(t_step, v_step, where='post', color='k', linewidth=0.8)
            ax1.set_yticks([0, 1, 2, 3, 4])
            ax1.set_yticklabels(['N3', 'N2', 'N1', 'REM', 'Wake'])
            ax1.set_ylabel('Stage')
            ax1.axvspan(analysis_start / 60, analysis_end / 60, color='0.9')

            # A3: NREM bouts + normalization interval
            for b_start, b_end in bouts:
                ax2.hlines(1, b_start / 60, b_end / 60, color='k', linewidth=6)
            for n_start, n_end in norm_intervals:
                ax2.axvspan(n_start / 60, n_end / 60, color='0.8')
            ax2.set_yticks([])
            ax2.set_ylabel('NREM')

            # A4: bouts >=120s
            for b_start, b_end in bouts_long:
                ax3.hlines(1, b_start / 60, b_end / 60, color='k', linewidth=6)
            ax3.set_yticks([])
            ax3.set_ylabel('>=120s')
            ax3.set_xlabel('Minutes from start')

            # ---- Save panel A early (panel B/C optional) ----
            out_root = out_dir if out_dir is not None else self.out_dir
            if not path.isabs(out_root):
                out_root = path.join(self.rootpath, out_root)
            if not path.exists(out_root):
                mkdir(out_root)
            out_sub = path.join(out_root, sub)
            if not path.exists(out_sub):
                mkdir(out_sub)
            out_ses = path.join(out_sub, ses)
            if not path.exists(out_ses):
                mkdir(out_ses)
            out_path = path.join(out_ses, 'figures')
            if not path.exists(out_path):
                mkdir(out_path)

            fig_a.tight_layout()
            fig_a.savefig(path.join(out_path, f'{sub}_{ses}_{fnamechan}_panel_A.png'), dpi=150)
            plt.close(fig_a)

            # ---- Panel B/C: select a bout ----
            if len(bouts_long) == 0:
                logger.warning(f'No bouts >=120s for {sub} {ses}; panel B/C skipped.')
                return

            random.seed(seed)
            sel_idx = int(random.randint(len(bouts_long)))
            sel_start, sel_end = bouts_long[sel_idx]

            def extract_bout_data(bout_start, bout_end):
                seg_match = next((s for s in seg_infos
                                  if s['start'] <= bout_start and s['end'] >= bout_end), None)
                if seg_match is None:
                    return None, None
                time_axis = seg_match['time']
                mask = (time_axis >= bout_start) & (time_axis <= bout_end)
                if not mask.any():
                    return None, None
                data = seg_match['data'][mask]
                fs = seg_match['fs']
                return data, fs

            eeg_bout, fs = extract_bout_data(sel_start, sel_end)
            if eeg_bout is None:
                logger.warning('Selected bout data not found; panel B/C skipped.')
                return

            # Downsample to target_fs (method uses 100 Hz)
            if target_fs is not None and fs > target_fs:
                n_samp = int(len(eeg_bout) * (target_fs / fs))
                eeg_bout = resample(eeg_bout, n_samp)
                fs = target_fs

            # ---- Panel B: EEG signal ----
            t_bout = arange(len(eeg_bout)) / fs
            fig_b, axb = plt.subplots(1, 1, figsize=(12, 3))
            axb.plot(t_bout, eeg_bout, color='k', linewidth=0.8)
            axb.set_xlabel('Time (s)')
            axb.set_ylabel('EEG (a.u.)')
            axb.set_title(f'{sub} {ses} {picked_chan} - Bout EEG')

            # ---- Panel C: normalized sigma power (raw and smoothed) ----
            freqs_high = arange(0.5, 24.1, 0.2)
            dt = 1 / fs
            scales_high = pywt.central_frequency('cmor4.0-1.0') / (freqs_high * dt)

            # Baseline mean from first 100 min of NREM
            baseline_vals = []
            baseline_remaining = 100 * 60
            for b_start, b_end in bouts:
                if baseline_remaining <= 0:
                    break
                data_seg, fs_seg = extract_bout_data(b_start, b_end)
                if data_seg is None:
                    continue
                if target_fs is not None and fs_seg > target_fs:
                    n_samp = int(len(data_seg) * (target_fs / fs_seg))
                    data_seg = resample(data_seg, n_samp)
                    fs_seg = target_fs
                dt_seg = 1 / fs_seg
                scales_seg = pywt.central_frequency('cmor4.0-1.0') / (freqs_high * dt_seg)
                coeffs_seg, freqs_used = pywt.cwt(data_seg, scales_seg, 'cmor4.0-1.0', sampling_period=dt_seg)
                power_seg = abs(coeffs_seg) ** 2
                band_idx = (freqs_used >= 10) & (freqs_used <= 15)
                band_power = power_seg[band_idx].mean(axis=0)
                win_samples = int(4 * fs_seg)
                band_power = Series(band_power).rolling(win_samples, center=True, min_periods=1).mean().to_numpy()
                step = int(0.5 * fs_seg)
                band_power = band_power[::step]
                n_take = int(builtins.min(len(band_power), builtins.max(0, baseline_remaining / 0.5)))
                if n_take > 0:
                    baseline_vals.extend(band_power[:n_take].tolist())
                    baseline_remaining -= n_take * 0.5

            baseline_mean = nanmean(array(baseline_vals)) if len(baseline_vals) > 0 else nan
            if isnan(baseline_mean):
                logger.warning('Baseline mean is NaN; normalization may be invalid.')

            coeffs, freqs_used = pywt.cwt(eeg_bout, scales_high, 'cmor4.0-1.0', sampling_period=dt)
            power = abs(coeffs) ** 2
            band_idx = (freqs_used >= 10) & (freqs_used <= 15)
            band_power_raw = power[band_idx].mean(axis=0)
            norm_raw = band_power_raw / baseline_mean

            win_samples = int(4 * fs)
            norm_smooth = Series(norm_raw).rolling(win_samples, center=True, min_periods=1).mean().to_numpy()

            fig_c, (axc1, axc2) = plt.subplots(2, 1, figsize=(12, 5), sharex=True)
            axc1.plot(t_bout, norm_raw, color='k', linewidth=0.8)
            axc1.set_ylabel('Norm sigma')
            axc1.set_title('Normalized sigma power (raw)')
            axc2.plot(t_bout, norm_smooth, color='k', linewidth=0.8)
            axc2.set_ylabel('Norm sigma')
            axc2.set_xlabel('Time (s)')
            axc2.set_title('Normalized sigma power (4 s smoothing)')

            # ---- Save panels B/C ----
            fig_b.tight_layout()
            fig_c.tight_layout()

            fig_b.savefig(path.join(out_path, f'{sub}_{ses}_{fnamechan}_panel_B.png'), dpi=150)
            fig_c.savefig(path.join(out_path, f'{sub}_{ses}_{fnamechan}_panel_C.png'), dpi=150)

            plt.close(fig_b)
            plt.close(fig_c)


def _plot_infraslow_panels_all(clam_obj, subs = None, sessions = None,
                               stage = None, xml_dir = None, out_dir = None,
                               seed = 1, target_fs = 100,
                               logger = create_logger('Infraslow panels')):
    """Helper to iterate through tracking sheet and plot per sub/ses/channel."""
    subs = clam_obj.subs if subs is None else subs
    sessions = clam_obj.sessions if sessions is None else sessions

    if not isinstance(clam_obj.chan, DataFrame):
        logger.warning('Channel tracking sheet not available; specify sub/ses and channel.')
        return

    chan_df = clam_obj.chan
    if subs != 'all':
        if isinstance(subs, str):
            subs = [subs]
        chan_df = chan_df[chan_df['sub'].isin(subs)]
    if sessions != 'all':
        if isinstance(sessions, str):
            sessions = [sessions]
        chan_df = chan_df[chan_df['ses'].isin(sessions)]

    subs_list = sorted(chan_df['sub'].unique())
    for sub in subs_list:
        ses_list = sorted(chan_df[chan_df['sub'] == sub]['ses'].unique())
        for ses in ses_list:
            flag, chanset = load_channels(sub, ses, chan_df, chan_df, flag = 0, logger = logger)
            if not chanset:
                continue
            for ch in list(chanset.keys()):
                clam_obj.plot_infraslow_panels(sub=sub, ses=ses,
                                               stage=stage, channel=ch,
                                               xml_dir=xml_dir, out_dir=out_dir,
                                               seed=seed, target_fs=target_fs,
                                               logger=logger)


            
            
            
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
    event_indices = sort(random.choice(arange(100), 
                                             size=10, 
                                             replace=False))  # Random non-adjacent events
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
    
    # Set the sample period and keep original length for reporting.
    dt = 1 / fs
    n = len(eeg_signal)

    # 1) CWT of the EEG signal to get time-resolved power (0.5-24 Hz, 0.2 Hz steps).
    freqs_high = arange(0.5, 24.1, 0.2)
    scales_high = pywt.central_frequency(wavelet) / (freqs_high * dt)
    coeffs, freqs_used = pywt.cwt(eeg_signal, scales_high, wavelet, sampling_period=dt)
    power = abs(coeffs) ** 2

    results = {}

    # 2) Baseline PSD used to normalize each band (Welch on baseline_data).
    f_fft, pxx = welch(baseline_data, fs=fs, nperseg=fs*4)

    for band_name, (fmin, fmax) in band_defs.items():

        # 3) Band-power time series: average CWT power across frequencies in the band.
        band_idx = (freqs_used >= fmin) & (freqs_used <= fmax)
        band_power = power[band_idx].mean(axis=0)

        # 4) Smooth the band-power time series (moving average, default 4 s).
        win_samples = int(smoothing_window_sec * fs)
        smooth = Series(band_power).rolling(win_samples, center=True, min_periods=1).mean().to_numpy()

        # 5) Normalize by mean band power from baseline PSD.
        idx_band = (f_fft >= fmin) & (f_fft <= fmax)
        band_baseline = nanmean(pxx[idx_band])
        normalized = smooth / band_baseline

        # 6) Downsample to the infraslow temporal resolution (default 0.5 s).
        step = int(downsample_step_sec * fs)
        smooth_ds = normalized[::step]
        dt_ds = downsample_step_sec

        # 7) CWT of the band-power trace in the infraslow range (e.g., 0.001-0.12 Hz).
        freqs_low = arange(infraslow_freq_range[0], infraslow_freq_range[1], infraslow_step)
        scales_low = pywt.central_frequency(wavelet) / (freqs_low * dt_ds)
        coeffs_low, _ = pywt.cwt(smooth_ds, scales_low, wavelet, sampling_period=dt_ds)
        power_low = abs(coeffs_low) ** 2

        # 8) Average power across time to yield the infraslow PSD estimate.
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

def gaussian_offset(x, offset, a, mu, sigma):
    return offset + gaussian(x, a, mu, sigma)

def multi_gaussian(x, *params):
    """Sum of multiple Gaussian components."""
    y = zeros_like(x, dtype=float)
    for i in range(0, len(params), 3):
        a, mu, sigma = params[i:i+3]
        y += gaussian(x, a, mu, sigma)
    return y

def multi_gaussian_offset(x, offset, *params):
    """Sum of multiple Gaussian components with a constant offset."""
    return offset + multi_gaussian(x, *params)

def fit_gaussian(freqs_low, psd_vals, num_components=1, full_freqs=None,
                 init_peak_range=(0.0075, 0.04), init_high_range=(0.04, 0.1),
                 return_params=False, r2_threshold=0.2):
    """
    Fit a Gaussian (or sum of Gaussians) to a 1D power spectrum and compute peak characteristics.

    Parameters
    ----------
    freqs_low : array_like
        1D array of frequency values (e.g., in the infraslow range).
    psd_vals : array_like
        1D array of power spectral density (PSD) values corresponding to `freqs_low`.

    num_components : int
        Number of Gaussian components to fit. Defaults to 1.

    Returns
    -------
    peak_freq : float
        Estimated peak frequency (global maximum) of the fitted Gaussian mixture [Hz].
    sigma : float
        Standard deviation (spread) of the component nearest the peak [Hz].
    avg_peak_power : float
        Fitted power at the peak frequency (MATLAB-style).
        Returns NaN for all outputs if the fit fails.
    fit_curve : ndarray
        Fitted curve evaluated at `full_freqs` if provided, otherwise at `freqs_low`.
    fit_params : dict, optional
        Returned only if return_params=True. Contains 'offset' and 'components'.
    """
    
    try:
        freqs_low = array(freqs_low, dtype=float)
        psd_vals = array(psd_vals, dtype=float)
        x_eval = full_freqs if full_freqs is not None else freqs_low
        components = []
        mu_bounds = None

        # MATLAB-style preprocessing
        y_data = psd_vals.copy()
        y_data[freqs_low < 0.0075] = nanmedian(y_data)
        offset = nanmean(y_data[freqs_low >= 0.059999])
        if isnan(offset):
            offset = 0.0
        y_fit = y_data - offset

        # MATLAB-style initialization: seed mu/a from peak in 0.0075-0.04 Hz
        init_mask = (freqs_low > init_peak_range[0]) & (freqs_low <= init_peak_range[1])
        if init_mask.any():
            idx_init = nanargmax(y_fit[init_mask])
            mu_init = freqs_low[init_mask][idx_init]
            a_init = y_fit[init_mask][idx_init]
        else:
            idx_init = nanargmax(y_fit)
            mu_init = freqs_low[idx_init]
            a_init = y_fit[idx_init]

        # Robustify initial amplitude / resolution
        if not isfinite(a_init) or a_init <= 0:
            a_init = nanmax(y_fit)
        if not isfinite(a_init) or a_init <= 0:
            a_init = 1.0

        df = nanmean(diff(freqs_low))
        if not isfinite(df) or df <= 0:
            df = 0.001
        minsd = 2 * df

        def _adj_r2(y_true, y_pred, p):
            y_true = array(y_true, dtype=float)
            y_pred = array(y_pred, dtype=float)
            mask = isfinite(y_true) & isfinite(y_pred)
            if mask.sum() < (p + 2):
                return nan
            y_true = y_true[mask]
            y_pred = y_pred[mask]
            ss_res = ((y_true - y_pred) ** 2).sum()
            ss_tot = ((y_true - nanmean(y_true)) ** 2).sum()
            if ss_tot == 0:
                return nan
            r2 = 1 - ss_res / ss_tot
            n = len(y_true)
            return 1 - (1 - r2) * (n - 1) / (n - p - 1)

        def _fit_gauss1():
            mu_lo = 0.0075 + minsd
            mu_hi = 0.1 - minsd
            mu0 = builtins.max(mu_lo, builtins.min(mu_init, mu_hi))
            bounds = ([0.01 * a_init, mu_lo, minsd], [2 * a_init, mu_hi, 0.2])
            popt, _ = curve_fit(gaussian, freqs_low, y_fit,
                                p0=[a_init, mu0, 0.02],
                                bounds=bounds,
                                maxfev=8000)
            comp = [(popt[0], popt[1], popt[2])]
            curve = gaussian(x_eval, *popt) + offset
            return comp, curve, (mu_lo, mu_hi), _adj_r2(y_fit, gaussian(freqs_low, *popt), 3)

        def _fit_gauss2(mu1_bounds, mu2_bounds):
            mu1_lo, mu1_hi = mu1_bounds
            mu2_lo, mu2_hi = mu2_bounds
            mu0 = builtins.max(mu1_lo, builtins.min(mu_init, mu1_hi))
            mu2_init = 0.06
            mu2_init = builtins.max(mu2_lo, builtins.min(mu2_init, mu2_hi))
            p0 = [a_init, mu0, 0.02,
                  0.5 * a_init, mu2_init, 0.01]
            bounds_low = [0.01 * a_init, mu1_lo, minsd,
                          0.01 * a_init, mu2_lo, minsd]
            bounds_high = [2 * a_init, mu1_hi, 0.2,
                           2 * a_init, mu2_hi, 0.2]
            popt, _ = curve_fit(multi_gaussian, freqs_low, y_fit,
                                p0=p0,
                                bounds=(bounds_low, bounds_high),
                                maxfev=12000)
            comp = []
            for k in range(2):
                a_k, mu_k, sig_k = popt[3*k:3*k+3]
                comp.append((a_k, mu_k, sig_k))
            curve = multi_gaussian(x_eval, *popt) + offset
            return comp, curve, (mu1_bounds, mu2_bounds), _adj_r2(y_fit, multi_gaussian(freqs_low, *popt), 6)

        # Option 1: fit gauss1 and gauss2, then select by adj R2 (MATLAB-style)
        comp1, curve1, bounds1, adjr1 = _fit_gauss1()
        mu1_bounds = (0.0075 + minsd, 0.0333 - minsd)
        mu2_bounds = (0.0267 + minsd, 0.1)
        comp2, curve2, bounds2, adjr2 = _fit_gauss2(mu1_bounds, mu2_bounds)

        if isnan(adjr2) or (not isnan(adjr1) and adjr1 >= adjr2):
            components = comp1
            fit_curve = curve1
            mu_bounds = [bounds1]
            adjr = adjr1
            model_used = 1
        else:
            components = comp2
            fit_curve = curve2
            mu_bounds = [bounds2[0], bounds2[1]]
            adjr = adjr2
            model_used = 2

        # Option 2: flag poor fits (adjR2 or boundary)
        boundary_tol = minsd
        mu_lo, mu_hi = mu_bounds[0]
        at_bounds = (components[0][1] <= mu_lo + boundary_tol) or (components[0][1] >= mu_hi - boundary_tol)
        poor_fit = (isnan(adjr) or adjr < r2_threshold or at_bounds)

        # Option 3/4: only apply when fit is poor
        peak_source = 'fit'
        if poor_fit:
            # PSD peak in 0.0075-0.04 Hz
            psd_mask = (freqs_low > init_peak_range[0]) & (freqs_low <= init_peak_range[1])
            if psd_mask.any():
                idx_psd = nanargmax(psd_vals[psd_mask])
                psd_peak_freq = freqs_low[psd_mask][idx_psd]
                psd_peak_power = psd_vals[psd_mask][idx_psd]
            else:
                idx_psd = nanargmax(psd_vals)
                psd_peak_freq = freqs_low[idx_psd]
                psd_peak_power = psd_vals[idx_psd]

            # Option 4: relaxed bounds re-fit (both components allowed in 0.0075-0.1)
            mu_rel_bounds = (0.0075 + minsd, 0.1 - minsd)
            try:
                comp_rel, curve_rel, bounds_rel, adjr_rel = _fit_gauss2(mu_rel_bounds, mu_rel_bounds)
            except Exception:
                comp_rel, curve_rel, bounds_rel, adjr_rel = None, None, None, nan

            if isfinite(adjr_rel) and (isnan(adjr) or adjr_rel > adjr):
                components = comp_rel
                fit_curve = curve_rel
                mu_bounds = [bounds_rel[0], bounds_rel[1]]
                adjr = adjr_rel
                model_used = 2
                # pick component closest to PSD peak
                peak_comp = builtins.min(components, key=lambda c: abs(c[1] - psd_peak_freq))
                peak_freq = peak_comp[1]
                sigma_sel = peak_comp[2]
                peak_source = 'relaxed_fit'
            else:
                # fallback to PSD peak
                peak_freq = psd_peak_freq
                sigma_sel = nan
                avg_peak_power = psd_peak_power
                peak_source = 'psd'

        if not poor_fit:
            # MATLAB-style: use the first Gaussian component as the peak
            peak_freq = components[0][1]
            sigma_sel = components[0][2]

        interp_fit = interp1d(x_eval, fit_curve, kind='linear',
                              bounds_error=False, fill_value='extrapolate')
        if 'avg_peak_power' not in locals() or isnan(avg_peak_power):
            avg_peak_power = float(interp_fit(peak_freq))

    except RuntimeError:
        peak_freq, sigma_sel, avg_peak_power = nan, nan, nan
        fit_curve = full_freqs * nan if full_freqs is not None else freqs_low * nan

    fit_params = {'offset': offset if 'offset' in locals() else nan,
                  'components': components if 'components' in locals() else [],
                  'mu_bounds': mu_bounds if 'mu_bounds' in locals() else None,
                  'adjr2': adjr if 'adjr' in locals() else nan,
                  'model_used': model_used if 'model_used' in locals() else nan,
                  'poor_fit': poor_fit if 'poor_fit' in locals() else False,
                  'peak_source': peak_source if 'peak_source' in locals() else 'fit'}
    if return_params:
        return peak_freq, sigma_sel if 'sigma_sel' in locals() else nan, avg_peak_power, fit_curve, fit_params
    return peak_freq, sigma_sel if 'sigma_sel' in locals() else nan, avg_peak_power, fit_curve


def plot_psd_fit(freqs_low, psd_vals, fit_curve, fit_params, mu_bounds, outpath, title=None):
    """Plot PSD with fitted curves and mu bounds."""
    if freqs_low is None or psd_vals is None or fit_curve is None:
        return

    freqs_low = array(freqs_low, dtype=float)
    psd_vals = array(psd_vals, dtype=float)
    fit_curve = array(fit_curve, dtype=float)

    mask = (freqs_low >= 0) & (freqs_low <= 0.1)
    x = freqs_low[mask]
    y = psd_vals[mask]
    fit_y = fit_curve[mask]

    fig, ax = plt.subplots(1, 1, figsize=(6, 4))
    ax.plot(x, y, color='0.5', linewidth=1.0, label='PSD')
    ax.plot(x, fit_y, color='r', linewidth=1.2, label='Fit')

    offset = fit_params.get('offset', 0.0) if isinstance(fit_params, dict) else 0.0
    comps = fit_params.get('components', []) if isinstance(fit_params, dict) else []
    for a, mu, sigma in comps:
        ax.plot(x, offset + gaussian(x, a, mu, sigma), linestyle=':', color='r', linewidth=0.8)

    if mu_bounds is not None:
        if isinstance(mu_bounds[0], (list, tuple)):
            for lo, hi in mu_bounds:
                ax.axvline(lo, color='k', linestyle=':', linewidth=0.8)
                ax.axvline(hi, color='k', linestyle=':', linewidth=0.8)
        else:
            ax.axvline(mu_bounds[0], color='k', linestyle=':', linewidth=0.8)
            ax.axvline(mu_bounds[1], color='k', linestyle=':', linewidth=0.8)

    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('Power (a.u.)')
    ax.set_xlim(0, 0.1)
    if title:
        ax.set_title(title)
    ax.legend(loc='best', fontsize=8)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    
            
