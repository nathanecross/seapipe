#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  7 15:47:20 2025

@author: ncro8394
"""

from numpy import (angle, arange, array, average, ceil, concatenate, degrees, diff, exp,
                   hamming, hanning, histogram,  inf, isnan, linspace, max, min, nan, nanargmax, nan_to_num,
                   nanargmin, nanmax, nanmean, nanmedian, nanstd, pi, random, sort, vstack, where, zeros, zeros_like)
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
                            min_psd_bout_sec = 300  # MATLAB script removes bouts shorter than 300 s
                            infraslow_step_hz = 0.001
                            infraslow_max_hz = 0.5  # compute up to 0.5, crop to 0.1 for fitting (Matlab style)

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

                        for seg_series in per_seg_series:
                            duration_sec = len(seg_series['times']) * step_sec
                            if duration_sec < min_psd_bout_sec:
                                continue
                            for band in freq_bands:
                                series = nan_to_num(array(seg_series[f'{band}_norm'], dtype=float),
                                                     nan=nanmean(seg_series[f'{band}_norm']))
                                
                                # High-pass to suppress slow drift (MATLAB-style FIR) or detrend if too short
                                hp_cut = 0.005
                                fs_band = 1 / step_sec
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
                                nperseg = int(300 / step_sec)
                                if len(centered) < nperseg:
                                    continue
                                noverlap = int(nperseg * 0.5)
                                freqs_pwel, power_pwel = welch(centered, fs=fs_band,
                                                               window='hann',
                                                               nperseg=nperseg,
                                                               noverlap=noverlap,
                                                               nfft=int(300 / step_sec),
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

                        stats = {}
                        for band in freq_bands:
                            arrays = array(psd_accum[band])
                            if len(arrays) == 0:
                                logger.warning(f'No valid bouts >= {min_psd_bout_sec}s for {band}; skipping.')
                                continue
                            weighted_avg = average(arrays, axis=0)
                            
                            # Normalize by mean power in 0.0075–0.1 Hz (Matlab style)
                            norm_mask = (freqs_low >= 0.0075) & (freqs_low <= 0.1)
                            mean_psd = weighted_avg / nanmean(weighted_avg[norm_mask])
                            
                            # No extra smoothing beyond PWELCH interpolation
                            fit_mask = (freqs_low >= 0.0075) & (freqs_low <= 0.1)
                            if not fit_mask.any():
                                logger.warning(f'No frequencies in fit window for {band}; skipping.')
                                continue
                            peak_freq, sigma_val, avg_peak_power, fit_curve = fit_gaussian(
                                freqs_low[fit_mask], mean_psd[fit_mask], num_components=2, full_freqs=freqs_low)

                            logger.debug(f"Chan: {fnamechan}: {chanset[ch]} - "
                                         f"{band} infraslow peak: {round(peak_freq,3)}; "
                                         f"Peak power: {round(avg_peak_power,3)}")

                            stats[f'{band}_peak_freq'] = round(peak_freq, 3)
                            stats[f'{band}_avg_peak_power'] = round(avg_peak_power, 3)
                            stats[f'{band}_sigma'] = round(sigma_val, 3)

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

def multi_gaussian(x, *params):
    """Sum of multiple Gaussian components."""
    y = zeros_like(x, dtype=float)
    for i in range(0, len(params), 3):
        a, mu, sigma = params[i:i+3]
        y += gaussian(x, a, mu, sigma)
    return y

def fit_gaussian(freqs_low, psd_vals, num_components=1, full_freqs=None):
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
        Average PSD value in the frequency window centered at the peak frequency
        and spanning ±0.5 * sigma, computed via linear interpolation.
        Returns NaN for all outputs if the fit fails.
    fit_curve : ndarray
        Fitted curve evaluated at `full_freqs` if provided, otherwise at `freqs_low`.
    """
    
    try:
        x_eval = full_freqs if full_freqs is not None else freqs_low
        components = []
        if num_components == 1:
            bounds = ([0, 0.001, 1e-4], [inf, 0.12, 0.1])  # [a, mu, sigma]
            popt, _ = curve_fit(gaussian, freqs_low, psd_vals, 
                                p0 = [nanmax(psd_vals), 0.02, 0.01],
                                bounds = bounds,
                                maxfev = 8000)
            a, mu, sigma = popt
            components.append((a, mu, sigma))
            fit_curve = gaussian(x_eval, *popt)
        else:
            centers = linspace(0.01, 0.09, num_components)
            p0 = []
            bounds_low = []
            bounds_high = []
            for c in centers:
                p0.extend([nanmax(psd_vals), c, 0.01])
                bounds_low.extend([0, 0.001, 1e-4])
                bounds_high.extend([inf, 0.12, 0.1])
            popt, _ = curve_fit(multi_gaussian, freqs_low, psd_vals,
                                p0 = p0,
                                bounds = (bounds_low, bounds_high),
                                maxfev = 12000)
            for k in range(num_components):
                a_k, mu_k, sig_k = popt[3*k:3*k+3]
                components.append((a_k, mu_k, sig_k))
            fit_curve = multi_gaussian(x_eval, *popt)

        # Peak from fitted curve (global max), then choose nearest component for sigma
        peak_idx = nanargmax(fit_curve)
        peak_freq = x_eval[peak_idx]
        if len(components) == 1:
            sigma_sel = components[0][2]
        else:
            mu_candidates = [c[1] for c in components]
            idx_near = nanargmin([abs(m - peak_freq) for m in mu_candidates])
            sigma_sel = components[idx_near][2]

        peak_band = (peak_freq - 0.5 * sigma_sel, peak_freq + 0.5 * sigma_sel)
        
        # Clip window to valid frequency range
        peak_band = (
                    builtins.max((peak_band[0], freqs_low.min())),
                    builtins.min((peak_band[1], freqs_low.max()))
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
        peak_freq, sigma_sel, avg_peak_power = nan, nan, nan
        fit_curve = full_freqs * nan if full_freqs is not None else freqs_low * nan

    return peak_freq, sigma_sel if 'sigma_sel' in locals() else nan, avg_peak_power, fit_curve
    
            
