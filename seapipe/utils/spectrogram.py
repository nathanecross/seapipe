#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
S.O.N.A.R. (Spectral Oscillation Normalisation And Representation)

Spectrogram utilities for seapipe.
"""

from __future__ import annotations

from os import listdir, mkdir, path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

from wonambi import Dataset
from wonambi.attr import Annotations
from wonambi.trans import fetch

from .logs import create_logger
from .load import load_channels, load_sessions, rename_channels
from ..spectrum.psa import default_event_opts


class SONAR:
    """S.O.N.A.R. spectrogram plotting."""

    def __init__(self, rootpath, rec_dir, xml_dir, out_dir, chan, ref_chan,
                 grp_name, stage, rater=None, subs='all', sessions='all',
                 event_opts=None, method='yasa', file_extension='.svg',
                 win_sec=30, fmin=0.5, fmax=25, trimperc=2.5,
                 cmap='RdBu_r', vmin=None, vmax=None, tracking=None):

        self.rootpath = rootpath
        self.rec_dir = rec_dir
        self.xml_dir = xml_dir
        self.out_dir = out_dir
        self.chan = chan
        self.ref_chan = ref_chan
        self.grp_name = grp_name
        self.stage = stage
        self.rater = rater
        self.subs = subs
        self.sessions = sessions
        self.method = method

        self.event_opts = event_opts if event_opts is not None else default_event_opts()
        self.file_extension = self._ensure_extension(file_extension)

        self.win_sec = win_sec
        self.fmin = fmin
        self.fmax = fmax
        self.trimperc = trimperc
        self.cmap = cmap
        self.vmin = vmin
        self.vmax = vmax

        if tracking is None:
            tracking = {'spectrogram': {}}
        self.tracking = tracking

    @staticmethod
    def _ensure_extension(ext):
        if not ext:
            return '.svg'
        return ext if ext.startswith('.') else f'.{ext}'

    @staticmethod
    def _segment_to_1d(seg):
        data = seg['data'].data[0][0]
        return np.asarray(data)

    def _compute_spectrogram(self, data, sf, method, logger):
        if method == 'yasa':
            try:
                from lspopt import spectrogram_lspopt
            except Exception as exc:
                logger.error(f"Failed to import lspopt for YASA spectrograms: {exc}")
                return None
            nperseg = int(self.win_sec * sf)
            if data.size <= 2 * nperseg:
                logger.warning('Data length too short for spectrogram window; skipping.')
                return None
            f, t, sxx = spectrogram_lspopt(data, sf, nperseg=nperseg, noverlap=0)
        elif method == 'scipy':
            try:
                from scipy.signal import spectrogram as scipy_spectrogram
            except Exception as exc:
                logger.error(f"Failed to import scipy.signal.spectrogram: {exc}")
                return None
            nperseg = int(self.win_sec * sf)
            if data.size <= 2 * nperseg:
                logger.warning('Data length too short for spectrogram window; skipping.')
                return None
            f, t, sxx = scipy_spectrogram(data, sf, nperseg=nperseg, noverlap=0)
        else:
            logger.error(f"Unknown spectrogram method: {method}")
            return None

        with np.errstate(divide='ignore'):
            sxx = 10 * np.log10(sxx)

        good = np.logical_and(f >= self.fmin, f <= self.fmax)
        if not np.any(good):
            logger.warning('No frequencies within the requested range; skipping.')
            return None
        f = f[good]
        sxx = sxx[good, :]
        return f, t, sxx

    def _get_norm_limits(self, sxx):
        if self.vmin is not None or self.vmax is not None:
            return self.vmin, self.vmax
        vmin, vmax = np.percentile(sxx, [0 + self.trimperc, 100 - self.trimperc])
        return vmin, vmax

    def _plot_matrix(self, f, t, sxx, logger):
        vmin, vmax = self._get_norm_limits(sxx)
        fig, ax = plt.subplots(nrows=1, figsize=(12, 4))
        im = ax.pcolormesh(t, f, sxx, norm=Normalize(vmin=vmin, vmax=vmax),
                           cmap=self.cmap, antialiased=True, shading='auto')
        if len(t) > 0:
            ax.set_xlim(t.min(), t.max())
        ax.set_ylabel('Frequency [Hz]')
        ax.set_xlabel('Time [sec]')
        cbar = fig.colorbar(im, ax=ax, shrink=0.95, fraction=0.1, aspect=25)
        cbar.ax.set_ylabel('Log Power (dB / Hz)', rotation=270, labelpad=20)
        return fig

    def _plot_spectrogram(self, data, sf, method, logger):
        if method == 'yasa':
            try:
                from yasa import plot_spectrogram as yasa_plot_spectrogram
            except Exception as exc:
                logger.error(f"Failed to import YASA plot_spectrogram: {exc}")
                return None
            try:
                fig = yasa_plot_spectrogram(
                    data,
                    sf,
                    win_sec=self.win_sec,
                    fmin=self.fmin,
                    fmax=self.fmax,
                    trimperc=self.trimperc,
                    cmap=self.cmap,
                    vmin=self.vmin,
                    vmax=self.vmax,
                )
            except Exception as exc:
                logger.error(f"Failed to plot YASA spectrogram: {exc}")
                return None
            return fig

        result = self._compute_spectrogram(data, sf, method, logger)
        if result is None:
            return None
        f, t, sxx = result
        return self._plot_matrix(f, t, sxx, logger)

    @staticmethod
    def _average_spectrograms(specs):
        min_f = min(spec[2].shape[0] for spec in specs)
        min_t = min(spec[2].shape[1] for spec in specs)
        f = specs[0][0][:min_f]
        t = specs[0][1][:min_t]
        stack = np.stack([spec[2][:min_f, :min_t] for spec in specs], axis=0)
        avg = np.nanmean(stack, axis=0)
        return f, t, avg

    def _group_segments(self, segments, cycle_idx, cat, stage_list):
        nsegs = []
        labels = []

        if cat[0] + cat[1] == 2:
            nsegs = [segments]
            labels = ["-".join(stage_list)]
        elif cat[0] + cat[1] == 0:
            for st in stage_list:
                for cy in cycle_idx:
                    segs = [s for s in segments if st in s['stage'] if cy in s['cycle']]
                    if segs:
                        nsegs.append(segs)
                        labels.append(f'{st}_cycle{cy}')
        elif cat[0] == 0:
            for cy in cycle_idx:
                segs = [s for s in segments if cy in s['cycle']]
                if segs:
                    nsegs.append(segs)
                    labels.append(f'cycle{cy}')
        elif cat[1] == 0:
            for st in stage_list:
                segs = [s for s in segments if st in s['stage']]
                if segs:
                    nsegs.append(segs)
                    labels.append(f'{st}')

        return nsegs, labels


    def spectrogram(self, cycle_idx=None, concat_stage=False, concat_cycle=True,
                    event_opts=None, method=None,
                    file_extension=None, filetype='.edf', progress=True,
                    logger=None):

        logger = logger or create_logger('SONAR')
        event_opts = event_opts if event_opts is not None else self.event_opts
        method = method or self.method
        file_extension = self._ensure_extension(file_extension or self.file_extension)

        buffer = event_opts.get('buffer', 0)
        if not isinstance(buffer, (int, float)):
            logger.critical("event_opts['buffer'] must be an int or float.")
            return
        
        if event_opts['evt_type']:
            cat = (int(concat_cycle), int(concat_stage), 0, 0)
        else:
            cat = (int(concat_cycle), int(concat_stage), 1, 0)
        if cat[0] == 0 and cycle_idx is None:
            logger.critical("To split cycles (concat_cycle=False), cycle_idx cannot be None.")
            return

        subs = self.subs
        if isinstance(subs, list):
            pass
        elif subs == 'all':
            subs = listdir(self.rec_dir)
            subs = [p for p in subs if not '.' in p]
        else:
            logger.error("'subs' must either be an array of subject ids or = 'all'")
            return

        subs.sort()
        for sub in subs:
            flag = 0
            if 'spectrogram' not in self.tracking:
                self.tracking['spectrogram'] = {}
            if sub not in self.tracking['spectrogram']:
                self.tracking['spectrogram'][sub] = {}

            flag, sessions = load_sessions(sub, self.sessions, self.rec_dir, flag,
                                           logger, verbose=2)
            for ses in sessions:
                logger.info('')
                logger.debug(f'Commencing {sub}, {ses}')
                if ses not in self.tracking['spectrogram'][sub]:
                    self.tracking['spectrogram'][sub][ses] = {}

                rdir = f'{self.rec_dir}/{sub}/{ses}/eeg/'
                try:
                    edf_file = [x for x in listdir(rdir) if x.endswith(filetype)]
                    dset = Dataset(rdir + edf_file[0])
                except Exception:
                    logger.warning(f'No input {filetype} file in {rdir}')
                    continue

                xdir = f'{self.xml_dir}/{sub}/{ses}/'
                try:
                    xml_file = [x for x in listdir(xdir) if x.endswith('.xml')]
                    annot = Annotations(xdir + xml_file[0], rater_name=self.rater)
                except Exception:
                    logger.warning(f'No input annotations file in {xdir}')
                    continue

                if cycle_idx is not None:
                    all_cycles = annot.get_cycles()
                    cycle = [all_cycles[y - 1] for y in cycle_idx if y <= len(all_cycles)]
                else:
                    cycle = None

                pflag = flag
                flag, chanset = load_channels(sub, ses, self.chan, self.ref_chan,
                                              flag, logger)
                if flag - pflag > 0:
                    logger.warning(f'Skipping {sub}, {ses}...')
                    continue

                newchans = rename_channels(sub, ses, self.chan, logger)

                for ch in chanset:
                    if newchans:
                        fnamechan = newchans[ch]
                    else:
                        fnamechan = ch

                    if not chanset[ch]:
                        logchan = ['(no re-referencing)']
                    else:
                        logchan = chanset[ch]

                    logger.debug(f"Reading EEG data for {sub}, {ses}, {str(ch)}:{'-'.join(logchan)}")
                    try:
                        if isinstance(event_opts['evt_type'], str):
                            event_opts['evt_type'] = [event_opts['evt_type']]
                        segments = fetch(dset, annot, cat=cat,
                                         evt_type=event_opts['evt_type'], stage=self.stage,
                                         cycle=cycle, buffer=buffer)
                        segments.read_data([ch], ref_chan=chanset[ch], grp_name=self.grp_name)
                    except Exception as exc:
                        logger.error(f'{exc}')
                        logger.warning(f'Skipping {sub}, {ses}, channel {str(ch)} ...')
                        continue

                    if len(segments) < 1:
                        logger.warning(f'No valid data found for {sub}, {ses}, {str(ch)}')
                        continue

                    stage_list = self.stage
                    if stage_list is None:
                        stages = []
                        for s in segments:
                            st = s['stage']
                            if isinstance(st, (list, tuple, set)):
                                stages.extend([x for x in st if x])
                            elif st:
                                stages.append(st)
                        stage_list = sorted(set(stages))

                    nsegs, labels = self._group_segments(segments, cycle_idx, cat, stage_list)
                    if not nsegs:
                        logger.warning(f'No segments to plot for {sub}, {ses}, {str(ch)}')
                        continue

                    outdir = self.out_dir
                    if not path.exists(outdir):
                        mkdir(outdir)
                    if not path.exists(f'{outdir}/{sub}'):
                        mkdir(f'{outdir}/{sub}')
                    if not path.exists(f'{outdir}/{sub}/{ses}'):
                        mkdir(f'{outdir}/{sub}/{ses}')
                    outpath = f'{outdir}/{sub}/{ses}'
                    
                    
                    for segs, label in zip(nsegs, labels):
                        if not segs:
                            continue

                        if event_opts['evt_type']:
                            specs = []
                            for seg in segs:
                                data = self._segment_to_1d(seg)
                                result = self._compute_spectrogram(data, dset.header['s_freq'],
                                                                   method, logger)
                                if result is not None:
                                    specs.append(result)
                            if not specs:
                                logger.warning(f'No spectrograms computed for {sub}, {ses}, {str(ch)}')
                                continue
                            f, t, avg = self._average_spectrograms(specs)
                            t = t - buffer
                            fig = self._plot_matrix(f, t, avg, logger)
                        else:
                            data = np.concatenate([self._segment_to_1d(s) for s in segs])
                            fig = self._plot_spectrogram(data, dset.header['s_freq'],
                                                         method, logger)

                        if fig is None:
                            continue

                        parts = [sub, ses, fnamechan]
                        if label:
                            parts.append(label)
                        if event_opts['evt_type']:
                            parts.append(str('_'.join(event_opts['evt_type'])))
                        fname = '_'.join(parts) + file_extension
                        logger.debug("Saving spectrogram file.")
                        fig.savefig(f'{outpath}/{fname}', dpi=300, bbox_inches='tight')
                        plt.close(fig)
                        
        logger.debug("Plot spectrogram finished without error.")

        return
