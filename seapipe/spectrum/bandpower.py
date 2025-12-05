#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compute band-limited power timecourses across subjects, sessions, and channels.

This module mirrors the iteration / filtering logic from psa.py but outputs a
time-resolved power estimate for a user-defined frequency band either via
sliding windows or per discrete sleep epochs.
"""
from copy import deepcopy
from dataclasses import dataclass
from math import isnan
from os import listdir, mkdir, path, walk
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from numpy import asarray, nan, trapz
from pandas import DataFrame
from scipy.signal import spectrogram, welch
from wonambi import Dataset
from wonambi.attr import Annotations
from wonambi.trans import fetch

from ..utils.load import infer_ref, load_channels, rename_channels
from ..utils.logs import create_logger
from ..utils.misc import bandpass_mne, laplacian_mne, notch_mne, notch_mne2
from .psa import (Spectrum, default_epoch_opts, default_event_opts,
                  default_filter_opts, default_general_opts)

Band = Tuple[float, float]


def default_bandpower_opts() -> Dict[str, float]:
    """Global defaults for the bandpower_timecourse method."""
    return {
        'mode': 'sliding',        # 'sliding' or 'epoch'
        'window_duration': 4.0,   # seconds, only used for sliding windows
        'window_step': 1.0,       # seconds between windows
        'epoch_duration': 30.0,   # seconds, convenience override for epoch mode
        'epoch_overlap': 0.0,     # seconds, when epoch mode piggybacks on fetch
        'min_segment': 2.0,       # minimum valid data duration in seconds
        'normalize': False,       # divide band power by baseline mean band power
        'baseline_cat': (0, 1, 1, 1),  # category tuple for baseline pull
        'baseline_duration': 100.0,    # seconds of baseline to use (from start)
    }


@dataclass
class BandWindow:
    """Container describing a single band-power estimate."""
    window_start: float
    window_end: float
    window_center: float
    band_power: float


class BandPowerTimecourse(Spectrum):

    """
    Iterate through the BIDS hierarchy and output band-limited power values
    across time either from sliding windows or discrete 30 s epochs.
    """

    def bandpower_timecourse(
        self,
        band: Band,
        general_opts: Optional[dict] = None,
        filter_opts: Optional[dict] = None,
        epoch_opts: Optional[dict] = None,
        bandpower_opts: Optional[dict] = None,
        event_opts: Optional[dict] = None,
        filetype: str = '.edf',
        split_by_stage: bool = False,
        logger=create_logger('Band power timecourse'),
    ):
        """
        Parameters
        ----------
        band : tuple(float, float)
            Lower and upper frequency limits whose broadband power will be
            quantified across time.
        general_opts / filter_opts / epoch_opts / event_opts :
            Same dictionaries used by Spectrum, see psa.py for full details.
            Only the subset of keys involved in IO + filtering is required.
        bandpower_opts :
            Dict returned by default_bandpower_opts controlling whether the
            power is computed on sliding windows or discrete epochs.
            When bandpower_opts['normalize'] is True, band power is divided by
            the mean band power computed from a baseline segment.
        filetype :
            EEG recording extension under <sub>/<ses>/eeg (defaults to '.edf').
        split_by_stage :
            If True, writes one CSV per sleep stage; otherwise one per channel/session.
        logger :
            Optional custom logger. When not provided a module-level logger is
            created for you.

        Outputs one CSV per channel containing the window start / end / center
        and band-limited power, allowing downstream plotting or statistics.
        """
        general_opts = deepcopy(general_opts or default_general_opts())
        filter_opts = deepcopy(filter_opts or default_filter_opts())
        epoch_opts = deepcopy(epoch_opts or default_epoch_opts())
        event_opts = deepcopy(event_opts or default_event_opts())
        bandpower_opts = deepcopy(bandpower_opts or default_bandpower_opts())

        _validate_band(band)
        mode = bandpower_opts['mode']
        if mode not in ('sliding', 'epoch'):
            raise ValueError("bandpower_opts['mode'] must be 'sliding' or 'epoch'")

        # Epoch mode is explicit about chunking the data, override defaults if asked.
        if mode == 'epoch':
            epoch_opts['epoch'] = True
            epoch_opts['epoch_dur'] = bandpower_opts['epoch_duration']
            epoch_opts['epoch_overlap'] = bandpower_opts['epoch_overlap']
            epoch_opts['epoch_step'] = None

        output_root = _prepare_output_root(self.out_dir, band, mode)
        subs = self.subs
        if subs == 'all':
            subs = next(walk(self.rec_dir))[1]
        elif not isinstance(subs, Iterable):
            raise ValueError("'subs' must either be a list of ids or 'all'")
        subs = sorted(subs)

        for sub in subs:
            sessions = self.sessions
            if sessions == 'all':
                sessions = [x for x in listdir(f'{self.rec_dir}/{sub}') if '.' not in x]
            sessions = sorted(sessions)

            for ses in sessions:
                logger.debug(f'Bandpower across time: {sub}, {ses}')
                dset, annot = _load_recording(self.rec_dir, self.xml_dir,
                                              sub, ses, filetype,
                                              self.rater, logger)
                if dset is None or annot is None:
                    continue

                cycle = _load_cycles(annot, self.cycle_idx)
                flag, chanset = load_channels(sub, ses, self.chan,
                                              self.ref_chan, 0, logger)
                if flag or not chanset:
                    logger.warning(f'Skipping {sub}, {ses} because channels could not be resolved')
                    continue

                rename_map = rename_channels(sub, ses, self.chan, logger) or {}
                channel_group = general_opts.get('chan_grp_name', 'eeg')

                for ch, ref in chanset.items():
                    channel_rows: List[dict] = []
                    fnamechan = rename_map.get(ch, ch)
                    channel_filter_opts = deepcopy(filter_opts)
                    if not channel_filter_opts.get('oREF'):
                        channel_filter_opts['oREF'] = infer_ref(sub, ses, self.chan, logger)
                    laplace_flag = bool(channel_filter_opts['laplacian'] and
                                        channel_filter_opts.get('lapchan'))
                    _baseline_band_mean = None

                    if isinstance(ref, list):
                        ref_label = '/'.join(ref)
                    elif ref:
                        ref_label = str(ref)
                    else:
                        ref_label = ''
                    logger.debug(f"Fetching data for {sub}, {ses}, {ch}:{ref_label or '[]'}")

                    if bandpower_opts.get('normalize'):
                        _baseline_band_mean = _compute_baseline_band_mean(
                            dset=dset,
                            annot=annot,
                            ch=ch,
                            ref=ref,
                            channel_filter_opts=channel_filter_opts,
                            laplace_flag=laplace_flag,
                            band=band,
                            band_opts=bandpower_opts,
                            stage=self.stage,
                            cycle=cycle,
                            logger=logger,
                        )
                    try:
                        segments = fetch(
                            dset,
                            annot,
                            cat=self.cat,
                            evt_type=event_opts['evt_type'],
                            stage=self.stage,
                            cycle=cycle,
                            epoch=epoch_opts['epoch'],
                            epoch_dur=epoch_opts['epoch_dur'],
                            epoch_overlap=epoch_opts['epoch_overlap'],
                            epoch_step=epoch_opts['epoch_step'],
                            reject_epoch=epoch_opts['reject_epoch'],
                            reject_artf=epoch_opts['reject_artf'],
                            min_dur=epoch_opts['min_dur'],
                            buffer=event_opts['buffer'],
                        )
                    except Exception as error:
                        logger.error(f'Failed to fetch data for {sub}, {ses}, {ch}: {error}')
                        continue

                    if len(segments) == 0:
                        logger.warning(f'No segments for {sub}, {ses}, {ch}')
                        continue

                    selectchans = channel_filter_opts['lapchan'] if laplace_flag else [ch]
                    if isinstance(selectchans, str):
                        selectchans = [selectchans]
                    try:
                        segments.read_data(selectchans, ref)
                    except Exception as error:
                        logger.error(f'Unable to read data for {ch}:{error}')
                        continue

                    for seg in segments:
                        data = seg['data']
                        segment_rows = _process_segment_bandpower(
                            seg=seg,
                            band=band,
                            data=data,
                            ch=ch,
                    ref_chan=ref,
                            channel_group=channel_group,
                            filter_opts=channel_filter_opts,
                            laplace_flag=laplace_flag,
                            mode=mode,
                            band_opts=bandpower_opts,
                            baseline_mean=_baseline_band_mean if _baseline_band_mean else None,
                            logger=logger,
                        )
                        if segment_rows:
                            for row in segment_rows:
                                row.update({
                                    'sub': sub,
                                    'ses': ses,
                                    'channel': fnamechan,
                                    'stage': seg.get('stage'),
                                    'cycle': seg.get('cycle'),
                                    'evt_type': seg.get('evt_type'),
                                })
                            channel_rows.extend(segment_rows)

                    if channel_rows:
                        outdir = _ensure_session_dir(output_root, sub, ses)
                        band_label = f"{band[0]:g}-{band[1]:g}Hz"

                        if split_by_stage:
                            rows_by_stage = {}
                            for row in channel_rows:
                                stage_label = str(row.get('stage') or 'NA')
                                rows_by_stage.setdefault(stage_label, []).append(row)
                            for stage_label, rows in rows_by_stage.items():
                                rows.sort(
                                    key=lambda row: (
                                        row.get('window_start', row.get('segment_start', 0.0)),
                                        row.get('window_end', row.get('segment_end', 0.0)),
                                    )
                                )
                                safe_stage = stage_label.replace(' ', '').replace('/', '-')
                                outfile = (f"{outdir}/{sub}_{ses}_{fnamechan}_"
                                           f"{safe_stage}_{mode}_{band_label}.csv")
                                DataFrame(rows).to_csv(outfile, index=False)
                                logger.debug(f'Wrote {len(rows)} rows to {outfile}')
                        else:
                            channel_rows.sort(
                                key=lambda row: (
                                    row.get('window_start', row.get('segment_start', 0.0)),
                                    row.get('window_end', row.get('segment_end', 0.0)),
                                )
                            )
                            if isinstance(self.stage, list):
                                stagename = '-'.join(self.stage)
                            else:
                                stagename = str(self.stage) if self.stage else 'all'
                            stagename = stagename.replace(' ', '')
                            outfile = f"{outdir}/{sub}_{ses}_{fnamechan}_{stagename}_{mode}_{band_label}.csv"
                            DataFrame(channel_rows).to_csv(outfile, index=False)
                            logger.debug(f'Wrote {len(channel_rows)} rows to {outfile}')


def _validate_band(band: Sequence[float]):
    if len(band) != 2 or band[0] >= band[1]:
        raise ValueError('band must be a tuple (low, high) with low < high')


def _prepare_output_root(base_out: str, band: Band, mode: str) -> str:
    label = f"{band[0]:g}-{band[1]:g}Hz_{mode}"
    out_root = f'{base_out}/bandpower_timeseries_{label}'
    if not path.exists(base_out):
        mkdir(base_out)
    if not path.exists(out_root):
        mkdir(out_root)
    return out_root


def _ensure_session_dir(out_root: str, sub: str, ses: str) -> str:
    sub_dir = f'{out_root}/{sub}'
    if not path.exists(sub_dir):
        mkdir(sub_dir)
    ses_dir = f'{sub_dir}/{ses}'
    if not path.exists(ses_dir):
        mkdir(ses_dir)
    return ses_dir


def _load_recording(rec_dir: str, xml_dir: str, sub: str, ses: str,
                    filetype: str, rater: Optional[str], logger):
    rdir = f'{rec_dir}/{sub}/{ses}/eeg/'
    xdir = f'{xml_dir}/{sub}/{ses}/'
    try:
        edf_file = [x for x in listdir(rdir) if x.endswith(filetype)][0]
        dset = Dataset(f'{rdir}{edf_file}')
    except Exception:
        logger.warning(f'No {filetype} file found in {rdir}')
        return None, None
    try:
        xml_file = [x for x in listdir(xdir) if x.endswith('.xml')][0]
        annot = Annotations(f'{xdir}{xml_file}', rater_name=rater)
    except Exception:
        logger.warning(f'No annotation file found in {xdir}')
        return None, None
    return dset, annot


def _load_cycles(annot: Optional[Annotations], cycle_idx: Optional[Sequence[int]]):
    if annot is None or cycle_idx is None:
        return None
    all_cycles = annot.get_cycles()
    return [all_cycles[i - 1] for i in cycle_idx if i <= len(all_cycles)]


def _process_segment_bandpower(
    seg: dict,
    band: Band,
    data,
    ch: str,
    ref_chan: Optional[Sequence[str]],
    channel_group: str,
    filter_opts: dict,
    laplace_flag: bool,
    mode: str,
    band_opts: dict,
    baseline_mean: Optional[float],
    logger,
) -> List[dict]:
    """Apply filters, compute PSDs, and summarise band power for one segment."""
    timeline = data.axis['time'][0]
    segment_start = float(timeline[0])
    segment_end = float(timeline[-1])
    duration = segment_end - segment_start

    if duration < band_opts['min_segment']:
        logger.debug('Skipping segment shorter than min_segment')
        return []

    selectchans = filter_opts['lapchan'] if laplace_flag else [ch]
    if isinstance(selectchans, str):
        selectchans = [selectchans]
    filtflag = _apply_filters(data, selectchans, ref_chan, filter_opts,
                              laplace_flag, ch, logger)
    if filtflag:
        return []

    signal = _extract_signal(data, ch, channel_group)
    if signal is None:
        logger.warning('Unable to locate requested channel in data block')
        return []

    if mode == 'sliding':
        windows = _sliding_band_power(
            signal,
            data.s_freq,
            segment_start,
            band_opts['window_duration'],
            band_opts['window_step'],
            band,
            baseline_mean=baseline_mean,
        )
    else:
        windows = [_epoch_band_power(signal, data.s_freq, segment_start,
                                     segment_end, band, baseline_mean)]
    rows: List[dict] = []
    for window in windows:
        if isnan(window.band_power):
            continue
        rows.append({
            'segment_start': segment_start,
            'segment_end': segment_end,
            'window_start': window.window_start,
            'window_end': window.window_end,
            'window_center': window.window_center,
            'band_power': window.band_power,
        })
    return rows


def _extract_signal(data, ch: str, grp: Optional[str]):
    """Return a 1-D numpy array containing the channel of interest."""
    chan_axis = list(data.axis['chan'][0])
    candidates = [ch]
    if grp:
        candidates.append(f'{ch} ({grp})')

    for label in candidates:
        if label in chan_axis:
            idx = chan_axis.index(label)
            return data.data[0][idx]
    if chan_axis:
        return data.data[0][0]
    return None


def _sliding_band_power(
    signal,
    s_freq: float,
    segment_start: float,
    window_duration: float,
    window_step: float,
    band: Band,
    baseline_mean: Optional[float] = None,
) -> List[BandWindow]:
    """Use scipy.signal.spectrogram to return band power windows."""
    nperseg = max(int(round(window_duration * s_freq)), 1)
    step_samples = max(int(round(window_step * s_freq)), 1)
    if len(signal) < nperseg:
        return []
    noverlap = max(nperseg - step_samples, 0)
    freqs, times, Sxx = spectrogram(
        signal,
        fs=s_freq,
        nperseg=nperseg,
        noverlap=noverlap,
        scaling='density',
        detrend='constant',
        mode='psd',
    )
    mask = (freqs >= band[0]) & (freqs <= band[1])
    if not mask.any():
        return []
    band_power = trapz(Sxx[mask], freqs[mask], axis=0)
    if baseline_mean and baseline_mean > 0:
        band_power = band_power / baseline_mean
    windows: List[BandWindow] = []
    for center, value in zip(times, band_power):
        abs_center = segment_start + float(center)
        window_start = abs_center - window_duration / 2
        window_end = window_start + window_duration
        windows.append(BandWindow(
            window_start=window_start,
            window_end=window_end,
            window_center=abs_center,
            band_power=float(value),
        ))
    return windows


def _epoch_band_power(signal, s_freq: float,
                      start: float, end: float, band: Band,
                      baseline_mean: Optional[float]) -> BandWindow:
    """One value per epoch using Welch PSD over the full segment."""
    nperseg = min(len(signal), max(int(s_freq * 2), 256))
    freqs, psd = welch(signal, fs=s_freq, nperseg=nperseg, scaling='density')
    mask = (freqs >= band[0]) & (freqs <= band[1])
    if not mask.any():
        value = nan
    else:
        value = float(trapz(psd[mask], freqs[mask]))
        if baseline_mean and baseline_mean > 0:
            value = value / baseline_mean
    center = start + (end - start) / 2
    return BandWindow(window_start=start, window_end=end,
                      window_center=center, band_power=value)


def _apply_filters(data, selectchans, ref_chan, filter_opts, laplace_flag, ch, logger):
    """Run notch/bandpass/laplacian stack; return filtflag (0 ok, >0 error)."""
    filtflag = 0
    if filter_opts['notch']:
        data.data[0], filtflag = notch_mne(
            data,
            channel=selectchans,
            freq=filter_opts['notch_freq'],
            oREF=filter_opts['oREF'],
            rename=filter_opts['laplacian_rename'],
            renames=filter_opts['renames'],
            montage=filter_opts['montage'],
        )
    if filtflag == 0 and filter_opts['notch_harmonics']:
        data.data[0], filtflag = notch_mne2(
            data,
            channel=selectchans,
            oREF=filter_opts['oREF'],
            rename=filter_opts['laplacian_rename'],
            renames=filter_opts['renames'],
            montage=filter_opts['montage'],
        )
    if filtflag == 0 and filter_opts['bandpass']:
        data.data[0], filtflag = bandpass_mne(
            data,
            channel=selectchans,
            highpass=filter_opts['highpass'],
            lowpass=filter_opts['lowpass'],
            oREF=filter_opts['oREF'],
            rename=filter_opts['laplacian_rename'],
            renames=filter_opts['renames'],
            montage=filter_opts['montage'],
        )
    if filtflag == 0 and laplace_flag:
        data.data[0], filtflag = laplacian_mne(
            data,
            channel=selectchans,
            ref_chan=ref_chan,
            oREF=filter_opts['oREF'],
            laplacian_rename=filter_opts['laplacian_rename'],
            renames=filter_opts['renames'],
            montage=filter_opts['montage'],
        )
        data.axis['chan'][0] = asarray([ch])
    if filtflag:
        logger.warning('Filtering failed, skipping segment')
    return filtflag


def _compute_baseline_band_mean(dset, annot, ch, ref, channel_filter_opts,
                                laplace_flag, band, band_opts, stage, cycle,
                                logger):
    """Fetch baseline data and return mean band power for normalization."""
    try:
        baseline_segments = fetch(
            dset,
            annot,
            cat=band_opts.get('baseline_cat', (0, 1, 1, 1)),
            stage=stage,
            cycle=cycle,
            reject_artf=True,
        )
    except Exception as error:
        logger.warning(f'Baseline fetch failed: {error}')
        return None
    if not baseline_segments:
        logger.debug('No baseline segments found for normalization')
        return None
    selectchans = channel_filter_opts['lapchan'] if laplace_flag else [ch]
    if isinstance(selectchans, str):
        selectchans = [selectchans]
    try:
        baseline_segments.read_data(selectchans, ref)
    except Exception as error:
        logger.warning(f'Baseline read failed: {error}')
        return None
    base_seg = baseline_segments[0]['data']
    base_copy = deepcopy(base_seg)
    filtflag = _apply_filters(base_copy, selectchans, ref, channel_filter_opts,
                              laplace_flag, ch, logger)
    if filtflag:
        return None
    sig = _extract_signal(base_copy, ch, None)
    if sig is None:
        return None
    dur_sec = band_opts.get('baseline_duration', 0) or 0
    if dur_sec > 0:
        max_samples = int(dur_sec * base_copy.s_freq)
        sig = sig[:max_samples]
    nperseg = min(len(sig), max(int(base_copy.s_freq * 2), 256))
    freqs, psd = welch(sig, fs=base_copy.s_freq, nperseg=nperseg,
                       scaling='density')
    mask = (freqs >= band[0]) & (freqs <= band[1])
    if not mask.any():
        return None
    mean_val = float(trapz(psd[mask], freqs[mask]))
    return mean_val if mean_val > 0 else None
