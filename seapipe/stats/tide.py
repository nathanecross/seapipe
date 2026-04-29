#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Transitions, Intervals, and Dynamics of Epochs (T.I.D.E.).
"""

from collections import Counter
from itertools import combinations
from os import listdir, makedirs, path

import numpy as np
from defusedxml import ElementTree as ET
from pandas import DataFrame

from ..utils.logs import create_logger


class tide:

    """Transitions, Intervals, and Dynamics of Epochs (T.I.D.E.)."""

    def __init__(self, xml_dir, out_dir=None, stage=None, rater=None,
                 subs='all', sessions='all', keyword=None,
                 subject_out_dir=None):

        self.xml_dir = xml_dir
        if out_dir is None:
            out_dir = path.expanduser('~/derivatives/datasets/hypnogram')
        self.out_dir = out_dir
        if subject_out_dir is None:
            subject_out_dir = self._default_subject_out_dir(out_dir)
        self.subject_out_dir = subject_out_dir
        self.stage = self._default_stage(stage)
        self.rater = rater
        self.subs = subs
        self.sessions = sessions
        self.keyword = keyword

    @staticmethod
    def _default_subject_out_dir(out_dir):
        out_dir = path.normpath(path.expanduser(out_dir))
        if (path.basename(out_dir) == 'hypnogram'
                and path.basename(path.dirname(out_dir)) == 'datasets'):
            return path.join(path.dirname(path.dirname(out_dir)), 'hypnogram')
        return out_dir

    @staticmethod
    def _default_stage(stage):
        if stage is None:
            return ['Wake', 'NREM1', 'NREM2', 'NREM3', 'REM']
        if isinstance(stage, str):
            return [stage]
        return stage

    @staticmethod
    def _stage_key(stage):
        key = {'W': 'Wake', 'N1': 'NREM1', 'N2': 'NREM2',
               'N3': 'NREM3', 'R': 'REM'}
        return key.get(stage, stage)

    @staticmethod
    def _sleep_stages():
        return ['NREM1', 'NREM2', 'NREM3', 'REM']

    @staticmethod
    def _nrem_stages():
        return ['NREM1', 'NREM2', 'NREM3']

    @staticmethod
    def _safe_divide(num, den):
        return np.nan if den == 0 else num / den

    @staticmethod
    def _skew(values):
        values = np.asarray(values, dtype=float)
        values = values[np.isfinite(values)]
        if values.size < 2:
            return np.nan
        mean = values.mean()
        std = values.std()
        if std == 0:
            return 0.0
        return np.mean(((values - mean) / std) ** 3)

    @staticmethod
    def _entropy(probabilities):
        probabilities = np.asarray(probabilities, dtype=float)
        probabilities = probabilities[np.isfinite(probabilities)]
        probabilities = probabilities[probabilities > 0]
        if probabilities.size == 0:
            return np.nan
        return float(-(probabilities * np.log2(probabilities)).sum())

    @staticmethod
    def _pearson(x, y):
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        mask = np.isfinite(x) & np.isfinite(y)
        x = x[mask]
        y = y[mask]
        if x.size < 2 or np.std(x) == 0 or np.std(y) == 0:
            return np.nan
        return float(np.corrcoef(x, y)[0, 1])

    def _subjects(self):
        if isinstance(self.subs, list):
            subs = self.subs
        elif self.subs == 'all':
            subs = [x for x in listdir(self.xml_dir) if not x.startswith('.')
                    and path.isdir(f'{self.xml_dir}/{x}')]
        else:
            raise TypeError("'subs' must be a list or 'all'.")
        subs.sort()
        return subs

    def _sessions(self, sub):
        if isinstance(self.sessions, list):
            sessions = self.sessions
        elif self.sessions == 'all':
            ses_dir = f'{self.xml_dir}/{sub}'
            sessions = [x for x in listdir(ses_dir) if not x.startswith('.')
                        and path.isdir(f'{ses_dir}/{x}')]
        else:
            raise TypeError("'sessions' must be a list or 'all'.")
        sessions.sort()
        return sessions

    def _xml_file(self, sub, ses, logger):
        xdir = f'{self.xml_dir}/{sub}/{ses}'
        if not path.exists(xdir):
            logger.warning(f'{xdir} does not exist. Skipping...')
            return None

        xml_files = [x for x in listdir(xdir) if x.endswith('.xml')
                     and not x.startswith('.')]
        if self.keyword:
            xml_files = [x for x in xml_files if self.keyword in x]

        if len(xml_files) == 0:
            logger.warning(f'No annotations file found for {sub}, {ses}. Skipping...')
            return None
        if len(xml_files) > 1:
            logger.warning(f'>1 annotations file found for {sub}, {ses}. '
                           "Use 'keyword' to select one. Skipping...")
            return None
        return f'{xdir}/{xml_files[0]}'

    def _out_ses_dir(self, sub, ses):
        out_ses = f'{self.subject_out_dir}/{sub}/{ses}'
        makedirs(out_ses, exist_ok=True)
        return out_ses

    def _load_hypnogram(self, xml_file):
        root = ET.parse(xml_file).getroot()
        rater = None

        if self.rater is None:
            rater = root.find('rater')
        else:
            for rat in root.iterfind('rater'):
                if rat.get('name') == self.rater:
                    rater = rat
                    break

        if rater is None:
            raise ValueError(f"Rater '{self.rater}' not found in {xml_file}")

        epochs = []
        for ep in rater.iterfind('stages/epoch'):
            start = int(float(ep.findtext('epoch_start')))
            end = int(float(ep.findtext('epoch_end')))
            stage = self._stage_key(ep.findtext('stage'))
            quality = ep.findtext('quality')
            epochs.append({'start': start, 'end': end,
                           'stage': stage, 'quality': quality})

        epochs.sort(key=lambda x: x['start'])
        return epochs

    @staticmethod
    def _epoch_stage_list(epochs):
        return [x['stage'] for x in epochs]

    @staticmethod
    def _epoch_durations(epochs):
        return [x['end'] - x['start'] for x in epochs]

    def _valid_transition_pairs(self, hypno, stage):
        stage = set(stage)
        return [(a, b) for a, b in zip(hypno[:-1], hypno[1:])
                if a in stage and b in stage]

    def _transition_counts(self, hypno, stage):
        counts = DataFrame(0, index=stage, columns=stage, dtype=float)
        for a, b in self._valid_transition_pairs(hypno, stage):
            counts.loc[a, b] += 1
        return counts

    @staticmethod
    def _transition_probabilities(counts):
        row_sums = counts.sum(axis=1)
        probs = counts.div(row_sums.replace(0, np.nan), axis=0)
        return probs.fillna(0)

    def _complete_transition_matrix(self, hypno, stage):
        counts = self._transition_counts(hypno, stage)
        return self._transition_probabilities(counts), counts

    def _reduced_transition_matrix(self, hypno):
        reduced = []
        for stg in hypno:
            if stg == 'Wake':
                reduced.append('Wake')
            elif stg in self._nrem_stages():
                reduced.append('NREM')
            elif stg == 'REM':
                reduced.append('REM')
            else:
                reduced.append(stg)

        stage = ['Wake', 'NREM', 'REM']
        counts = self._transition_counts(reduced, stage)
        return self._transition_probabilities(counts), counts

    def _transition_entropy(self, hypno, stage):
        counts = self._transition_counts(hypno, stage)
        total = counts.values.sum()
        if total == 0:
            return np.nan
        return self._entropy((counts / total).values.ravel())

    def _stage_shift_rate(self, hypno, stage):
        valid_epoch_count = sum(x in stage for x in hypno)
        if valid_epoch_count == 0:
            return np.nan
        transition_pairs = self._valid_transition_pairs(hypno, stage)
        stage_changes = sum(a != b for a, b in transition_pairs)
        return stage_changes / valid_epoch_count

    @staticmethod
    def _bouts(epochs):
        bouts = []
        if not epochs:
            return bouts

        cur_stage = epochs[0]['stage']
        cur_start = epochs[0]['start']
        cur_end = epochs[0]['end']
        for ep in epochs[1:]:
            if ep['stage'] == cur_stage and ep['start'] == cur_end:
                cur_end = ep['end']
            else:
                bouts.append({'stage': cur_stage, 'start': cur_start,
                              'end': cur_end, 'duration': cur_end - cur_start})
                cur_stage = ep['stage']
                cur_start = ep['start']
                cur_end = ep['end']
        bouts.append({'stage': cur_stage, 'start': cur_start,
                      'end': cur_end, 'duration': cur_end - cur_start})
        return bouts

    def _bout_duration_metrics(self, durations):
        durations = np.asarray(durations, dtype=float)
        return {
            'stage_mean_bout_dur_min': np.nan if durations.size == 0 else durations.mean(),
            'stage_median_bout_dur_min': np.nan if durations.size == 0 else np.median(durations),
            'stage_p75_bout_dur_min': np.nan if durations.size == 0 else np.percentile(durations, 75),
            'stage_skew_bout_dur': self._skew(durations),
            'stage_num_bouts': int(durations.size),
            'stage_prop_short_bouts': np.nan if durations.size == 0 else (durations < 2).sum() / durations.size,
        }

    @staticmethod
    def _stage_suffix(stage):
        return str(stage).replace(' ', '_')

    def _sleep_period_epochs(self, epochs):
        sleep = self._sleep_stages()
        sleep_idx = [i for i, ep in enumerate(epochs) if ep['stage'] in sleep]
        if not sleep_idx:
            return []
        return epochs[sleep_idx[0]:sleep_idx[-1] + 1]

    def _sleep_onset_epochs(self, epochs):
        sleep = self._sleep_stages()
        for i, ep in enumerate(epochs):
            if ep['stage'] in sleep:
                return epochs[i:]
        return []

    def _count_sleep_cycles(self, bouts, min_nrem_sec=15 * 60, min_rem_sec=2 * 60):
        bouts = [b for b in bouts if b['stage'] in self._sleep_stages()]
        cycles = 0
        nrem_since_rem = 0
        for bout in bouts:
            if bout['stage'] in self._nrem_stages():
                nrem_since_rem += bout['duration']
            elif (bout['stage'] == 'REM'
                  and bout['duration'] >= min_rem_sec
                  and nrem_since_rem >= min_nrem_sec):
                cycles += 1
                nrem_since_rem = 0
        return cycles

    def _proportion_in_half(self, epochs, target_stage, half):
        if not epochs:
            return np.nan
        midpoint = int(np.ceil(len(epochs) / 2))
        subset = epochs[:midpoint] if half == 'first' else epochs[midpoint:]
        if not subset:
            return np.nan
        return sum(ep['stage'] == target_stage for ep in subset) / len(subset)

    def transition_matrix(self, stage=None, resolution='complete',
                          logger=create_logger('TIDE transition matrix')):
        """Export per-subject/session transition matrices."""

        flag = 0
        stage = self.stage if stage is None else self._default_stage(stage)
        resolution = str(resolution).lower()
        if resolution not in ['complete', 'reduced']:
            logger.warning("resolution must be 'complete' or 'reduced'. "
                           "Defaulting to 'complete'.")
            resolution = 'complete'

        makedirs(self.out_dir, exist_ok=True)
        makedirs(self.subject_out_dir, exist_ok=True)
        summary = []

        for sub in self._subjects():
            for ses in self._sessions(sub):
                logger.debug(f'Calculating {resolution} transition matrix for {sub}, {ses}')
                xml_file = self._xml_file(sub, ses, logger)
                if not xml_file:
                    flag += 1
                    continue

                try:
                    epochs = self._load_hypnogram(xml_file)
                except Exception as err:
                    logger.warning(f'Could not read hypnogram for {sub}, {ses}: {err}')
                    flag += 1
                    continue

                hypno = self._epoch_stage_list(epochs)
                out_ses = self._out_ses_dir(sub, ses)

                if resolution == 'complete':
                    matrix, counts = self._complete_transition_matrix(hypno, stage)
                    matrix.to_csv(f'{out_ses}/{sub}_{ses}_tide_transition_matrix_complete.csv')
                    counts.to_csv(f'{out_ses}/{sub}_{ses}_tide_transition_counts_complete.csv')
                    row = {'sub': sub, 'ses': ses}
                    for from_stage in stage:
                        for to_stage in stage:
                            row[f'p_{from_stage}_to_{to_stage}'] = matrix.loc[from_stage, to_stage]
                    summary.append(row)
                else:
                    matrix, counts = self._reduced_transition_matrix(hypno)
                    matrix.to_csv(f'{out_ses}/{sub}_{ses}_tide_transition_matrix_reduced.csv')
                    counts.to_csv(f'{out_ses}/{sub}_{ses}_tide_transition_counts_reduced.csv')
                    row = {'sub': sub, 'ses': ses}
                    for from_stage in matrix.index:
                        for to_stage in matrix.columns:
                            row[f'p_{from_stage}_to_{to_stage}'] = matrix.loc[from_stage, to_stage]
                    summary.append(row)

        if summary:
            DataFrame(summary).to_csv(
                    f'{self.out_dir}/tide_transition_matrix_{resolution}_summary.csv',
                    index=False)

        if flag == 0:
            logger.debug('TIDE transition matrix export finished without ERROR.')
        else:
            logger.warning(f'TIDE transition matrix export finished with {flag}+ WARNINGS.')
        return

    def stage_duration_distributions(self, stage=None,
                                     logger=create_logger('TIDE stage durations')):
        """Export bout duration and global hypnogram dynamics metrics."""

        flag = 0
        stage = self.stage if stage is None else self._default_stage(stage)
        makedirs(self.out_dir, exist_ok=True)
        makedirs(self.subject_out_dir, exist_ok=True)
        summary = []

        for sub in self._subjects():
            for ses in self._sessions(sub):
                logger.debug(f'Calculating stage duration distributions for {sub}, {ses}')
                xml_file = self._xml_file(sub, ses, logger)
                if not xml_file:
                    flag += 1
                    continue

                try:
                    epochs = self._load_hypnogram(xml_file)
                except Exception as err:
                    logger.warning(f'Could not read hypnogram for {sub}, {ses}: {err}')
                    flag += 1
                    continue

                hypno = self._epoch_stage_list(epochs)
                bouts = self._bouts(epochs)
                sleep_period = self._sleep_period_epochs(epochs)
                transition_pairs = self._valid_transition_pairs(hypno, stage)
                p_stay = self._safe_divide(
                        sum(a == b for a, b in transition_pairs),
                        len(transition_pairs))
                stage_shift_rate = self._stage_shift_rate(hypno, stage)
                transition_entropy = self._transition_entropy(hypno, stage)
                num_sleep_cycles = self._count_sleep_cycles(bouts)
                rem_first = self._proportion_in_half(sleep_period, 'REM', 'first')
                rem_second = self._proportion_in_half(sleep_period, 'REM', 'second')
                n3_early = self._proportion_in_half(sleep_period, 'NREM3', 'first')
                n3_late = self._proportion_in_half(sleep_period, 'NREM3', 'second')
                n3_ratio = self._safe_divide(n3_early, n3_late)

                row = {'sub': sub, 'ses': ses}
                for stg in stage:
                    durations = [b['duration'] / 60 for b in bouts
                                 if b['stage'] == stg]
                    suffix = self._stage_suffix(stg)
                    for metric, value in self._bout_duration_metrics(durations).items():
                        row[f'{metric}_{suffix}'] = value

                all_stage_durations = [b['duration'] / 60 for b in bouts
                                       if b['stage'] in stage]
                for metric, value in self._bout_duration_metrics(all_stage_durations).items():
                    row[f'{metric}_all_stages'] = value

                row.update({
                    'p_stay_same_stage': p_stay,
                    'stage_shift_rate': stage_shift_rate,
                    'transition_entropy': transition_entropy,
                    'num_sleep_cycles': num_sleep_cycles,
                    'rem_first_half_prop': rem_first,
                    'rem_second_half_prop': rem_second,
                    'delta_n3_early_late_ratio': n3_ratio,
                })
                summary.append(row.copy())

                out_ses = self._out_ses_dir(sub, ses)
                DataFrame([row]).to_csv(
                        f'{out_ses}/{sub}_{ses}_tide_stage_duration_distributions.csv',
                        index=False)

        if summary:
            DataFrame(summary).to_csv(
                    f'{self.out_dir}/tide_stage_duration_distributions_summary.csv',
                    index=False)

        if flag == 0:
            logger.debug('TIDE stage duration export finished without ERROR.')
        else:
            logger.warning(f'TIDE stage duration export finished with {flag}+ WARNINGS.')
        return

    def _collect_hypnograms(self, stage, logger):
        hypnograms = []
        for sub in self._subjects():
            for ses in self._sessions(sub):
                xml_file = self._xml_file(sub, ses, logger)
                if not xml_file:
                    continue
                try:
                    epochs = self._load_hypnogram(xml_file)
                except Exception as err:
                    logger.warning(f'Could not read hypnogram for {sub}, {ses}: {err}')
                    continue

                sleep_onset_epochs = self._sleep_onset_epochs(epochs)
                if not sleep_onset_epochs:
                    logger.warning(f'No sleep onset found for {sub}, {ses}. Skipping...')
                    continue

                hypno = self._epoch_stage_list(sleep_onset_epochs)
                matrix, _ = self._complete_transition_matrix(hypno, stage)
                hypnograms.append({
                    'id': f'{sub}_{ses}',
                    'sub': sub,
                    'ses': ses,
                    'hypno': hypno,
                    'transition': matrix.values.ravel(),
                })
        return hypnograms

    def _align_pair(self, hyp_a, hyp_b, stage):
        n = min(len(hyp_a), len(hyp_b))
        if n == 0:
            return [], []
        stage = set(stage)
        a = hyp_a[:n]
        b = hyp_b[:n]
        valid = [(x, y) for x, y in zip(a, b) if x in stage and y in stage]
        if not valid:
            return [], []
        return [x for x, _ in valid], [y for _, y in valid]

    def _epoch_similarity(self, hyp_a, hyp_b, stage):
        a, b = self._align_pair(hyp_a, hyp_b, stage)
        if not a:
            return np.nan
        return sum(x == y for x, y in zip(a, b)) / len(a)

    def _kappa_similarity(self, hyp_a, hyp_b, stage):
        a, b = self._align_pair(hyp_a, hyp_b, stage)
        if not a:
            return np.nan
        observed = sum(x == y for x, y in zip(a, b)) / len(a)
        counts_a = Counter(a)
        counts_b = Counter(b)
        expected = sum((counts_a[stg] / len(a)) * (counts_b[stg] / len(b))
                       for stg in stage)
        if expected == 1:
            return 1.0 if observed == 1 else np.nan
        return (observed - expected) / (1 - expected)

    def hypnogram_similarity(self, stage=None,
                             logger=create_logger('TIDE hypnogram similarity')):
        """Export group-level pairwise hypnogram similarity matrices."""

        stage = self.stage if stage is None else self._default_stage(stage)
        makedirs(self.out_dir, exist_ok=True)
        hypnograms = self._collect_hypnograms(stage, logger)
        ids = [x['id'] for x in hypnograms]

        epoch = DataFrame(np.eye(len(ids)), index=ids, columns=ids)
        kappa = DataFrame(np.eye(len(ids)), index=ids, columns=ids)
        transition = DataFrame(np.eye(len(ids)), index=ids, columns=ids)

        for i, j in combinations(range(len(hypnograms)), 2):
            a = hypnograms[i]
            b = hypnograms[j]
            epoch_value = self._epoch_similarity(a['hypno'], b['hypno'], stage)
            kappa_value = self._kappa_similarity(a['hypno'], b['hypno'], stage)
            transition_value = self._pearson(a['transition'], b['transition'])

            epoch.iloc[i, j] = epoch.iloc[j, i] = epoch_value
            kappa.iloc[i, j] = kappa.iloc[j, i] = kappa_value
            transition.iloc[i, j] = transition.iloc[j, i] = transition_value

        epoch.to_csv(f'{self.out_dir}/hyp_sim_epoch.csv')
        kappa.to_csv(f'{self.out_dir}/hyp_sim_kappa.csv')
        transition.to_csv(f'{self.out_dir}/hyp_sim_transition_corr.csv')
        DataFrame([{'id': x['id'], 'sub': x['sub'], 'ses': x['ses'],
                    'n_epochs_from_sleep_onset': len(x['hypno'])}
                   for x in hypnograms]).to_csv(
                           f'{self.out_dir}/hypnogram_similarity_manifest.csv',
                           index=False)

        logger.debug('TIDE hypnogram similarity export finished.')
        return
