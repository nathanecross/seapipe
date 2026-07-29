"""Event synchrony utilities.

Author: Jordan O'Byrne & formatted by Nathan Cross
Refactored into CORAL by Codex (2026)
"""

from copy import deepcopy
from numpy import asarray, nan, ones
from os import listdir, mkdir, path
from pandas import DataFrame
from shutil import copy
from wonambi import Dataset
from wonambi.attr import Annotations
from wonambi.detect import match_events
from wonambi.trans import fetch
from seapipe.utils.logs import create_logger
from ..utils.load import load_channels, load_sessions, rename_channels


class CORAL:

    """ Co-Occurrence and Recall Analysis of Labels (CORAL).

    This module divides a target event type from XML into events with and
    without a probe (co-occurrence based).
    """

    def __init__(self, rootpath=None, rec_dir=None, xml_dir=None, out_dir=None,
                 chan=None, ref_chan=None, stage=None, grp_name='eeg',
                 rater=None, subs='all', sessions='all', reject_artf=None,
                 filetype=None, tracking=None):

        self.rootpath = rootpath
        self.rec_dir = rec_dir
        self.xml_dir = xml_dir
        self.out_dir = out_dir
        self.chan = chan
        self.ref_chan = ref_chan
        self.stage = stage
        self.grp_name = grp_name
        self.rater = rater
        self.subs = subs
        self.sessions = sessions
        self.reject = reject_artf if reject_artf is not None else [
            'Artefact', 'Arou', 'Arousal'
        ]
        if filetype is None:
            filetype = ('.edf', '.rec', '.eeg')
        elif isinstance(filetype, list):
            filetype = tuple(filetype)
        self.filetype = filetype

        if tracking is None:
            tracking = {'sync': {}}
        self.tracking = tracking

    def _resolve_rec_dir(self, sub, ses):
        rdir = f'{self.rec_dir}/{sub}/{ses}/'
        if path.exists(f'{rdir}/eeg'):
            rdir = f'{rdir}/eeg/'
        return rdir

    def event_sync(self, evttype_target, evttype_probe, iu_thresh,
                   evttype_tp_target, evttype_fn=None, cat=(0, 0, 0, 0),
                   subs=None, sessions=None, chan=None, stage=None,
                   grp_name=None, rater=None, logger=None):

        '''
        This function divides a target event type from XML into events with and
        without a probe.
        '''
        flag = 0 
        if not logger:
            logger = create_logger('Event syncrony')

        if path.exists(self.out_dir):
            logger.debug("Output directory: " + self.out_dir + " exists")
        else:
            mkdir(self.out_dir)

        # Check input list
        if subs is None:
            subs = self.subs
        if isinstance(subs, list):
            None
        elif subs == 'all':
            subs = listdir(self.rec_dir)
            subs = [p for p in subs if '.' not in p]
        else:
            logger.error("ERROR: 'subs' must either be an array of subject ids or = 'all'")
            return

        if sessions is None:
            sessions = self.sessions

        logger.debug(
            rf""" Merging events...
                                 ⣀
                                ⠘⢷⣤⣿⡇
                        ⢰⡗⠀⠀⢠⡀⣠⡄⠀⠈⣿
                   ⠸⢶⣤⣄⢿⡇⠀⠀⠈⣿⠏⠀⠀⠀⣿⡀⠀⣴⠟⠁
                      ⠙⠻⣿⣦⡀⢸⡏⠀⠀⠀⠀⢹⣇⣼⣏⣀⣀⣠⣤⡦
                  ⢰⣶⡄⠀⠀⠘⢿⣿⣾⣧⠀⠀⠀⠀⣼⣿⠟⠉⠉⠉⢉⡀
                ⠘⠷⠶⢿⣿⡄⠀⠀⠀⠙⠿⣿⣦⣄⡀⣼⣿⠃⠀⠰⣦⣀⣸⡇
                   ⠈⣿⣷⠀⠀⠀⢀⠀⠈⠛⠿⣿⣿⡏⠀⠀⠀⠈⣹⡟⠀⠀⣀⣤⣄
                ⢠⡶⠟⠻⣿⣧⡀⣰⠏⠀⠀⠀⠀⢸⣿⡇⠀⣀⣠⣾⣯⣶⣶⣾⡏⠁
                     ⠈⠻⣿⣿⣀⠀⠀⠀⠀⠸⣿⣷⣿⣿⣿⣯⣉⡉⠉⠙⢷⡄
                       ⠈⠙⢿⣷⣦⣄⣠⣾⣿⠋⠀⠀⠀⠈⣩⡿⠷⣤
                    ⠶⠶⣶⡶⠟⠛⠛⢿⣿⡿⠁⠀⠀⠀⠀⣰⡟⠁
                      ⠹⠇⠀⠀⠀⣼⣿⠃⠀⠀⠀⠀⠀⠉⠁
                           ⣿⣿
                           ⠉⠉

⠀⠀⠀⠀⠀⠀⠀⠀⠀Co-Occurrence and Recall Analysis of Labels 
                       (C.O.R.A.L)
            
          Combining events "{evttype_target}" with "{evttype_probe}" into "{evttype_tp_target}"
          Consesus threshold = {iu_thresh}
            

                 """,
        )



        # Check for channel specification
        if chan is None:
            chan = self.chan
        if isinstance(chan, str):
            chan = [chan]

        if stage is None:
            stage = self.stage
        if isinstance(stage, str):
            stage = [stage]
        if not stage:
            stage = [None]

        if grp_name is None:
            grp_name = self.grp_name
        if rater is None:
            rater = self.rater

        # Concatenation
        cat = list(cat)
        cat[2] = 0
        cat[3] = 0  # force event type non-concatenation (required for analysis)
        cat = tuple(cat)

        multi_stage = len(stage) > 1
        if cat[1] == 1 and multi_stage:
            stage = [stage]

        # Begin loop through participants
        subs.sort()
        for sub in subs:
            if 'sync' not in self.tracking:
                self.tracking['sync'] = {}
            if sub not in self.tracking['sync'].keys():
                self.tracking['sync'][sub] = {}

            # Begin loop through sessions
            pflag = deepcopy(flag)
            flag, sub_sessions = load_sessions(sub, sessions, self.rec_dir,
                                               flag, logger, verbose=2)
            if flag - pflag > 0 or not sub_sessions:
                logger.warning(f'Skipping {sub}...')
                continue

            for ses in sub_sessions:
                logger.info('')
                logger.debug(f'Commencing {sub}, {ses}')
                if ses not in self.tracking['sync'][sub].keys():
                    self.tracking['sync'][sub][ses] = {}

                # Load recording
                rdir = self._resolve_rec_dir(sub, ses)
                try:
                    edf_file = [x for x in listdir(rdir)
                                if x.endswith(self.filetype)]
                    dset = Dataset(rdir + edf_file[0])
                except Exception:
                    logger.warning(f' No input {self.filetype} file in {rdir}')
                    flag += 1
                    continue

                # Load annotations
                xdir = f'{self.xml_dir}/{sub}/{ses}/'
                try:
                    xml_file = [x for x in listdir(xdir) if x.endswith('.xml')]
                    if not path.exists(f'{self.out_dir}/{sub}'):
                        mkdir(f'{self.out_dir}/{sub}')
                    if not path.exists(f'{self.out_dir}/{sub}/{ses}'):
                        mkdir(f'{self.out_dir}/{sub}/{ses}')

                    outpath = f'{self.out_dir}/{sub}/{ses}'
                    backup_file = f'{outpath}/{sub}_{ses}_sync.xml'
                    if not path.exists(backup_file):
                        copy(xdir + xml_file[0], backup_file)
                    else:
                        logger.debug(f'Using annotations file: {backup_file}')
                except Exception:
                    logger.warning(f' No input annotations file in {xdir}')
                    flag += 1
                    continue

                annot = Annotations(backup_file, rater_name=rater)

                # Channel setup
                pflag = deepcopy(flag)
                flag, chanset = load_channels(sub, ses, chan, self.ref_chan,
                                              flag, logger)
                if flag - pflag > 0:
                    logger.warning(f'Skipping {sub}, {ses}...')
                    flag += 1
                    continue

                newchans = rename_channels(sub, ses, chan, logger)

                # Loop through channels
                for ch in chanset:
                    chan_full = f'{ch} ({grp_name})'
                    fnamechan = newchans.get(ch, ch) if newchans else ch
                    logger.debug(f'Channel {fnamechan}')

                    for stg in stage:
                        logger.debug(f'Stage {stg}')
                        stage_arg = stg if isinstance(stg, list) else [stg] if stg is not None else None

                        # Get events - target
                        segments = fetch(
                            dset,
                            annot,
                            cat=cat,
                            evt_type=[evttype_target],
                            stage=stage_arg,
                            chan_full=[chan_full],
                            reject_epoch=True,
                            reject_artf=self.reject,
                        )

                        if len(segments) < 1:
                            logger.warning(
                                f"No valid targets found for {sub}, {ses}, {ch}. "
                                "Skipping ..."
                            )
                            flag += 1
                            continue

                        evt_target_seg = segments.segments
                        logger.debug(f"#targets '{evt_target_seg}' = {len(evt_target_seg)}")
                        for evt_idx, seg in enumerate(evt_target_seg):
                            evt_target_seg[evt_idx]['start'] = seg['times'][0][0]
                            evt_target_seg[evt_idx]['end'] = seg['times'][0][1]
                            evt_target_seg[evt_idx]['chan'] = [seg['chan']]
                            evt_target_seg[evt_idx]['quality'] = 'Good'

                        # Get events - probe
                        segments = fetch(
                            dset,
                            annot,
                            cat=cat,
                            evt_type=[evttype_probe],
                            stage=stage_arg,
                            chan_full=[chan_full],
                            reject_epoch=True,
                            reject_artf=self.reject,
                        )

                        if len(segments) < 1:
                            logger.warning(
                                f"No valid probes found for {sub}, {ses}, {ch}. "
                                "Skipping ..."
                            )
                            flag += 1
                            continue

                        evt_probe_seg = segments.segments
                        logger.debug(f"#probes '{evt_probe_seg}' = {len(evt_probe_seg)}")
                        for evt_idx, seg in enumerate(evt_probe_seg):
                            evt_probe_seg[evt_idx]['start'] = seg['times'][0][0]
                            evt_probe_seg[evt_idx]['end'] = seg['times'][0][1]
                            evt_probe_seg[evt_idx]['chan'] = [seg['chan']]
                            evt_probe_seg[evt_idx]['quality'] = 'Good'

                        # Assess co-occurrence and write output events to XML.
                        matched = match_events(
                            evt_probe_seg, evt_target_seg, iu_thresh
                        )
                        matched.to_annot(annot, 'tp_std', evttype_tp_target)

                        if evttype_fn is not None:
                            matched.to_annot(annot, 'fn', evttype_fn)

        ### 3. Check completion status and print
        if flag == 0:
            logger.debug('Event synchrony finished without error.')  
        else:
            logger.warning(f'Event synchrony finished with {flag} WARNINGS. See log for details.')
            
        return

    def event_sync_dataset(self, outfile_suffix, evttype_target, evttype_probe,
                           iu_thresh, evttype_tp_target, evttype_fn=None,
                           cat=(0, 0, 0, 0), subs=None, sessions=None,
                           chan=None, stage=None, grp_name=None, rater=None,
                           logger=None):

        '''
        This function accepts 2 event types: {target} and {probe}
        The {target} is split into 2 depending on the presence of the probe at
        the same time.
        '''
        if not logger:
            logger = create_logger('Event syncrony dataset')

        if path.exists(self.out_dir):
            logger.info(self.out_dir + " already exists")
        else:
            mkdir(self.out_dir)

        # Loop through records
        if subs is None:
            subs = self.subs
        if isinstance(subs, list):
            None
        elif subs == 'all':
            subs = listdir(self.xml_dir)
            subs = [p for p in subs if '.' not in p]
        else:
            logger.error("ERROR: 'subs' must either be an array of subject ids or = 'all'")
            return

        if sessions is None:
            sessions = self.sessions
        if sessions == 'all':
            sessions = [x for x in listdir(self.xml_dir + '/' + subs[0])
                        if '.' not in x]

        subs.sort()
        sessions.sort()

        logger.info("Extracting events and creating dataset...")

        # Create output
        out_base = self.out_dir if self.out_dir.endswith('/') else self.out_dir + '/'
        stats_header = ['Recall', 'Precision', 'F1 score', evttype_tp_target,
                        f'{evttype_probe}-', evttype_fn]
        stats = ones((len(subs) * len(sessions), len(stats_header))) * nan

        # Check for channel specification
        if chan is None:
            chan = self.chan
        if not chan:
            chan = [None]

        if stage is None:
            stage = self.stage
        if not stage:
            stage = [None]

        if grp_name is None:
            grp_name = self.grp_name
        if rater is None:
            rater = self.rater

        # Concatenation
        cat = list(cat)
        cat[2] = 0
        cat[3] = 0  # force event type non-concatenation (required for analysis)
        cat = tuple(cat)

        multi_stage = len(stage) > 1
        if cat[1] == 1 and multi_stage:
            stage = [stage]

        for stg in stage:
            for channel in chan:
                ids = []
                stage_arg = stg if isinstance(stg, list) else [stg] if stg is not None else None
                stage_label = None if stage_arg is None else '-'.join(stage_arg)
                logger.info(f'Channel {channel}; Stage {stage_arg}')

                if channel:
                    if stage_label:
                        stats_file = (out_base + channel + '_' + stage_label + '_'
                                      + outfile_suffix)
                    else:
                        stats_file = out_base + channel + '_' + outfile_suffix
                    chan_full = channel + ' (' + grp_name + ')'
                else:
                    if stage_label:
                        stats_file = out_base + stage_label + '_' + outfile_suffix
                    else:
                        stats_file = out_base + outfile_suffix
                    chan_full = channel

                for sub_idx, p in enumerate(subs):
                    for ses_idx, vis in enumerate(sessions):

                        logger.info(f'Subject: {p}, Visit: {vis}')
                        logger.info(f'{chan_full}')
                        ids.append(f'{p}_{vis}')

                        # Define files
                        rdir = self._resolve_rec_dir(p, vis)
                        xdir = f'{self.xml_dir}/{p}/{vis}/'
                        edf_file = [x for x in listdir(rdir)
                                    if x.endswith(self.filetype)]
                        xml_file = [x for x in listdir(xdir)
                                    if x.endswith('.xml')]

                        # Import data
                        dset = Dataset(rdir + edf_file[0])
                        annot = Annotations(xdir + xml_file[0], rater_name=rater)

                        # Get events - target
                        segments = fetch(
                            dset,
                            annot,
                            cat=cat,
                            evt_type=[evttype_target],
                            stage=stage_arg,
                            chan_full=[chan_full],
                            reject_epoch=True,
                            reject_artf=self.reject,
                        )
                        evt_target_seg = segments.segments
                        logger.info(f'#targets = {len(evt_target_seg)}')
                        for evt_idx, seg in enumerate(evt_target_seg):
                            evt_target_seg[evt_idx]['start'] = seg['times'][0][0]
                            evt_target_seg[evt_idx]['end'] = seg['times'][0][1]
                            evt_target_seg[evt_idx]['chan'] = [seg['chan']]
                            evt_target_seg[evt_idx]['quality'] = 'Good'

                        # Get events - probe
                        segments = fetch(
                            dset,
                            annot,
                            cat=cat,
                            evt_type=[evttype_probe],
                            stage=stage_arg,
                            chan_full=[chan_full],
                            reject_epoch=True,
                            reject_artf=self.reject,
                        )
                        evt_probe_seg = segments.segments
                        logger.info(f'#probes = {len(evt_probe_seg)}')
                        for evt_idx, seg in enumerate(evt_probe_seg):
                            evt_probe_seg[evt_idx]['start'] = seg['times'][0][0]
                            evt_probe_seg[evt_idx]['end'] = seg['times'][0][1]
                            evt_probe_seg[evt_idx]['chan'] = [seg['chan']]
                            evt_probe_seg[evt_idx]['quality'] = 'Good'

                        # Assess co-occurrence
                        matched = match_events(
                            evt_probe_seg, evt_target_seg, iu_thresh
                        )

                        # Wonambi's fp/fn counters can be inconsistent for some
                        # overlap patterns, so derive summary stats directly
                        # from the TP matrix and event-list lengths.
                        tp_matrix = asarray(
                            matched.tp, dtype=bool
                        ).reshape((len(evt_probe_seg), len(evt_target_seg)))
                        n_tp = int(tp_matrix.sum())
                        n_fp = int(len(evt_probe_seg) - tp_matrix.any(axis=1).sum())
                        n_fn = int(len(evt_target_seg) - tp_matrix.any(axis=0).sum())

                        if n_tp + n_fn == 0:
                            recall = 0.0
                        else:
                            recall = n_tp / (n_tp + n_fn)

                        if n_tp + n_fp == 0:
                            precision = 0.0
                        else:
                            precision = n_tp / (n_tp + n_fp)

                        if precision + recall == 0:
                            f1score = 0.0
                        else:
                            f1score = 2 * precision * recall / (precision + recall)

                        # Store stats in master table
                        row_idx = (sub_idx * len(sessions)) + ses_idx
                        stats[row_idx, 0] = recall
                        stats[row_idx, 1] = precision
                        stats[row_idx, 2] = f1score
                        stats[row_idx, 3] = n_tp
                        stats[row_idx, 4] = n_fp
                        stats[row_idx, 5] = n_fn

                # Export stats table
                logger.info(f'Saving {stats_file}')
                df = DataFrame(data=stats, index=ids, columns=stats_header)
                df.to_csv(stats_file)


def event_sync(rec_dir, xml_dir, out_dir, part, visit, cat,
               evttype_target, evttype_probe, iu_thresh, evttype_tp_target,
               evttype_fn, chan=None, stage=None, grp='eeg', rater=None,
               reject_artf=None, filetype=None, logger=None):

    coral = CORAL(
        rootpath=None,
        rec_dir=rec_dir,
        xml_dir=xml_dir,
        out_dir=out_dir,
        chan=chan,
        stage=stage,
        grp_name=grp,
        rater=rater,
        subs=part,
        sessions=visit,
        reject_artf=reject_artf,
        filetype=filetype,
    )
    coral.event_sync(evttype_target, evttype_probe, iu_thresh,
                     evttype_tp_target, evttype_fn, cat=cat, logger=logger)


def event_sync_dataset(rec_dir, xml_dir, out_dir, part, visit, cat, outfile_suffix,
                       evttype_target, evttype_probe, iu_thresh, evttype_tp_target,
                       evttype_fn, chan=None, stage=None, grp='eeg', rater=None,
                       reject_artf=None, filetype=None, logger=None):

    coral = CORAL(
        rootpath=None,
        rec_dir=rec_dir,
        xml_dir=xml_dir,
        out_dir=out_dir,
        chan=chan,
        stage=stage,
        grp_name=grp,
        rater=rater,
        subs=part,
        sessions=visit,
        reject_artf=reject_artf,
        filetype=filetype,
    )
    coral.event_sync_dataset(outfile_suffix, evttype_target, evttype_probe,
                             iu_thresh, evttype_tp_target, evttype_fn,
                             cat=cat, logger=logger)
