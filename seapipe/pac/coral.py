"""Event synchrony utilities.

Author: Jordan O'Byrne & formatted by Nathan Cross
Refactored into CORAL by Codex (2026)
"""

from numpy import nan, ones
from os import listdir, mkdir, path
from pandas import DataFrame
from shutil import copy
from wonambi import Dataset
from wonambi.attr import Annotations
from wonambi.detect import match_events
from wonambi.trans import fetch


class CORAL:

    """ Co-Occurrence and Recall Analysis of Labels (CORAL).

    This module divides a target event type from XML into events with and
    without a probe (co-occurrence based).
    """

    def __init__(self, rootpath=None, rec_dir=None, xml_dir=None, out_dir=None,
                 chan=None, stage=None, grp_name='eeg', rater=None, subs='all',
                 sessions='all', reject_artf=None, filetype=None, tracking=None):

        self.rootpath = rootpath
        self.rec_dir = rec_dir
        self.xml_dir = xml_dir
        self.out_dir = out_dir
        self.chan = chan
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

    def _log(self, msg, logger=None, level='info'):
        if logger is None:
            print(msg)
            return
        if level == 'error':
            logger.error(msg)
        elif level == 'warning':
            logger.warning(msg)
        else:
            logger.info(msg)

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
        
        if path.exists(self.out_dir):
            self._log(self.out_dir + " already exists", logger)
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
            self._log("ERROR: 'subs' must either be an array of subject ids or = 'all'",
                      logger, level='error')
            return

        if sessions is None:
            sessions = self.sessions
        if sessions == 'all':
            sessions = [x for x in listdir(self.xml_dir + '/' + subs[0])
                        if '.' not in x]

        subs.sort()
        sessions.sort()

        self._log(
            f""" Merging events...
            
                    ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀     ⠀⣀⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
                        ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠘⢷⣤⣿⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
                        ⠀⠀⠀⠀⠀⠀⠀⢰⡗⠀⠀⢠⡀⣠⡄⠀⠈⣿⠀⠀⠀⢀⠀⠀⠀⠀⠀⠀⠀⠀
                        ⠀⠀⠀⠸⢶⣤⣄⢿⡇⠀⠀⠈⣿⠏⠀⠀⠀⣿⡀⠀⣴⠟⠁⠀⠀⠀⠀⠀⠀⠀
                        ⠀⠀⠀⠀⠀⠀⠙⠻⣿⣦⡀⢸⡏⠀⠀⠀⠀⢹⣇⣼⣏⣀⣀⣠⣤⡦⠀⠀⠀⠀
                        ⠀⠀⠀⢰⣶⡄⠀⠀⠘⢿⣿⣾⣧⠀⠀⠀⠀⣼⣿⠟⠉⠉⠉⢉⡀⠀⠀⠀⠀⠀
                        ⠀⠘⠷⠶⢿⣿⡄⠀⠀⠀⠙⠿⣿⣦⣄⡀⣼⣿⠃⠀⠰⣦⣀⣸⡇⠀⠀⠀⠀⠀
                        ⠀⠀⠀⠀⠈⣿⣷⠀⠀⠀⢀⠀⠈⠛⠿⣿⣿⡏⠀⠀⠀⠈⣹⡟⠀⠀⣀⣤⣄⠀
                        ⠀⠀⢠⡶⠟⠻⣿⣧⡀⣰⠏⠀⠀⠀⠀⢸⣿⡇⠀⣀⣠⣾⣯⣶⣶⣾⡏⠁⠀⠀
                        ⠀⠀⠀⠀⠀⠀⠈⠻⣿⣿⣀⠀⠀⠀⠀⠸⣿⣷⣿⣿⣿⣯⣉⡉⠉⠙⢷⡄⠀⠀
                        ⠀⠀⠀⠀⠀⠀⠀⠀⠈⠙⢿⣷⣦⣄⣠⣾⣿⠋⠀⠀⠀⠈⣩⡿⠷⣤⠀⠀⠀⠀
                        ⠀⠀⠀⠀⠀⠀⠶⠶⣶⡶⠟⠛⠛⢿⣿⡿⠁⠀⠀⠀⠀⣰⡟⠁⠀⠀⠀⠀⠀⠀
                        ⠀⠀⠀⠀⠀⠀⠀⠀⠹⠇⠀⠀⠀⣼⣿⠃⠀⠀⠀⠀⠀⠉⠁⠀⠀⠀⠀⠀⠀⠀
                        ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣿⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
        ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀               ⠉⠉

⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀Co-Occurrence and Recall Analysis of Labels 
                            (C.O.R.A.L)
            
                Combining events "{evttype_target}" with "{evttype_probe}" into "{evttype_tp_target}"
                Consesus threshold = {iu_thresh}
            

                 """,
            logger,
        )

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

        if cat[1] == 1 and len(stage) > 1:
            stage = [stage]

        for stg in stage:
            for channel in chan:
                ids = []
                self._log(f'Channel {channel}; Stage {stg}', logger)

                chan_full = channel
                if channel:
                    chan_full = channel + ' (' + grp_name + ')'

                for i, sub in enumerate(subs):
                    for v, ses in enumerate(sessions):

                        self._log(f'Subject: {sub}, Visit: {ses}', logger)
                        self._log(f'{chan_full}', logger)
                        ids.append(f'{sub}_{ses}')

                        # Define files
                        rdir = self._resolve_rec_dir(sub, ses)
                        xdir = f'{self.xml_dir}/{sub}/{ses}/'
                        edf_file = [x for x in listdir(rdir)
                                    if x.endswith(self.filetype)]
                        xml_file = [x for x in listdir(xdir)
                                    if x.endswith('.xml')]

                        # Copy annotations file before beginning
                        if not path.exists(self.out_dir):
                            mkdir(self.out_dir)
                        if not path.exists(f'{self.out_dir}/{sub}'):
                            mkdir(f'{self.out_dir}/{sub}')
                        if not path.exists(f'{self.out_dir}/{sub}/{ses}'):
                            mkdir(f'{self.out_dir}/{sub}/{ses}')
                            backup = f'{self.out_dir}/{sub}/{ses}/'
                            backup_file = (f'{backup}{sub}_{ses}_spindles.xml')
                            copy(xdir + xml_file[0], backup_file)
                        else:
                            backup = f'{self.out_dir}/{sub}/{ses}/'
                            backup_file = (f'{backup}{sub}_{ses}_spindles.xml')

                        # Import data
                        try:
                            dset = Dataset(rdir + edf_file[0])
                            annot = Annotations(backup_file, rater_name=rater)
                        except Exception as exc:  
                            logger.warning(f"Failed to load data: {exc}")
                            flag +=1
                            continue
    
                        # Get events - target
                        segments = fetch(
                            dset,
                            annot,
                            cat=cat,
                            evt_type=[evttype_target],
                            stage=stg,
                            chan_full=[chan_full],
                            reject_epoch=True,
                            reject_artf=self.reject,
                        )

                        if len(segments) < 1:
                            logger.warning(f"No valid targets found for {sub}, {ses}, {channel}"
                                           f"Skipping ... ")
                            flag +=1
                            continue
                        
                        evt_target_seg = segments.segments
                        self._log(f'#targets = {len(evt_target_seg)}', logger)
                        for i, seg in enumerate(evt_target_seg):
                            evt_target_seg[i]['start'] = seg['times'][0][0]
                            evt_target_seg[i]['end'] = seg['times'][0][1]
                            evt_target_seg[i]['chan'] = [seg['chan']]
                            evt_target_seg[i]['quality'] = 'Good'

                        # Get events - probe
                        segments = fetch(
                            dset,
                            annot,
                            cat=cat,
                            evt_type=[evttype_probe],
                            stage=stg,
                            chan_full=[chan_full],
                            reject_epoch=True,
                            reject_artf=self.reject,
                        )
                        
                        if len(segments) < 1:
                            logger.warning(f"No valid probes found for {sub}, {ses}, {channel}"
                                           f"Skipping ... ")
                            flag +=1
                            continue
                        
                        evt_probe_seg = segments.segments
                        self._log(f'#probes = {len(evt_probe_seg)}', logger)
                        for i, seg in enumerate(evt_probe_seg):
                            evt_probe_seg[i]['start'] = seg['times'][0][0]
                            evt_probe_seg[i]['end'] = seg['times'][0][1]
                            evt_probe_seg[i]['chan'] = [seg['chan']]
                            evt_probe_seg[i]['quality'] = 'Good'

                        # Assess co-occurrence
                        matched = match_events(
                            evt_probe_seg, evt_target_seg, iu_thresh
                        )

                        # Write true positive target events to xml
                        matched.to_annot(annot, 'tp_std', evttype_tp_target)

                        if evttype_fn is not None:
                            matched.to_annot(annot, 'fn', evttype_fn)
                            
                            
        ### 3. Check completion status and print
        if flag == 0:
            logger.debug('Spindle detection finished without error.')  
        else:
            logger.warning(f'Spindle detection finished with {flag} WARNINGS. See log for details.')
            
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

        if path.exists(self.out_dir):
            self._log(self.out_dir + " already exists", logger)
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
            self._log("ERROR: 'subs' must either be an array of subject ids or = 'all'",
                      logger, level='error')
            return

        if sessions is None:
            sessions = self.sessions
        if sessions == 'all':
            sessions = [x for x in listdir(self.xml_dir + '/' + subs[0])
                        if '.' not in x]

        subs.sort()
        sessions.sort()

        self._log("Extracting events and creating dataset...", logger)

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

        if cat[1] == 1 and len(stage) > 1:
            stage = [stage]

        for stg in stage:
            for channel in chan:
                ids = []
                self._log(f'Channel {channel}; Stage {stg}', logger)

                if channel:
                    if not stg:
                        stats_file = out_base + channel + '_' + outfile_suffix
                    elif len(stg) < 2:
                        stats_file = (out_base + channel + '_' + stg[0]
                                      + '_' + outfile_suffix)
                    else:
                        if type(stg) is not list:
                            stg = [stg]
                        stages = '_'.join(stg)
                        stats_file = (out_base + channel + '_' + stages + '_'
                                      + outfile_suffix)
                    chan_full = channel + ' (' + grp_name + ')'
                else:
                    stats_file = out_base + outfile_suffix
                    chan_full = channel

                for i, p in enumerate(subs):
                    for v, vis in enumerate(sessions):

                        self._log(f'Subject: {p}, Visit: {vis}', logger)
                        self._log(f'{chan_full}', logger)
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
                            stage=stg,
                            chan_full=[chan_full],
                            reject_epoch=True,
                            reject_artf=self.reject,
                        )
                        evt_target_seg = segments.segments
                        self._log(f'#targets = {len(evt_target_seg)}', logger)
                        for i, seg in enumerate(evt_target_seg):
                            evt_target_seg[i]['start'] = seg['times'][0][0]
                            evt_target_seg[i]['end'] = seg['times'][0][1]
                            evt_target_seg[i]['chan'] = [seg['chan']]
                            evt_target_seg[i]['quality'] = 'Good'

                        # Get events - probe
                        segments = fetch(
                            dset,
                            annot,
                            cat=cat,
                            evt_type=[evttype_probe],
                            stage=stg,
                            chan_full=[chan_full],
                            reject_epoch=True,
                            reject_artf=self.reject,
                        )
                        evt_probe_seg = segments.segments
                        self._log(f'#probes = {len(evt_probe_seg)}', logger)
                        for i, seg in enumerate(evt_probe_seg):
                            evt_probe_seg[i]['start'] = seg['times'][0][0]
                            evt_probe_seg[i]['end'] = seg['times'][0][1]
                            evt_probe_seg[i]['chan'] = [seg['chan']]
                            evt_probe_seg[i]['quality'] = 'Good'

                        # Assess co-occurrence
                        matched = match_events(
                            evt_probe_seg, evt_target_seg, iu_thresh
                        )

                        # Store stats in master table
                        stats[(i * len(sessions)) + v, 0] = matched.recall
                        stats[(i * len(sessions)) + v, 1] = matched.precision
                        stats[(i * len(sessions)) + v, 2] = matched.f1score
                        stats[(i * len(sessions)) + v, 3] = matched.n_tp
                        stats[(i * len(sessions)) + v, 4] = matched.n_fp
                        stats[(i * len(sessions)) + v, 5] = matched.n_fn

                # Export stats table
                self._log(f'Saving {stats_file}', logger)
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
