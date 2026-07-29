#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 13:36:12 2023

@author: nathancross
"""
from copy import deepcopy
from datetime import datetime
from json import dump
from os import listdir, mkdir, makedirs, path, walk
from numpy import array, ceil, delete, zeros
from pandas import DataFrame
from pyedflib import highlevel
from shutil import copy2
import re
from statistics import mode
from wonambi import Dataset
from wonambi.attr import Annotations
from .logs import create_logger
from .load import load_channels, load_stages, rename_channels, read_tracking_sheet



def check_dataset(rootpath, datapath, outfile = False, filetype = '.edf',
                  tracking = False, logger = create_logger("Audit")):

    """ Audits the directory specified by <in_dir> to check if the dataset is
        BIDS compatible, how many sessions, recordings (e.g. edfs) and annotations
        files there are per participant.
        You can specify  an optional output filename that will contain the printout.
    """


    if not path.exists(datapath):
        logger.critical(f"PATH: {datapath} does not exist. Check documentation for "
                     "how to arrange data:"
                     "\nhttps://seapipe.readthedocs.io/en/latest/index.html\n")
        return DataFrame()

    logger.debug(f'Checking dataset in directory: {datapath}')

    # Extract participants to check
    if not tracking:
        subs = [x for x in listdir(datapath) if path.isdir(path.join(datapath, x))]
        subs.sort()
    else:
        logger.debug('Reading participant list from tracking sheet.')
        tracking = read_tracking_sheet(rootpath, logger)
        subs = tracking['sub'].drop_duplicates().to_list()
    subs.sort()

    # Initialise certain reporting metrics
    nsd = [] #num subject dirs
    nedf = [] #num subject files
    bids = []
    finalbids = 0
    filesize = 0

    if isinstance(filetype, str):
        filetype = [filetype]

    if len(subs) == 0:
        logger.critical(f"{datapath} doesn't contain any directories.\n")
        finalbids += 1

    for sub in subs:
        real_files = [x for x in listdir(path.join(datapath, sub)) if not x.startswith('.')]
        sessions = [x for x in real_files if path.isdir(path.join(datapath, sub, x))]
        files = [x for x in real_files if path.isfile(path.join(datapath, sub, x))]

        nsd.append(len(sessions))

        edfs = 0

        if len(sessions) < 1:
            finalbids += 1
            if len(files) > 0:
                edfs = len([x for x in files if any(ft in x for ft in filetype)])
                logger.critical(f"{sub} doesn't have sessions directories.\n")
                logger.info('Check documentation for how to setup data in BIDS:')
                logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
                logger.info('-' * 10)
            else:
                logger.critical(f'{sub} has no files\n')
        else:
            for ses in sessions:
                eeg_dir = path.join(datapath, sub, ses, 'eeg')
                if path.exists(eeg_dir):
                    real_files = [x for x in listdir(eeg_dir) if not x.startswith('.')]
                    files = [f for f in real_files if any(ft in f for ft in filetype)]
                    if len(files) == 1:
                        edfs += 1
                        filesize += sum([path.getsize(path.join(eeg_dir, f)) for
                                         f in files if any(ft in f for ft in filetype)])
                    elif len(files) > 1:
                        finalbids += 1
                        logger.critical("BIDS incompatibility. >1 recording file "
                                        f"found for {sub}, {ses}. There should only "
                                        "be 1 recording per session directory\n")
                    else:
                        logger.warning(f'{sub}, {ses} has no files\n')
                else:
                    finalbids += 1
                    logger.critical("BIDS incompatibility. No 'eeg' directory found "
                                    f"for {sub}, {ses}\n")
                    logger.info('Check documentation for how to setup data in BIDS:')
                    logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
                    logger.info('-' * 10)

        bids.append(all([len(dirs2) < 1 for dirs2 in sessions]))
        nedf.append(edfs)

    if len(set(nsd)) > 1:
        logger.warning('Not all participants have the same number of sessions\n')

    subdirs = DataFrame({'BIDS?': bids, '#sessions': nsd, '#recordings': nedf},
                        index=subs)
    subdirs[''] = ['!!' if c1 != c2
                   else '!!' if c1 == 0
                   else '!!' if c2 == 0
                   else ''
                   for c1, c2 in zip(subdirs['#sessions'], subdirs['#recordings'])]

    if outfile:
        if isinstance(outfile, str):
            subdirs.to_csv(outfile)
        elif isinstance(outfile, str):
            subdirs.to_csv(f'{rootpath}/derivatives/seapipe/audit/audit.csv')
        else:
            logger.warning("'outfile' should be set to an instance of boolean or str, not {type(outfile)}. No log will be saved. \n")

    if finalbids == 0:
        logger.info('\n                      Summary:')
        logger.info(f"                      {sum(subdirs['#recordings'])} files, {filesize / (10**9):,.2f} GB")
        logger.info(f"                      Subjects: {subdirs.shape[0]}")
        logger.info(f"                      Sessions: {max(subdirs['#sessions'])}\n")
        logger.debug('The dataset appears compatible for SEAPIPE analysis.\n')
    else:
        logger.critical('The dataset DOES NOT appear compatible for SEAPIPE analysis.\n')

    return subdirs



def make_bids(root_dir, subs = 'all', origin = 'SCN', filetype = '.edf',
              indir = 'sourcedata', logger = create_logger("Make bids")):

    """Copy data from a known source layout into a BIDS-style rawdata tree.

    Files are copied from ``<root_dir>/<indir>`` and written under
    ``<root_dir>/rawdata``. Supported origin values are ``SCN``, ``Woolcock``,
    ``MASS``, and ``SleepProfiler``.
    """

    def _ensure_dir(directory):
        makedirs(directory, exist_ok=True)

    def _copy_file(src, dst):
        if path.exists(dst):
            logger.warning(f"Destination already exists, not overwriting: {dst}")
            return False
        _ensure_dir(path.dirname(dst))
        copy2(src, dst)
        return True

    def _visible_dirs(directory):
        return sorted([x for x in listdir(directory) if not x.startswith('.')
                       if path.isdir(path.join(directory, x))])

    def _visible_files(directory):
        return sorted([x for x in listdir(directory) if not x.startswith('.')
                       if path.isfile(path.join(directory, x))])

    def _prefix(value, prefix):
        value = str(value).strip().rstrip('/')
        return value if value.startswith(prefix) else f'{prefix}{value}'

    def _normalise_sub(value):
        value = str(value).strip()
        match = re.search(r'(\d{4})', value)
        if match is not None and 'ABC' in value:
            return f"sub-ABC{match.group(1)}"
        return _prefix(value, 'sub-')

    def _normalise_ses(value):
        return _prefix(value, 'ses-')

    def _normalise_run(value):
        value = str(value).strip()
        if value.endswith('.0'):
            value = value[:-2]
        return int(value)

    def _session_dirs(src_sub_dir):
        return _visible_dirs(src_sub_dir)

    def _tracking_map(root_dir):
        tracking = read_tracking_sheet(root_dir, logger)
        if isinstance(tracking, str) and tracking == 'error':
            return None
        required_cols = ['sub', 'ses', 'run']
        missing_cols = [x for x in required_cols if x not in tracking.columns]
        if len(missing_cols) > 0:
            logger.critical("Tracking sheet is missing required columns: "
                            f"{', '.join(missing_cols)}")
            return None
        mapping = {}
        for _, row in tracking.iterrows():
            try:
                key = (_normalise_sub(row['sub']), _normalise_run(row['run']))
            except (TypeError, ValueError):
                logger.warning(f"Skipping tracking row with invalid sub/run: {row}")
                continue
            if key in mapping:
                logger.warning(f"Duplicate tracking row for {key[0]}, run {key[1]}. "
                               "Using the first entry.")
                continue
            mapping[key] = _normalise_ses(row['ses'])
        return mapping

    def _sleep_profiler_sub(value):
        match = re.search(r'\b(\d{4})\s+\(', value)
        if match is None:
            return None
        return f"sub-ABC{match.group(1)}"

    def _normalise_sleep_profiler_sub(value):
        match = re.search(r'(\d{4})', str(value))
        if match is None:
            return _normalise_sub(value)
        return f"sub-ABC{match.group(1)}"

    def _sleep_profiler_night(filename):
        night_match = re.search(r'_Export_N(\d+)\.edf$', filename)
        if night_match:
            return int(night_match.group(1))
        if re.search(r'_Export\.edf$', filename):
            return 1
        return None

    root_dir = path.normpath(root_dir)
    if path.basename(root_dir) in ['rawdata', 'DATA']:
        logger.warning("make_bids() received a data directory rather than the "
                       "project root. Using its parent directory as root_dir.")
        root_dir = path.dirname(root_dir)
    in_dir = indir if path.isabs(indir) else path.join(root_dir, indir)
    out_dir = path.join(root_dir, 'rawdata')
    derivs_dir = path.join(root_dir, 'derivatives')
    copied = 0
    skipped = 0
    origin_key = origin.lower()
    _ensure_dir(out_dir)

    if origin_key == 'scn':
        sub_dirs = _visible_dirs(in_dir) if subs == 'all' else subs
        for sub_dir in sub_dirs:
            src_sub_dir = path.join(in_dir, sub_dir)
            if not path.isdir(src_sub_dir):
                logger.warning(f"Subject directory does not exist: {src_sub_dir}")
                skipped += 1
                continue
            sub = _normalise_sub(sub_dir)
            for ses_dir in _session_dirs(src_sub_dir):
                src_ses_dir = path.join(src_sub_dir, ses_dir)
                ses = _normalise_ses(ses_dir)
                eeg_dir = path.join(out_dir, sub, ses, 'eeg')
                stage_dir = path.join(derivs_dir, 'staging_manual', sub, ses)
                for file in _visible_files(src_ses_dir):
                    src = path.join(src_ses_dir, file)
                    if file.endswith(filetype):
                        dst = path.join(eeg_dir, f'{sub}_{ses}_eeg{filetype}')
                    elif file.endswith('.xml'):
                        dst = path.join(stage_dir, f'{sub}_{ses}_eeg.xml')
                    else:
                        continue
                    copied += int(_copy_file(src, dst))

    elif origin_key == 'woolcock':
        sub_dirs = _visible_dirs(in_dir) if subs == 'all' else subs
        for sub_dir in sub_dirs:
            src_sub_dir = path.join(in_dir, sub_dir)
            if not path.isdir(src_sub_dir):
                logger.warning(f"Subject directory does not exist: {src_sub_dir}")
                skipped += 1
                continue
            sub = _normalise_sub(sub_dir)
            for ses_dir in _session_dirs(src_sub_dir):
                src_ses_dir = path.join(src_sub_dir, ses_dir)
                ses = _normalise_ses(ses_dir)
                eeg_dir = path.join(out_dir, sub, ses, 'eeg')
                stage_dir = path.join(derivs_dir, 'staging_manual', sub, ses)
                for file in _visible_files(src_ses_dir):
                    src = path.join(src_ses_dir, file)
                    ext = path.splitext(file)[1]
                    if file.endswith('.xml'):
                        dst = path.join(stage_dir, f'{sub}_{ses}_eeg.xml')
                    elif file.endswith(filetype):
                        dst = path.join(eeg_dir, f'{sub}_{ses}_task-psg_eeg{ext}')
                    else:
                        continue
                    copied += int(_copy_file(src, dst))

    elif origin_key in ['sleepprofiler', 'sleepprofileroriginal']:
        tracking = _tracking_map(root_dir)
        if tracking is None:
            return
        source_dirs = _visible_dirs(in_dir)
        if subs != 'all':
            wanted = [_normalise_sleep_profiler_sub(x) for x in subs]
            source_dirs = [x for x in source_dirs if _sleep_profiler_sub(x) in wanted]
        for sub_dir in source_dirs:
            sub = _sleep_profiler_sub(sub_dir)
            if sub is None:
                logger.warning(f"Could not determine ABC subject ID from {sub_dir}.")
                skipped += 1
                continue
            src_dir = path.join(in_dir, sub_dir)
            edfs = [x for x in _visible_files(src_dir) if x.endswith(filetype)]
            if len(edfs) == 0:
                logger.warning(f"No {filetype} files found for {sub} in {src_dir}.")
                skipped += 1
                continue
            for edf in edfs:
                night = _sleep_profiler_night(edf)
                if night is None:
                    logger.warning(f"Could not determine night/run from {edf}.")
                    skipped += 1
                    continue
                if (sub, night) not in tracking:
                    logger.warning(f"No tracking row found for {sub}, run {night}.")
                    skipped += 1
                    continue
                ses = f"{tracking[(sub, night)]}_{night}"
                src = path.join(src_dir, edf)
                dst = path.join(out_dir, sub, ses, 'eeg',
                                f'{sub}_{ses}_task-sleep_acq-EEG_eeg{filetype}')
                copied += int(_copy_file(src, dst))

    elif origin_key == 'mass':
        dir_check = _visible_dirs(root_dir)
        data_dir = path.join(root_dir, 'Biosignals') if 'Biosignals' in dir_check else root_dir
        annot_dir = path.join(root_dir, 'Annotations') if 'Annotations' in dir_check else None
        stage_dir = path.join(derivs_dir, 'staging')
        _ensure_dir(stage_dir)

        files = [x for x in _visible_files(data_dir) if 'PSG' in x]
        stages = [x for x in _visible_files(data_dir) if 'Base' in x]
        if subs == 'all':
            sublist = []
            for file in files:
                match = re.search(r'SS\d+\s*-\s*(\d+)', file)
                if match is not None:
                    sublist.append(match.group(1))
            sublist = sorted(set(sublist))
        else:
            sublist = [str(x).replace('sub-', '') for x in subs]

        if len(sublist) == 0:
            logger.critical(f'No {filetype} files in {data_dir}. Check paths are correct.')
            return

        stagekey = {'Sleep stage 1' : 1,
                    'Sleep stage 2' : 2,
                    'Sleep stage 3' : 3,
                    'Sleep stage 4' : 3,
                    'Sleep stage ?' : 0,
                    'Sleep stage R' : 4,
                    'Sleep stage W' : 0}

        for sub_id in sublist:
            sub = _normalise_sub(sub_id)
            ses = 'ses-1'
            eeg_dir = path.join(out_dir, sub, ses, 'eeg')
            psg = [x for x in files if sub_id in x]
            if len(psg) == 0:
                logger.warning(f"No PSG file found for {sub}.")
                skipped += 1
                continue
            src = path.join(data_dir, psg[0])
            dst = path.join(eeg_dir, f'{sub}_{ses}_acq-PSG.edf')
            if _copy_file(src, dst):
                copied += 1

            hd = Dataset(dst).header
            s_freq = hd['s_freq']
            dur = hd['n_samples'] / s_freq
            dictionary = {
                "TaskName": "Sleep",
                "SamplingFrequency": s_freq,
                "EEGReference": "Unknown",
                "PowerLineFrequency": "Unknown",
                "SoftwareFilters": "n/a",
                "InstitutionName": "University of Sydney",
                "InstitutionalDepartmentName": "Woolcock Institute of Medical Research",
                "RecordingDuration": dur}
            json_file = '.'.join(dst.split('.')[0:-1]) + '.json'
            with open(json_file, "w") as outfile:
                dump(dictionary, outfile)

            stage_files = [x for x in stages if sub_id in x]
            if len(stage_files) == 0:
                logger.warning(f"No staging EDF found for {sub}.")
            else:
                stage_file = stage_files[0]
                _, _, header = highlevel.read_edf(path.join(data_dir, stage_file))
                epochs = [x for x in header['annotations']]
                length = int(epochs[-1][0])
                hypno = zeros(length)
                for e, epoch in enumerate(epochs):
                    if e == 0:
                        start = 0
                        end = int(epoch[0]) + int(ceil(epoch[1]))
                    else:
                        start = int(epoch[0])
                        end = start + int(ceil(epoch[1]))
                    hypno[start:end] = stagekey[epoch[2]]

                stage_df = DataFrame(columns=['onset', 'duration', 'staging'])
                for row, onset in enumerate(range(0, length, 30)):
                    stage_df.loc[row, 'onset'] = onset
                    stage_df.loc[row, 'duration'] = 30
                    stage_df.loc[row, 'staging'] = int(mode(hypno[onset:onset+30]))
                stage_df.to_csv(path.join(eeg_dir,
                                f'{sub}_{ses}_acq-PSGScoring_events.tsv'),
                                sep='\t', header=True, index=False)

                dst_stage_dir = path.join(stage_dir, sub, ses)
                _copy_file(path.join(data_dir, stage_file),
                           path.join(dst_stage_dir, stage_file))

            if annot_dir is not None:
                annot_files = [x for x in _visible_files(annot_dir) if sub_id in x]
                for afile in annot_files:
                    _copy_file(path.join(annot_dir, afile),
                               path.join(stage_dir, sub, ses, afile))

        load_stages(out_dir, stage_dir, subs=subs)

    else:
        logger.critical(f"Unknown origin '{origin}'. Expected one of: "
                        "SCN, Woolcock, MASS, SleepProfiler.")
        return

    logger.info(f"{origin} BIDS conversion complete: copied {copied} file(s), "
                f"skipped {skipped} file(s).")
    check_dataset(root_dir, out_dir, filetype=filetype)

def extract_channels(in_dir, exclude=None, quality=False):

    """Reads channel information from the files in the directory specified by
    <in_dir> and writes them to the BIDS compatible channels.tsv file per participant
    and session.
    You can specify whether to exclude any channels, if preferrable.
    """

    if not exclude:
        exclude = ['A1','A2','M1','M2']

    parts = [x for x in listdir(in_dir) if '.' not in x]

    for p, part in enumerate(parts):
        ppath = f'{in_dir}/{part}'
        sess = [x for x in listdir(ppath) if '.' not in x]

        for s, ses in enumerate(sess):
            spath = f'{ppath}/{ses}/eeg/'
            files = [x for x in listdir(spath) if '.edf' in x]

            for f, file in enumerate(files):
                src = f'{spath}/{file}'

                data = Dataset(src)
                chans = data.header['orig']['label'] #data.header['chan_name']
                types = array([x.split('-')[0] for x in data.header['orig']['transducer']])
                units = array(data.header['orig']['physical_dim'])

                if exclude:
                    ex = [chans.index(x) for x in exclude if x in chans]
                    chans = delete(array(chans), ex)
                    types = delete(types, ex)
                    units = delete(units, ex)
                else:
                    chans = array(chans)

                # Save dataframe
                df = DataFrame(chans)
                df.columns = ['name']
                df['type'] = types
                df['units'] = units
                df['status'] = 'N/A'
                df['status_description'] = 'N/A'
                df.to_csv(f"{spath}{part}_{ses}_channels.tsv", sep = "\t",
                          header=True, index=False)


def track_processing(self, step, subs, tracking, df, chan, stage, show=False,
                     log=True):

    ## Set up logging
    lg = create_logger('Tracking')

    ## Ensure correct format of chan and stage
    if isinstance(chan, str):
        chan = [chan]
    if isinstance(stage, str):
        stage = [stage]

    ## Track sleep staging
    if 'staging' in step or 'stage' in step:
         stage_df = []
         stage_dict = {}
         spath = self.outpath + '/staging/'
         for sub in subs:
            try:
                stage_ses = next(walk(f'{spath}/{sub}'))[1]
                stage_dict[sub] = dict([(x,[]) if x in stage_ses else (x,'-')
                                     for x in tracking['ses'][sub]])
                stage_df.append([x if x in stage_ses else '-'
                              for x in tracking['ses'][sub]])
            except:
                stage_df.append(['-'])

         # Update tracking
         tracking['staging'] = stage_dict
         df['staging'] = stage_df

         # Check for Artefact or Arousal events
         if list(map(list, list(set(map(tuple, stage_dict.values()))))) == [['-']]:
            lg.warning('Staging has NOT been run.')
         else:
            for sub in stage_dict.keys():
                stage_ses = [x for x in stage_dict[sub].keys()]
                for ses in stage_ses:
                    tracking['staging'][sub][ses] = ['stage']
                    try:
                        xml = [x for x in listdir(f'{spath}/{sub}/{ses}') if '.xml' in x]
                        if len(xml) == 0:
                            if log:
                                lg.warning(f'No staging found for {sub}, {ses}')
                        elif len(xml) > 2:
                            if log:
                                lg.warning(f'>1 staging files found for {sub}, {ses} - only 1 staging file is allowed.')
                        else:
                            xml = xml[0]
                            annot = Annotations(f'{spath}/{sub}/{ses}/{xml}')
                            events = sorted(set([x['name'] for x in annot.get_events()]))
                            for event in events:
                                if event in ['Arou', 'Arousal', 'Artefact']:
                                    tracking['staging'][sub][ses].append(event)
                    except:
                        if log:
                            lg.warning(f'No staging found for {sub}, {ses}')


    ## Track spindle detection
    if 'spindles' in step or 'spindle' in step:
        spin_dict = {}
        spath = self.outpath + '/spindle/'
        df['spindle'] = [['-']] * len(df)
        spin_df = deepcopy(df['spindle'])

        for sub in subs:
           try:
               stage_ses = next(walk(f'{spath}/{sub}'))[1]
               spin_dict[sub] = dict([(x,{}) if x in stage_ses else (x,'-')
                                    for x in tracking['ses'][sub]])
               spin_df.loc[sub] = [x if x in stage_ses else '-'
                             for x in tracking['ses'][sub]]
           except:
               spin_dict[sub] = {'-':'-'}


        # Update tracking
        tracking['spindle'] = spin_dict

        # Check for events
        if list(map(list, list(set(map(tuple, spin_dict.values()))))) == [['-']]:
            if log:
                lg.debug('Spindle detection has NOT been run.')
        else:
            methods = ['Lacourse2018','Moelle2011','Ferrarelli2007','Nir2011',
                       'Wamsley2012','Martin2013','Ray2015','FASST','FASST2',
                       'UCSD','Concordia','Lacourse2018_adap','Moelle2011_adap',
                       'Ferrarelli2007_adap','Nir2011_adap','Wamsley2012_adap',
                       'Martin2013_adap','Ray2015_adap','FASST_adap','FASST2_adap',
                       'UCSD_adap','Concordia_adap']
            for sub in spin_dict.keys():
                for ses in spin_dict[sub]:
                    if not spin_dict[sub][ses] == '-':
                        xml = [x for x in listdir(f'{spath}/{sub}/{ses}') if '.xml' in x]
                        if len(xml) == 0:
                            if log:
                                lg.warning(f'No spindle annotations found for {sub}, {ses}')
                        elif len(xml) > 2:
                            if log:
                                lg.warning(f'>1 spindle annotation files found for {sub}, {ses}..')
                        else:
                            xml = xml[0]
                            try:
                                annot = Annotations(f'{spath}/{sub}/{ses}/{xml}')
                                events = [x for x in annot.get_events() if x['name'] in methods]
                                chans = sorted(set([x['chan'][0] for x in events]))
                                if chan:
                                    chans = [x for x in chans for y in chan if y in x]
                                if len(chans) == 0:
                                    if log:
                                        lg.warning(f'Spindles have NOT been detected for {sub}, {ses}.')
                                    spin_dict[sub][ses] = ('-')
                                    spin_df.loc[sub] = list(map(lambda x: x.replace(ses,'-'),spin_df.loc[sub]))
                                    break
                                else:
                                    for spinchan in chans:
                                        tracking['spindle'][sub][ses][spinchan] = []
                                        methlist = sorted(set([x['name'] for x in events]))
                                        if len(methlist) > 0:
                                            for method in methlist:
                                                update = datetime.fromtimestamp(path.getmtime(f'{spath}/{sub}/{ses}/{xml}')).strftime("%m-%d-%Y, %H:%M:%S")
                                                tracking['spindle'][sub][ses][spinchan].append({'Method':method,
                                                                                     'Stage':'',      # FLAG FOR UPDATE
                                                                                     'Cycle':'',      # FLAG FOR UPDATE
                                                                                     'File':f'{spath}/{sub}/{ses}/{xml}',
                                                                                     'Updated':update})
                            except:
                                lg.warning(f'Error loading nnotations found for {sub}, {ses}')


        df['spindle'] = spin_df


    ## Track slow oscillation detection
    if 'slow wave' in step or 'slow oscillation' in step or 'so' in step:
        so_dict = {}
        spath = self.outpath + '/slow_oscillation/'
        df['slow_osc'] = [['-']] * len(df)
        so_df = deepcopy(df['slow_osc'])

        for sub in subs:
           try:
               stage_ses = next(walk(f'{spath}/{sub}'))[1]
               so_dict[sub] = dict([(x,{}) if x in stage_ses else (x,'-')
                                    for x in tracking['ses'][sub]])
               so_df.loc[sub] = [x if x in stage_ses else '-'
                             for x in tracking['ses'][sub]]
           except:
               so_dict[sub] = {'-':'-'}

        # Update tracking
        tracking['slow_osc'] = so_dict

        # Check for events
        if list(map(list, list(set(map(tuple, so_dict.values()))))) == [['-']]:
            if log:
                lg.debug('Slow oscillation detection has NOT been run.')
        else:
            methods = ['Massimini2004','AASM/Massimini2004','Ngo2015','Staresina2015',]
            for sub in so_dict.keys():
                for ses in so_dict[sub]:
                    if not so_dict[sub][ses] == '-':
                        xml = [x for x in listdir(f'{spath}/{sub}/{ses}') if '.xml' in x]
                        if len(xml) == 0:
                            if log:
                                lg.warning(f'No slow oscillation annotations found for {sub}, {ses}')
                        elif len(xml) > 2:
                            if log:
                                lg.warning(f'>1 slow oscillation annotation files found for {sub}, {ses}..')
                        else:
                            xml = xml[0]
                            try:
                                annot = Annotations(f'{spath}/{sub}/{ses}/{xml}')
                                events = [x for x in annot.get_events() if x['name'] in methods]
                                chans = sorted(set([x['chan'][0] for x in events]))
                                if chan:
                                    chans = [x for x in chans for y in chan if y in x]
                                if len(chans) == 0:
                                    if log:
                                        lg.warning(f'Slow oscillations have NOT been detected for {sub}, {ses}.')
                                    so_dict[sub][ses] = ('-')
                                    so_df.loc[sub] = list(map(lambda x: x.replace(ses,'-'),so_df.loc[sub]))
                                    break
                                else:
                                    for sochan in chans:
                                        tracking['slow_osc'][sub][ses][sochan] = []
                                        methlist = sorted(set([x['name'] for x in events]))
                                        if len(methlist) > 0:
                                            for method in methlist:
                                                update = datetime.fromtimestamp(path.getmtime(f'{spath}/{sub}/{ses}/{xml}')).strftime("%m-%d-%Y, %H:%M:%S")
                                                tracking['slow_osc'][sub][ses][sochan].append({'Method':method,
                                                                                     'Stage':xml.split(f'_{sochan}_')[1].split('_')[0],      # FLAG FOR UPDATE
                                                                                     'Cycle':'',      # FLAG FOR UPDATE
                                                                                     'File':f'{spath}/{sub}/{ses}/{xml}',
                                                                                     'Updated':update})
                            except:
                                lg.warning(f'Error loading Annotations found for {sub}, {ses}')

        df['slow_osc'] = so_df


    ## Track fooof detection
    if 'fooof' in step or 'specparams' in step:
        fooof_dict = {}
        spath = self.outpath + '/fooof/'
        df['fooof'] = [['-']] * len(df)
        fooof_df = deepcopy(df['fooof'])

        for sub in subs:
           try:
               stage_ses = next(walk(f'{spath}/{sub}'))[1]
               fooof_dict[sub] = dict([(x,{}) if x in stage_ses else (x,'-')
                                    for x in tracking['ses'][sub]])
               fooof_df.loc[sub] = [x if x in stage_ses else '-'
                             for x in tracking['ses'][sub]]
           except:
               fooof_dict[sub] = dict([(ses,'-') for ses in tracking['ses'][sub]])


        # Update tracking
        tracking['fooof'] = fooof_dict

        # Check for events
        if list(map(list, list(set(map(tuple, fooof_dict.values()))))) == [['-']]:
            if log:
                lg.debug('FOOOF detection has NOT been run.')
        else:
            for sub in fooof_dict.keys():
                for ses in fooof_dict[sub]:
                    if not fooof_dict[sub][ses] == '-':
                        files = [x for x in listdir(f'{spath}/{sub}/{ses}') if '.csv' in x]

                        chans = sorted(set([file.split(ses)[1].split('_')[1] for file in files]))
                        if chan:
                            chans = [x for x in chans for y in chan if y in x]
                        if len(chans) == 0:
                            if log:
                                lg.warning(f'FOOOF has NOT been run for {sub}, {ses}.')
                            fooof_dict[sub][ses] = ('-')
                            fooof_df.loc[sub] = list(map(lambda x: x.replace(ses,'-'),fooof_df.loc[sub]))
                            break
                        else:
                            for fooofchan in chans:
                                tracking['fooof'][sub][ses][fooofchan] = []
                                chan_files = [file for file in files if f'_{fooofchan}_' in file]
                                for chanfile in chan_files:
                                    update = datetime.fromtimestamp(path.getmtime(f'{spath}/{sub}/{ses}/{chanfile}')).strftime("%m-%d-%Y, %H:%M:%S")
                                    tracking['fooof'][sub][ses][fooofchan].append({'Stage':chanfile.split(f'_{fooofchan}_')[1].split('_')[0],
                                                                              'Cycle':'',      # FLAG FOR UPDATE
                                                                              'Bandwidth':chanfile.split('specparams_')[1].split('.csv')[0],
                                                                              'File':f'{spath}/{sub}/{ses}/{chanfile}',
                                                                              'Updated':update})
        df['fooof'] = fooof_df

    return df, tracking


def check_fooof(self, frequency, chan, ref_chan, stage, cat, cycle_idx, logger):

    bandwidth = f'{frequency[0]}-{frequency[1]}Hz'
    review = []
    for sub in self.tracking['fooof']:
        sessions = list(self.tracking['fooof'][sub].keys())
        for ses in sessions:
            if not self.tracking['fooof'][sub][ses] == '-':
                flag, chanset = load_channels(sub, ses, chan, ref_chan, 0,
                                              logger, verbose=0)
                if flag>0:
                    return 'error', None, None, None
                newchans = rename_channels(sub, ses, chan, logger)

                for c, ch in enumerate(chanset):
                    if newchans:
                        fnamechan = newchans[ch]
                    else:
                        fnamechan = ch
                    try:
                        fooof = self.tracking['fooof'][sub][ses][fnamechan]
                    except:
                        logger.warning(f'No fooof for {sub}, {ses}, {ch}')
                        break
                    if cat[0] + cat[1] == 2: # whole night
                        num_files = 1
                        stagename = '-'.join(stage)
                        files = [x['File'] for x in fooof if stagename in x['Stage']
                                 if bandwidth in x['Bandwidth']]
                    elif cat[0] + cat[1] == 0: # stage*cycle
                        # num_files = len(stage)*len(cycle_idx)
                        # files = []
                        # for stg in stage:
                        #     for cyc in cycle_idx:
                        #         files.append([x['File'] for x in fooof
                        #                       if stage in x['Stage']
                        #                       if cyc in x['Cycle']
                        #                       if bandwidth in x['Bandwidth']])
                        logger.error('Adapted bands for stage*cycle has not yet been implemented')
                        return 'error', None, None, None
                    elif cat[0] == 0:
                        # num_files = len(cycle_idx)
                        # files = []
                        # for cyc in cycle_idx:
                        #     files.append([x['File'] for x in fooof
                        #                   if stage in x['Stage']
                        #                   if cyc in x['Cycle']
                        #                   if bandwidth in x['Bandwidth']])
                        logger.error('Adapted bands for per_cycle has not yet been implemented')
                        return 'error', None, None, None
                    elif cat[1] == 0:
                        num_files = len(stage)
                        files = []
                        for stg in stage:
                            [files.append(x['File']) for x in fooof if stg in x['Stage']
                                          if bandwidth in x['Bandwidth']]

                    if num_files != len(files):
                        flag +=1
            else:
                flag = 1

            if flag>0:
                review.append([sub, ses])

    for row in chan.index:
        subses = [chan['sub'].loc[row], chan['ses'].loc[row]]
        if not subses in review:
            chan = chan.drop([row])

    sub = list(chan['sub'])
    ses = list(chan['ses'])

    return 'review', chan, sub, ses



