"""Regression tests isolated from optional seapipe pipeline dependencies."""
import ast
import json
import logging
import tempfile
import unittest
from unittest.mock import patch
from pathlib import Path
from datetime import datetime
from os import path, makedirs, walk, listdir, mkdir
from shutil import copy2
import re
import numpy as np
from pandas import DataFrame, read_csv
from pyedflib import highlevel

AUDIT = Path(__file__).resolve().parents[1] / 'seapipe/utils/audit.py'
if not AUDIT.exists():
    AUDIT = Path(__file__).with_name('audit.py')
namespace = dict(path=path, makedirs=makedirs, walk=walk, copy2=copy2,
                 re=re, DataFrame=DataFrame, dump=json.dump)
tree = ast.parse(AUDIT.read_text())
names = {'_ensure_dir', '_copy_file', '_prefix', '_sleep_edfx_copy', '_sleep_edfx_write',
         '_sleep_edfx_metadata', '_make_bids_sleep_edfx'}
exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n, ast.FunctionDef)
                              and n.name in names], type_ignores=[]), str(AUDIT), 'exec'), namespace)
convert = namespace['_make_bids_sleep_edfx']
logger = logging.getLogger('sleep-edfx-test')
load_namespace = dict(path=path, listdir=listdir, mkdir=mkdir, read_csv=read_csv,
                      create_logger=lambda *args: logger)
load_tree = ast.parse(AUDIT.with_name('load.py').read_text())
exec(compile(ast.Module(body=[n for n in load_tree.body if isinstance(n, ast.FunctionDef)
                              and n.name in {'load_stages', '_stage_interval_epochs'}],
                        type_ignores=[]), 'load.py', 'exec'), load_namespace)

class SleepEDFTests(unittest.TestCase):
    def test_telemetry_pairing_timing_and_rerun(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); src = root/'sleep-telemetry'; src.mkdir()
            out = root/'rawdata'; out.mkdir()
            signals = [np.zeros(1000), np.zeros(1000), np.zeros(10)]
            headers = [highlevel.make_signal_header('EEG Fpz-Cz', sample_frequency=100),
                       highlevel.make_signal_header('EEG Pz-Oz', sample_frequency=100),
                       highlevel.make_signal_header('Event marker', sample_frequency=1)]
            header = highlevel.make_header(startdate=datetime(2000, 1, 1))
            highlevel.write_edf(str(src/'ST7011J0-PSG.edf'), signals, headers, header)
            annot = highlevel.make_header(startdate=datetime(2000, 1, 1, 0, 0, 2))
            annot['annotations'] = [[0, 3, 'Sleep stage 4'], [3, 5, 'Movement time']]
            highlevel.write_edf(str(src/'ST7011JP-Hypnogram.edf'), signals, headers, annot)
            self.assertEqual(convert(str(root), str(out), 'all', logger), (1, 0))
            stem = out/'sub-ST701/ses-X1/eeg/sub-ST701_ses-X1_task-sleep_acq-PSG'
            events = read_csv(str(stem)+'_events.tsv', sep='\t')
            self.assertEqual(events.onset.tolist(), [2, 5])
            self.assertEqual(events.duration.tolist(), [3, 5])
            self.assertEqual(events.trial_type.tolist(), ['Sleep stage 4', 'Movement time'])
            channels = read_csv(str(stem)+'_channels.tsv', sep='\t')
            self.assertEqual(channels.sampling_frequency.tolist(), [100, 100, 1])
            self.assertEqual(channels.reference.fillna('n/a').tolist(), ['Cz', 'Oz', 'n/a'])
            edf = Path(str(stem)+'_eeg.edf')
            self.assertEqual(edf.read_bytes(), (src/'ST7011J0-PSG.edf').read_bytes())
            modified = edf.stat().st_mtime_ns
            self.assertEqual(convert(str(root), str(out), 'all', logger), (0, 1))
            self.assertEqual(edf.stat().st_mtime_ns, modified)
            # An old interrupted copy must stop before missing sidecars are created.
            edf.write_bytes(b'truncated')
            metadata = Path(str(stem)+'_eeg.json'); metadata.unlink()
            with self.assertRaisesRegex(FileExistsError, 'differs from source'):
                convert(str(src), str(out), 'all', logger)
            self.assertFalse(metadata.exists())
            (src/'ST7011JP-Hypnogram.edf').unlink()
            with self.assertRaisesRegex(ValueError, 'expected one PSG'):
                convert(str(src), str(out), 'all', logger)

class CopyRecoveryTests(unittest.TestCase):
    def test_failed_copy_is_not_published_and_retry_succeeds(self):
        with tempfile.TemporaryDirectory() as tmp:
            src, dst = Path(tmp)/'source.edf', Path(tmp)/'out/data.edf'
            src.write_bytes(b'complete recording')
            def interrupted(source, destination):
                Path(destination).write_bytes(b'partial')
                raise OSError('copy interrupted')
            with patch.dict(namespace, copy2=interrupted):
                with self.assertRaisesRegex(OSError, 'copy interrupted'):
                    namespace['_sleep_edfx_copy'](str(src), str(dst), logger)
            self.assertFalse(dst.exists())
            self.assertEqual(list(dst.parent.iterdir()), [])
            self.assertTrue(namespace['_sleep_edfx_copy'](str(src), str(dst), logger))
            self.assertEqual(dst.read_bytes(), src.read_bytes())
            self.assertFalse(namespace['_sleep_edfx_copy'](str(src), str(dst), logger))

    def test_existing_truncated_and_same_size_corrupt_copies_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            src, dst = Path(tmp)/'source.edf', Path(tmp)/'data.edf'
            src.write_bytes(b'complete recording')
            for content in (b'partial', b'x' * src.stat().st_size):
                dst.write_bytes(content)
                with self.assertRaisesRegex(FileExistsError, 'differs from source'):
                    namespace['_sleep_edfx_copy'](str(src), str(dst), logger)
                self.assertEqual(dst.read_bytes(), content)

    def test_failed_sidecar_write_does_not_leave_partial_json(self):
        with tempfile.TemporaryDirectory() as tmp:
            dst = Path(tmp)/'sidecar.json'
            with self.assertRaises(ValueError):
                namespace['_sleep_edfx_write'](str(dst), {'invalid': float('nan')}, logger)
            self.assertEqual(list(Path(tmp).iterdir()), [])
            namespace['_sleep_edfx_write'](str(dst), {'TaskName': 'sleep'}, logger)
            self.assertEqual(json.loads(dst.read_text()), {'TaskName': 'sleep'})


class StageLoaderTests(unittest.TestCase):
    def test_interval_gaps_mixed_epochs_and_recording_boundaries(self):
        frame = DataFrame({'onset': [-30, 60, 120, 135, 160, 180],
                           'duration': [90, 30, 15, 15, 20, 60],
                           'trial_type': ['W', 'N2', 'N2', 'R', 'W', 'R']})
        epochs = [{'start': x, 'end': x + 30} for x in range(0, 210, 30)]
        stages, uncertain = load_namespace['_stage_interval_epochs'](
            frame, epochs, 200, {'W': 'Wake', 'N2': 'NREM2', 'R': 'REM'})
        self.assertEqual([s for _, s in stages],
                         ['Wake', 'Wake', 'NREM2', 'Undefined', 'Undefined', 'Undefined', 'REM'])
        self.assertEqual(uncertain, 3)
        for start, duration in [(float('nan'), 30), (0, -1), (0, float('inf'))]:
            with self.assertRaises(ValueError):
                load_namespace['_stage_interval_epochs'](
                    DataFrame({'onset': [start], 'duration': [duration], 'trial_type': ['W']}),
                    epochs, 200, {})
        with self.assertRaisesRegex(ValueError, 'Overlapping'):
            load_namespace['_stage_interval_epochs'](
                DataFrame({'onset': [0, 30], 'duration': [60, 60], 'trial_type': ['W', 'W']}),
                epochs, 200, {})

    def test_bids_intervals_and_legacy_stages_through_wonambi(self):
        from wonambi import Dataset
        from wonambi.attr import Annotations, create_empty_annotations
        load_namespace.update(Dataset=Dataset, Annotations=Annotations,
                              create_empty_annotations=create_empty_annotations)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); raw = root/'rawdata'; xml = root/'xml'; xml.mkdir()
            eeg = raw/'sub-ST701/ses-X1/eeg'; eeg.mkdir(parents=True)
            stem = eeg/'sub-ST701_ses-X1_task-sleep_acq-PSG'
            highlevel.write_edf(str(stem)+'_eeg.edf', [np.zeros(27300)],
                               [highlevel.make_signal_header('EEG Fpz-Cz', sample_frequency=100)],
                               highlevel.make_header(startdate=datetime(2000, 1, 1)))
            events = Path(str(stem)+'_events.tsv')
            DataFrame({'onset': [0, 60, 90, 120, 150, 180, 210, 240, 270],
                       'duration': [60, 30, 30, 30, 30, 30, 30, 30, 60],
                       'trial_type': ['Sleep stage W', 'Sleep stage 1', 'Sleep stage 2',
                                      'Sleep stage 3', 'Sleep stage 4', 'Sleep stage R',
                                      'Movement time', 'Sleep stage ?', 'Sleep stage W']}).to_csv(
                                          events, sep='\t', index=False)
            self.assertEqual(load_namespace['load_stages'](str(raw), str(xml), logger=logger), 0)
            xmlfile = xml/'sub-ST701/ses-X1/sub-ST701_ses-X1_staging.xml'
            self.assertEqual([e['stage'] for e in Annotations(str(xmlfile)).get_epochs()],
                             ['Wake', 'Wake', 'NREM1', 'NREM2', 'NREM3', 'NREM3',
                              'REM', 'Movement', 'Undefined', 'Wake'])
            # A rerun must clear stages where the source now contains a gap.
            frame = read_csv(events, sep='\t'); frame = frame[frame.onset != 90]
            frame.to_csv(events, sep='\t', index=False)
            load_namespace['load_stages'](str(raw), str(xml), logger=logger)
            self.assertEqual(Annotations(str(xmlfile)).get_epochs()[3]['stage'], 'Undefined')
            events.unlink()
            legacy = Path(str(stem)+'Scoring_events.tsv')
            DataFrame({'onset': list(range(0, 300, 30)), 'staging': [0, 1, 2, 3, 4, 9, 0, 1, 2, 3]}).to_csv(
                legacy, sep='\t', index=False)
            load_namespace['load_stages'](str(raw), str(xml), logger=logger)
            self.assertEqual([e['stage'] for e in Annotations(str(xmlfile)).get_epochs()],
                             ['Wake', 'NREM1', 'NREM2', 'NREM3', 'REM', 'Undefined',
                              'Wake', 'NREM1', 'NREM2', 'NREM3'])


if __name__ == '__main__':
    unittest.main()
