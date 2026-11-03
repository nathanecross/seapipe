"""CSD numerical and pipeline regressions without importing optional seapipe modules."""
import ast
import copy
import inspect
import logging
from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import Mock
import mne
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
LOGGER = logging.getLogger('laplacian-test')
NS = {'mne': mne, 'np': np, 'create_logger': logging.getLogger}
NAMES = {'_prepare_mne_data', '_laplacian_montage', 'laplacian_mne',
         '_temporal_filter_mne', 'notch_mne', 'notch_mne2', 'bandpass_mne'}
tree = ast.parse((ROOT / 'seapipe/utils/misc.py').read_text())
exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n, ast.FunctionDef)
                             and n.name in NAMES], type_ignores=[]), 'misc.py', 'exec'), NS)
laplacian = NS['laplacian_mne']
CHANNELS = ['Fpz', 'F3', 'F4', 'C3', 'C4', 'Cz', 'P3', 'P4', 'Pz', 'Oz']
PIPELINES = ['seapipe/spectrum/psa.py', 'seapipe/pac/pacats.py', 'seapipe/pac/octopus.py']


def recording(names=None):
    names = CHANNELS if names is None else names
    return SimpleNamespace(s_freq=256, axis={'chan': [np.array(names)]},
                           data=[np.random.RandomState(12).normal(size=(len(names), 1280))])


def compile_body(body, namespace, loop=False):
    """Execute production control flow with the surrounding channel loop intact."""
    header = ast.parse('def run():\n    flag = 0\n    for ch in chanset:\n        pass\n    return True')
    if loop:
        header.body[0].body[1].body = copy.deepcopy(body)
    else:
        header.body[0].body = copy.deepcopy(body)
    ast.fix_missing_locations(header)
    exec(compile(header, 'pipeline-regression', 'exec'), namespace)
    return namespace['run']


class LaplacianTests(unittest.TestCase):
    def test_matches_mne_and_keeps_input(self):
        data = recording()
        before = data.data[0].copy()
        raw = mne.io.RawArray(before.copy(), mne.create_info(CHANNELS, 256, 'eeg'), verbose=40)
        raw.set_montage('standard_alphabetic', verbose=40)
        expected = mne.preprocessing.compute_current_source_density(raw, verbose=40).get_data()
        actual, flag = laplacian(data, ['Oz', 'Fpz'], ref_chan=['unloaded M1'])
        self.assertEqual(flag, 0)
        np.testing.assert_allclose(actual, expected[[9, 0]], rtol=1e-10, atol=1e-8)
        np.testing.assert_array_equal(data.data[0], before)
        self.assertEqual(data.axis['chan'][0].tolist(), CHANNELS)

    def test_case_alias_and_online_reference(self):
        expected, _ = laplacian(recording(), ['Fpz', 'Cz'])
        for name, kwargs in [('FpZ', {}), ('EEG', {'laplacian_rename': True,
                                                  'renames': {'FpZ': 'EEG'}})]:
            names = [name] + CHANNELS[1:]
            names[5] = '_REF'
            data = recording(names)
            result, _ = laplacian(data, [name, '_REF'], oREF='cZ', **kwargs)
            np.testing.assert_allclose(result, expected, rtol=1e-10, atol=1e-8)
        # An oREF mapping has no effect when the _REF channel was not loaded.
        actual, _ = laplacian(recording(), ['Fpz', 'Cz'], oREF='not used')
        np.testing.assert_allclose(actual, expected)

    def test_bad_coordinates_and_mappings_are_logged(self):
        cases = [
            (['unknown'] + CHANNELS[1:], {}, 'Montage label'),
            (['EEG'] + CHANNELS[1:], {'laplacian_rename': True}, 'nonempty renames'),
            (['EEG'] + CHANNELS[1:], {'laplacian_rename': True,
                                     'renames': {'EEG': 'Fpz'}}, 'Recording label'),
            (['FpZ'] + CHANNELS[1:], {'laplacian_rename': True,
                                     'renames': {'Fpz': 'F3'}}, 'share a montage'),
            (['_REF'] + CHANNELS[1:], {'oREF': 'bogus'}, 'Online reference'),
            (['_REF'] + CHANNELS[1:], {'oREF': 'Fpz', 'laplacian_rename': True,
                                      'renames': {'F3': '_REF'}}, 'different positions'),
            (CHANNELS[:3], {}, 'at least four'),
        ]
        for names, kwargs, message in cases:
            with self.subTest(kwargs=kwargs, names=names):
                data = recording(names)
                before = data.data[0].copy()
                with self.assertLogs('Laplacian filter', level='ERROR') as logs:
                    with self.assertRaisesRegex(ValueError, message):
                        laplacian(data, names[0], **kwargs)
                self.assertIn('montage=', logs.output[0])
                np.testing.assert_array_equal(data.data[0], before)

    def test_bandpower_stack_has_consistent_rows_and_labels(self):
        tree = ast.parse((ROOT / 'seapipe/spectrum/bandpower.py').read_text())
        fn = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == '_apply_filters')
        ns = dict(NS, asarray=np.asarray)
        exec(compile(ast.Module(body=[fn], type_ignores=[]), 'bandpower.py', 'exec'), ns)
        opts = {'notch': True, 'notch_freq': 50, 'notch_harmonics': False,
                'bandpass': True, 'highpass': 2, 'lowpass': 30, 'oREF': None,
                'laplacian_rename': False, 'renames': None, 'montage': 'standard_alphabetic'}
        data = recording()
        expected = recording()
        expected.data[0], _ = NS['notch_mne'](expected, CHANNELS, 50)
        expected.data[0], _ = NS['bandpass_mne'](expected, CHANNELS, 2, 30)
        result, _ = laplacian(expected, 'F3')
        self.assertEqual(ns['_apply_filters'](data, CHANNELS, ['M1'], opts, True, 'F3', LOGGER), 0)
        self.assertEqual(data.axis['chan'][0].tolist(), ['F3'])
        self.assertEqual(data.data[0].shape, (1, 1280))
        np.testing.assert_allclose(data.data[0], result)

    def test_psa_and_pac_csd_outputs_have_one_matching_label(self):
        for name in PIPELINES:
            tree = ast.parse((ROOT / name).read_text())
            block = next(n for n in ast.walk(tree) if isinstance(n, ast.If)
                         and ast.unparse(n.test) == "filter_opts['laplacian'] and laplace_flag")
            data = recording()
            expected, _ = laplacian(recording(), 'F3')
            ns = dict(NS, data=data, ch='F3', chanset={'F3': ['M1'], 'F4': []},
                      laplace_flag=True, logger=LOGGER, asarray=np.asarray,
                      filter_opts={'laplacian': True, 'oREF': None,
                                   'laplacian_rename': False, 'renames': None,
                                   'montage': 'standard_alphabetic'})
            with self.subTest(name=name):
                body = [block]
                if '/pac/' in name:
                    body += ast.parse('return dat').body
                result = compile_body(body, ns)()
                np.testing.assert_allclose(data.data[0], expected)
                self.assertEqual(data.axis['chan'][0].tolist(), ['F3'])
                if '/pac/' in name:
                    np.testing.assert_allclose(result, expected[0])

    def test_all_production_calls_match_new_signatures(self):
        for path in (ROOT / 'seapipe').rglob('*.py'):
            for node in ast.walk(ast.parse(path.read_text())):
                if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id in NAMES:
                    if node.func.id.startswith('_'): continue
                    with self.subTest(path=path, line=node.lineno):
                        inspect.signature(NS[node.func.id]).bind(*[None] * len(node.args),
                                                               **{kw.arg: None for kw in node.keywords})

    def test_missing_channels_skip_channel_without_fallback(self):
        for name in PIPELINES:
            tree = ast.parse((ROOT / name).read_text())
            # Select production read block followed by its channel-loop continuation.
            candidate = next(n for n in ast.walk(tree) if isinstance(n, ast.Try)
                             and n.body and isinstance(n.body[0], ast.Expr)
                             and ast.unparse(n.body[0]).startswith('segments.read_data(readchans'))
            read = Mock(side_effect=ValueError('missing F3'))
            ns = {'segments': SimpleNamespace(read_data=read), 'readchans': ['F3'],
                  'chanset': {'EEG': []}, 'sub': 'sub-test', 'ses': 'ses-test',
                  'flag': 0, 'logger': LOGGER}
            with self.subTest(name=name), self.assertLogs(LOGGER, level='ERROR') as logs:
                run = compile_body([candidate, ast.parse('raise AssertionError("channel not skipped")').body[0]], ns, loop=True)
                self.assertTrue(run())
            self.assertEqual(read.call_count, 1)
            self.assertIn('missing F3', logs.output[0])

    def test_psa_baseline_with_and_without_laplacian(self):
        tree = ast.parse((ROOT / 'seapipe/spectrum/psa.py').read_text())
        baseline = next(n for n in ast.walk(tree) if isinstance(n, ast.If)
                        and ast.unparse(n.test) == "norm == 'baseline'")
        class ReachedSpectrum(Exception): pass
        for enabled in [True, False]:
            data = recording()
            class Segments(list):
                def read_data(self, channels, ref):
                    if not enabled:
                        data.data[0] = data.data[0][[1]]
                        data.axis['chan'][0] = np.asarray(['F3'])
            segments = Segments([{'data': data}])
            def frequency(actual, **kwargs):
                self.assertEqual(actual.data[0].shape, (1, 1280))
                self.assertEqual(actual.axis['chan'][0].tolist(), ['F3'])
                raise ReachedSpectrum
            opts = {'laplacian': enabled, 'lapchan': CHANNELS, 'oREF': None,
                    'laplacian_rename': False, 'renames': None,
                    'montage': 'standard_alphabetic'}
            norm_opts = dict.fromkeys(['norm_cat', 'norm_evt_type', 'norm_stage', 'norm_epoch'])
            freq_opts = dict.fromkeys(['output', 'scaling', 'sides', 'taper', 'halfbandwidth', 'NW',
                                      'duration', 'overlap', 'step', 'detrend', 'n_fft', 'log_trans', 'centend'])
            ns = dict(NS, logger=LOGGER, norm='baseline', norm_opts=norm_opts,
                      filter_opts=opts, frequency_opts=freq_opts, dset=None, annot=None,
                      chanset={'F3': ['M1']}, fetch=lambda *a, **k: segments,
                      frequency=frequency, asarray=np.asarray)
            with self.subTest(enabled=enabled), self.assertRaises(ReachedSpectrum):
                compile_body([baseline], ns, loop=True)()


if __name__ == '__main__':
    unittest.main()
