"""Temporal filter regressions using real MNE, without optional pipeline imports."""
import ast
import logging
from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import mne
import numpy as np

SOURCE = Path(__file__).resolve().parents[1] / 'seapipe/utils/misc.py'
if not SOURCE.exists():
    SOURCE = Path(__file__).with_name('misc.py')
namespace = {'mne': mne, 'np': np, 'create_logger': logging.getLogger}
names = {'_prepare_mne_data', '_temporal_filter_mne', 'notch_mne', 'notch_mne2', 'bandpass_mne'}
nodes = [node for node in ast.parse(SOURCE.read_text()).body
         if isinstance(node, ast.FunctionDef) and node.name in names]
exec(compile(ast.Module(body=nodes, type_ignores=[]), str(SOURCE), 'exec'), namespace)
notch = namespace['notch_mne']
harmonics = namespace['notch_mne2']
bandpass = namespace['bandpass_mne']


class TemporalFilterTests(unittest.TestCase):
    def data(self, names=('EEG',), seconds=20):
        t = np.arange(int(256 * seconds)) / 256
        signal = np.sin(2 * np.pi * 10 * t) + np.sin(2 * np.pi * 50 * t)
        return SimpleNamespace(axis={'chan': [np.array(names)]}, s_freq=256,
                               data=[np.array([signal * (i + 1) for i in range(len(names))])])

    def assert_fails(self, data, channel, freq, message):
        original = data.data[0].copy()
        with self.assertLogs('Notch filter', level='ERROR') as logs:
            with self.assertRaisesRegex(ValueError, message):
                notch(data, channel, freq)
        self.assertIn('channel=', logs.output[0])
        np.testing.assert_array_equal(data.data[0], original)

    def test_recording_labels_need_no_montage(self):
        for label in ['EEG', 'FpZ', 'Fpz', '_REF', 'custom lead']:
            with self.subTest(label=label):
                data = self.data((label,))
                original = data.data[0].copy()
                with patch.object(mne.channels, 'make_standard_montage',
                                  side_effect=AssertionError('Montage is unnecessary')):
                    result, flag = notch(data, label, 50)
                self.assertEqual(flag, 0)
                np.testing.assert_array_equal(data.data[0], original)
                self.assertEqual(data.axis['chan'][0].tolist(), [label])
                # Ignore FIR edge transients; preserve 10 Hz and attenuate 50 Hz.
                t = np.arange(result.shape[1]) / 256
                middle = slice(1280, -1280)
                amplitude = lambda hz: abs(2 * np.mean(result[0, middle] *
                                            np.exp(-2j * np.pi * hz * t[middle])))
                self.assertGreater(amplitude(10), 0.9)
                self.assertLess(amplitude(50), 0.05)

    def test_selection_and_order(self):
        data = self.data(('EEG', 'eeg', 'FpZ'))
        all_rows, _ = notch(data, ['EEG', 'eeg', 'FpZ'], [50, 100])
        reordered, _ = notch(data, ['FpZ', 'EEG'], [50, 100])
        np.testing.assert_allclose(reordered, all_rows[[2, 0]])
        only_eeg, _ = notch(data, 'eeg', 50)
        self.assertEqual(only_eeg.shape, (1, 5120))

    def test_bad_selection(self):
        for selection in ['Fpz', [], None, ['EEG', 'EEG'], [0]]:
            with self.subTest(selection=selection):
                self.assert_fails(self.data(), selection, 50,
                                  'absent|nonempty|iterable|Duplicate')

    def test_bad_frequencies(self):
        for freq in [None, [], 0, -50, 128, 200, np.nan, np.inf, [[50]], [50, np.nan]]:
            with self.subTest(freq=freq):
                self.assert_fails(self.data(), 'EEG', freq, 'Notch frequencies')

    def test_bad_sampling_rate(self):
        for rate in [0, -256, np.nan, np.inf]:
            data = self.data()
            data.s_freq = rate
            self.assert_fails(data, 'EEG', 50, 'Sampling frequency')

    def test_nonfinite_or_complex_samples(self):
        for value in [np.nan, np.inf]:
            data = self.data()
            data.data[0][0, 12] = value
            self.assert_fails(data, 'EEG', 50, 'finite, real')
        data = self.data()
        data.data[0] = data.data[0].astype(complex)
        self.assert_fails(data, 'EEG', 50, 'finite, real')

    def test_bad_shape(self):
        for values in [np.zeros((2, 10)), np.zeros((1, 0)), np.zeros(10)]:
            data = self.data()
            data.data[0] = values
            self.assert_fails(data, 'EEG', 50, 'channels-by-samples')

    def test_duplicate_labels(self):
        self.assert_fails(self.data(('EEG', 'EEG')), 'EEG', 50, 'Duplicate channel')

    def test_mne_failure_has_context_and_preserves_input(self):
        data = self.data()
        with patch.object(mne.io.RawArray, 'notch_filter',
                          side_effect=RuntimeError('backend failed')):
            self.assert_fails(data, 'EEG', 50, 'backend failed')

    def test_other_temporal_filters_without_coordinates(self):
        data = self.data(('EEG', 'FpZ'))
        original = data.data[0].copy()
        with patch.object(mne.channels, 'make_standard_montage',
                          side_effect=AssertionError('Montage is unnecessary')):
            for function, kwargs in [(harmonics, {}),
                                     (bandpass, {'highpass': 1, 'lowpass': 30})]:
                filtered, flag = function(data, ['FpZ', 'EEG'], **kwargs)
                self.assertEqual(flag, 0)
                self.assertEqual(filtered.shape, original.shape)
                self.assertTrue(np.isfinite(filtered).all())
                np.testing.assert_array_equal(data.data[0], original)

    def test_invalid_bandpass_logged(self):
        with self.assertLogs('Bandpass filter', level='ERROR'):
            with self.assertRaisesRegex(ValueError, 'MNE filter failed'):
                bandpass(self.data(), 'EEG', 1, 150)


if __name__ == '__main__':
    unittest.main()
