"""Focused regressions; load functions directly to avoid optional pipeline imports."""
import ast
import logging
from pathlib import Path
from types import SimpleNamespace
import unittest
import numpy as np
from pandas import DataFrame
from scipy.signal import find_peaks, butter, sosfiltfilt

ROOT = Path(__file__).resolve().parents[1] / 'seapipe'
def load_function(module, name, namespace):
    tree = ast.parse((ROOT / module).read_text())
    node = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == name)
    exec(compile(ast.Module(body=[node], type_ignores=[]), str(module), 'exec'), namespace)
    return namespace[name]

logger = logging.getLogger('test')
read = load_function('utils/load.py', 'read_inversion', {
    'DataFrame': DataFrame, 'search_chans': lambda rows: rows.filter(regex='^chan')})
def infer(arrays):
    segments = SimpleNamespace(segments=[{'data': SimpleNamespace(data=[[x]])} for x in arrays],
                               read_data=lambda *args: None)
    def transform(x, fs, method, method_opt):
        high = butter(method_opt['order'], method_opt['freq'][0],
                      btype='highpass', fs=fs, output='sos')
        low = butter(method_opt['order'], method_opt['freq'][1],
                     btype='lowpass', fs=fs, output='sos')
        return sosfiltfilt(low, sosfiltfilt(high, x))
    function = load_function('utils/misc.py', 'infer_polarity', {
        'np': np, 'find_peaks': find_peaks, 'transform_signal': transform,
        'create_logger': lambda *args: logger, 'fetch': lambda *args, **kwargs: segments})
    return function(SimpleNamespace(header={'s_freq': 100}), None, 'C3', [], logger=logger)

class PolarityTests(unittest.TestCase):
    def test_boolean_values(self):
        for value, expected in [('False', False), (' TRUE ', True), (False, False),
                                (True, True), ('0', False), ('1', True), (0.0, False),
                                ('No', False), ('yes', True), ('unknown', None)]:
            with self.subTest(value=value):
                frame = DataFrame([{'sub':'s', 'ses':'v', 'chan':'C3', 'invert':value}])
                self.assertIs(read('s', 'v', frame, 'C3', logger), expected)
    def test_channel_mapping_and_missing_rows(self):
        frame = DataFrame([{'sub':'s', 'ses':'v', 'chan':'C3, C4', 'invert':'False, True'}])
        self.assertIs(read('s','v',frame,'C3',logger), False)
        self.assertIs(read('s','v',frame,'C4',logger), True)
        self.assertIsNone(read('absent','v',frame,'C3',logger))
        self.assertIsNone(read('s','v',frame,'missing',logger))
    def test_waveforms(self):
        t = np.arange(3000) / 100
        symmetric = np.cos(2*np.pi*t)
        positive = symmetric + 0.4*np.cos(4*np.pi*t)
        self.assertIs(infer([positive] * 5), True)
        self.assertIs(infer([-positive] * 5), False)
        self.assertIs(infer([symmetric] * 5), False)
        self.assertIs(infer([positive] * 2), False)
        self.assertIs(infer([positive] * 3 + [-positive] * 2), False)
        self.assertIs(infer([positive * 0] * 5), False)
        self.assertIs(infer([np.full(3000, np.nan)] * 5), False)
        artifact = symmetric.copy(); artifact[1500] = 1000
        self.assertIs(infer([artifact] * 5), False)
        self.assertIs(infer([]), False)

if __name__ == '__main__':
    unittest.main()
