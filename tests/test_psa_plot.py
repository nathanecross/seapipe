"""Real SVG output and export-order regressions without optional pipeline imports."""
import ast
import builtins
import copy
import logging
import warnings
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch
from xml.etree import ElementTree

import numpy as np
from matplotlib.figure import Figure
from wonambi import ChanFreq
from wonambi.trans import export_freq, get_descriptives

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'seapipe/spectrum/psa.py'
LOGGER = logging.getLogger('psa-plot-test')
ns = dict(deepcopy=copy.deepcopy, where=np.where, reshape=np.reshape, asarray=np.asarray,
          get_descriptives=get_descriptives, export_freq=export_freq)
tree = ast.parse(SOURCE.read_text())
nodes = [n for n in tree.body if isinstance(n, ast.FunctionDef)
         and n.name in {'save_plot', '_save_plot', 'save_freq', 'default_general_opts'}]
exec(compile(ast.Module(body=nodes, type_ignores=[]), str(SOURCE), 'exec'), ns)
save_plot, save_freq = ns['save_plot'], ns['save_freq']


def spectrum(multiplier=1, values=None):
    data = ChanFreq()
    data.s_freq = 128
    data.axis['freq'] = np.empty(1, dtype=object)
    data.axis['chan'] = np.empty(1, dtype=object)
    data.data = np.empty(1, dtype=object)
    data.axis['freq'][0] = np.arange(0, 41, 0.25)
    data.axis['chan'][0] = np.array(['EEG'])
    data.data[0] = (multiplier * (1 + np.exp(-((data.axis['freq'][0] - 12) / 2)**2)))[None, :]
    if values is not None:
        data.data[0] = np.asarray(values)[None, :]
    return dict(data=data, start=0, end=20, duration=20, n_stitch=0,
                stage='NREM3', cycle=None, name=None)


class PlotTests(unittest.TestCase):
    def setUp(self):
        self.temp = TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.output = str(Path(self.temp.name)/'sub-01_ses-01_EEG_NREM3**0.25-40Hz.csv')
        self.opts = ns['default_general_opts']()
        self.opts['freq_plot'] = True
        self.frequency = {'frequency': (0.25, 40), 'log_trans': False, 'scaling': 'power'}

    def capture(self, out, **kwargs):
        captured = []
        original = Figure.savefig
        def save(figure, *args, **opts):
            captured.append([dict(x=a.lines[0].get_xdata().copy(),
                                  y=a.lines[0].get_ydata().copy(),
                                  scale=a.get_yscale(), limits=a.get_xlim(), ylabel=a.get_ylabel(),
                                  title=a.get_title()) for a in figure.axes])
            return original(figure, *args, **opts)
        with patch.object(Figure, 'savefig', save):
            result = save_plot(out, self.opts, self.output, LOGGER,
                               frequency_opts=self.frequency, **kwargs)
        return result, captured

    def test_svg_output_defaults_and_no_input_changes(self):
        out = [spectrum()]
        before = copy.deepcopy(out)
        success, captured = self.capture(out)
        self.assertTrue(success)
        path = Path(self.temp.name)/'sub-01_ses-01_EEG_NREM3_powerspectrum_0.25-40Hz.svg'
        root = ElementTree.parse(path).getroot()
        self.assertTrue(root.tag.endswith('svg'))
        self.assertIn('Frequency (Hz)', path.read_text())
        self.assertIn('Power spectral density', path.read_text())
        self.assertEqual(captured[0][0]['scale'], 'log')
        self.assertEqual(captured[0][0]['x'][0], 0.25)
        self.assertEqual(captured[0][0]['x'][-1], 35)
        self.assertEqual(captured[0][0]['limits'], (0.25, 35))
        np.testing.assert_array_equal(out[0]['data'].data[0], before[0]['data'].data[0])
        np.testing.assert_array_equal(out[0]['data'].axis['freq'][0], before[0]['data'].axis['freq'][0])

    def test_plot_is_independent_of_csv_output(self):
        out = [spectrum()]
        success, first = self.capture(out)
        self.assertTrue(success)
        save_freq(out, (0.25, 40), self.opts, self.output, LOGGER)
        success, second = self.capture(out)
        self.assertTrue(success)
        np.testing.assert_array_equal(first[0][0]['x'], second[0][0]['x'])
        np.testing.assert_array_equal(first[0][0]['y'], second[0][0]['y'])
        self.assertEqual(out[0]['data'].axis['freq'][0][-1], 40.75)
        self.assertTrue((Path(self.temp.name)/'sub-01_ses-01_EEG_NREM3_freq_full_0.25-40Hz.csv').exists())

    def test_maximum_and_no_cap(self):
        for maximum, expected in [(20, 20), (None, 39.75), (200, 39.75)]:
            with self.subTest(maximum=maximum):
                self.opts['max_freq_plot'] = maximum
                success, axes = self.capture([spectrum()])
                self.assertTrue(success)
                self.assertEqual(axes[0][0]['x'][-1], expected)

    def test_segments_mean_and_multiple_channels(self):
        out = [spectrum(), spectrum(3)]
        for seg in out:
            seg['data'].axis['chan'][0] = np.array(['EEG', 'Cz'])
            seg['data'].data[0] = np.vstack([seg['data'].data[0], 2*seg['data'].data[0]])
        success, axes = self.capture(out)
        self.assertTrue(success)
        self.assertEqual(len(axes[0]), 2)
        chosen = (out[0]['data'].axis['freq'][0] > 0) & (out[0]['data'].axis['freq'][0] <= 35)
        np.testing.assert_allclose(axes[0][0]['y'], 2*out[0]['data'].data[0][0, chosen])
        np.testing.assert_allclose(axes[0][1]['y'], 2*axes[0][0]['y'])
        self.assertEqual(axes[0][1]['title'], 'Cz')

    def test_log_spectra_are_plotted_on_linear_axis(self):
        self.frequency['log_trans'] = True
        values = np.linspace(-4, 1, 164)
        success, axes = self.capture([spectrum(values=values)])
        self.assertTrue(success)
        self.assertEqual(axes[0][0]['scale'], 'linear')
        self.assertLess(axes[0][0]['y'].min(), 0)
        self.assertEqual(axes[0][0]['ylabel'], 'ln(Power spectral density)')

    def test_normalized_and_laplacian_labels(self):
        for kwargs, label in [({'norm': 'integral'}, 'Relative power spectral density'),
                              ({'norm': 'baseline'}, 'Spectrum / baseline'),
                              ({'laplacian': True}, 'CSD power spectral density')]:
            with self.subTest(kwargs=kwargs):
                success, axes = self.capture([spectrum()], **kwargs)
                self.assertTrue(success)
                self.assertEqual(axes[0][0]['ylabel'], label)

    def test_invalid_samples_leave_gaps(self):
        out = [spectrum()]
        out[0]['data'].data[0][0, 10:13] = [0, -1, np.nan]
        success, axes = self.capture(out)
        self.assertTrue(success)
        self.assertTrue(np.isnan(axes[0][0]['y'][9:12]).all())

    def test_invalid_limits_or_data_logged_without_crashing(self):
        cases = [([], 35, 'No spectra'), ([spectrum()], 0, 'max_freq_plot'),
                 ([spectrum()], np.nan, 'max_freq_plot'),
                 ([spectrum()], 0.01, 'No positive frequency'),
                 ([spectrum(values=np.zeros(164))], 35, 'No finite plottable')]
        different = [spectrum(), spectrum()]
        different[1]['data'].axis['freq'][0] += 0.01
        cases.append((different, 35, 'different channels or frequency grids'))
        for out, maximum, message in cases:
            with self.subTest(message=message):
                self.opts['max_freq_plot'] = maximum
                with self.assertLogs(LOGGER, level='ERROR') as logs:
                    self.assertFalse(save_plot(out, self.opts, self.output, LOGGER,
                                               frequency_opts=self.frequency))
                self.assertIn(message, logs.output[0])

    def test_plot_deprecations_suppressed_and_other_warnings_preserved(self):
        from pyparsing import PyparsingDeprecationWarning
        original = Figure.savefig

        def rendering(figure, *args, **kwargs):
            for module in ['matplotlib._mathtext', 'matplotlib._fontconfig_pattern', 'pyparsing.util']:
                warnings.warn_explicit('plot deprecation', PyparsingDeprecationWarning,
                                       filename='dependency.py', lineno=1, module=module)
            warnings.warn_explicit('numerical warning', RuntimeWarning,
                                   filename='dependency.py', lineno=2, module='matplotlib._mathtext')
            warnings.warn_explicit('other deprecation', DeprecationWarning,
                                   filename='dependency.py', lineno=3, module='matplotlib._mathtext')
            warnings.warn_explicit('application warning', PyparsingDeprecationWarning,
                                   filename='application.py', lineno=4, module='application')
            return original(figure, *args, **kwargs)

        with warnings.catch_warnings(record=True) as recorded:
            warnings.simplefilter('always')
            before = list(warnings.filters)
            with patch.object(Figure, 'savefig', rendering):
                self.assertTrue(save_plot([spectrum()], self.opts, self.output, LOGGER))
            self.assertEqual(warnings.filters, before)
            warnings.warn_explicit('outside plotting', PyparsingDeprecationWarning,
                                   filename='dependency.py', lineno=5, module='matplotlib._mathtext')
        messages = [str(w.message) for w in recorded]
        self.assertNotIn('plot deprecation', messages)
        self.assertIn('numerical warning', messages)
        self.assertIn('other deprecation', messages)
        self.assertIn('application warning', messages)
        self.assertIn('outside plotting', messages)
        # Includes real font/mathtext warnings emitted by this Matplotlib version.
        dependency_warnings = [w for w in recorded if issubclass(w.category, PyparsingDeprecationWarning)
                               and ('matplotlib/' in w.filename or 'pyparsing/' in w.filename)]
        self.assertEqual(dependency_warnings, [])

    def test_warning_filters_restored_when_plotting_fails(self):
        with warnings.catch_warnings():
            warnings.simplefilter('always')
            before = list(warnings.filters)
            with patch.object(Figure, 'savefig', side_effect=OSError('cannot write')):
                with self.assertLogs(LOGGER, level='ERROR'):
                    self.assertFalse(save_plot([spectrum()], self.opts, self.output, LOGGER))
            self.assertEqual(warnings.filters, before)

    def test_missing_matplotlib_is_actionable(self):
        original = builtins.__import__
        def unavailable(name, *args, **kwargs):
            if name.startswith('matplotlib'):
                raise ImportError('matplotlib missing')
            return original(name, *args, **kwargs)
        with patch('builtins.__import__', unavailable), self.assertLogs(LOGGER, level='ERROR') as logs:
            self.assertFalse(save_plot([spectrum()], self.opts, self.output, LOGGER))
        self.assertIn('Install matplotlib', logs.output[0])

    def test_plot_failure_does_not_suppress_csv(self):
        with patch.object(Figure, 'savefig', side_effect=OSError('disk not writable')):
            with self.assertLogs(LOGGER, level='ERROR') as logs:
                self.assertFalse(save_plot([spectrum()], self.opts, self.output, LOGGER))
        self.assertIn('disk not writable', logs.output[0])
        save_freq([spectrum()], (0.25, 40), self.opts, self.output, LOGGER)
        self.assertTrue(list(Path(self.temp.name).glob('*.csv')))

    def test_production_plot_only_exports_for_stage_and_cycle_groups(self):
        block = next(n for n in ast.walk(tree) if isinstance(n, ast.If)
                     and ast.unparse(n.test) == 'not ev_model')
        code = compile(ast.Module(body=[block], type_ignores=[]), str(SOURCE), 'exec')
        for model in ['whole_night', 'per_stage', 'per_cycle', 'stage*cycle']:
            with self.subTest(model=model):
                self.opts.update(freq_full=False, freq_band=False)
                scope = dict(ns, out=spectrum(), ev_model=False, model=model,
                             self=SimpleNamespace(stage=['NREM3'], cycle_idx=[1]), sg=0,
                             outpath=self.temp.name, sub='sub-01', ses='ses-01', fnamechan='EEG',
                             suffix='0.25-40Hz', general_opts=self.opts, freq_bw=(0.25, 40),
                             frequency_opts=self.frequency, logger=LOGGER, norm=None,
                             laplace_flag=False, flag=0)
                scope['save_freq'] = lambda *args: self.fail('CSV output is disabled')
                scope['save_bands'] = lambda *args: self.fail('CSV output is disabled')
                exec(code, scope)
                self.assertEqual(scope['flag'], 0)
        self.assertTrue(list(Path(self.temp.name).glob('*_powerspectrum_*.svg')))
        self.assertFalse(list(Path(self.temp.name).glob('*.csv')))

    def test_both_pipeline_plot_hooks_pass_transform_context(self):
        calls = [n for n in ast.walk(tree) if isinstance(n, ast.Call)
                 and isinstance(n.func, ast.Name) and n.func.id == 'save_plot']
        self.assertEqual(len(calls), 2)
        for call in calls:
            self.assertEqual({k.arg for k in call.keywords}, {'frequency_opts', 'norm', 'laplacian'})


if __name__ == '__main__':
    unittest.main()
