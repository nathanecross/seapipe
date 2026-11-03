.. _Power_spectrum:

Power spectral analysis
=======================

.. _psa_overview:

PSA overview
------------

Power spectral analysis (PSA) describes how the power in an EEG signal is
distributed across frequencies. Seapipe uses Wonambi to estimate spectra from
selected sleep stages, cycles or annotation events, and exports frequency-bin
values and integrated power in frequency bands. The default is Welch estimation
with 4-second Hann windows and 50% overlap.

Prepare the dataset as described in :doc:`usage`, including a Wonambi
:ref:`Annotations file <Annotations_file>` with sleep staging. Channel and
reference names can be supplied explicitly or read from the
:ref:`tracking file <tracking file>`. Use recording-channel labels for data
selection; output aliases and montage names have separate roles.

Run power spectral analysis
---------------------------

Create a pipeline, then call ``power_spectrum()``:

.. code-block:: python

    from seapipe import pipeline

    project = pipeline('/path/to/my_project')
    project.power_spectrum(
        subs=['sub-01'],
        sessions=['ses-01'],
        chan=['C3'],
        ref_chan=['M2'],
        stage=['NREM2', 'NREM3'],
        concat_stage=True,
        concat_cycle=True,
    )

This example combines NREM2 and NREM3 across the night. Use the exact recording
labels for ``C3`` and ``M2``. To use channels and references from the tracking
file, leave ``chan`` and ``ref_chan`` as ``None``.

*Function arguments:*

.. code-block:: python

    project.power_spectrum(
        xml_dir=None, out_dir=None,
        subs='all', sessions='all', filetype='.edf',
        chan=None, ref_chan=None, grp_name='eeg', rater=None,
        stage=None, cycle_idx=None,
        concat_cycle=True, concat_stage=False,
        general_opts=None, frequency_opts=None, filter_opts=None,
        epoch_opts=None, event_opts=None,
        norm=None, norm_opts=None, outfile=True,
    )

*Argument descriptions:*

``xml_dir`` / ``out_dir`` (str or None)
    Input annotation directory and output directory, each containing
    ``sub-XXX/ses-XXX`` subdirectories. With ``xml_dir=None``, the pipeline
    searches the derivatives directory for a staging folder. With
    ``out_dir=None``, results are written to
    ``<root_dir>/derivatives/powerspectrum/``.

``subs`` / ``sessions`` (str or list)
    ``'all'`` selects all available subjects or sessions. Supply lists such as
    ``['sub-01']`` and ``['ses-01']`` to select a subset.

``filetype`` (str)
    Recording filename extension, default ``'.edf'``. Other extensions depend
    on the recording formats supported by Wonambi.

``chan`` / ``ref_chan`` (list or None)
    Recording-channel names and reference-channel names, e.g. ``['C3']`` and
    ``['M2']``. ``None`` reads the corresponding tracking-sheet columns.
    ``ref_chan=[]`` requests no rereferencing. Channels are processed separately.

``grp_name`` (str) / ``rater`` (str or None)
    Annotation channel-group name (default ``'eeg'``) and annotation rater.
    ``rater=None`` uses Wonambi's default rater selection.

``stage`` (list or None)
    Sleep stages to analyse. ``None`` becomes
    ``['NREM1', 'NREM2', 'NREM3', 'REM']``. Use ``['Wake']`` for wakefulness.

``cycle_idx`` (list or None)
    One-based indices of annotated sleep cycles, e.g. ``[1, 2, 3]``. Required
    when ``concat_cycle=False``; the annotation file must contain these cycles.

``concat_cycle`` / ``concat_stage`` (bool)
    Combine cycles (default ``True``) and combine selected stages (default
    ``False``). See the grouping table below.

``general_opts`` / ``frequency_opts`` / ``filter_opts`` / ``epoch_opts`` / ``event_opts`` (dict or None)
    Output, spectral-estimation, filtering, segment-selection and event
    settings. ``None`` loads the corresponding default dictionary.

``norm`` (str or None) / ``norm_opts`` (dict or None)
    ``None`` exports unnormalised spectra. ``'integral'`` divides each spectrum
    by its integral over the computed frequency grid, before export-range
    cropping. ``'baseline'`` divides it, frequency by frequency, by an averaged
    baseline spectrum selected using ``norm_opts``. Use matching frequency
    grids for analysis and baseline spectra.

``outfile`` (bool or str)
    ``True`` saves a timestamped log under
    ``<root_dir>/derivatives/audit/logs/``. A filename string selects a custom
    log filename in that directory. ``False`` disables file logging.

Grouping stages and cycles
--------------------------

.. list-table::
   :header-rows: 1
   :widths: 20 20 60

   * - ``concat_cycle``
     - ``concat_stage``
     - Result
   * - ``True``
     - ``True``
     - Combine selected stages and cycles.
   * - ``True``
     - ``False``
     - Separate selected stages, combining cycles (the default).
   * - ``False``
     - ``True``
     - Separate cycles, combining selected stages.
   * - ``False``
     - ``False``
     - Separate each stage and cycle combination.

Discontinuous selected signal is concatenated within each condition. Spectral
windows can therefore cross a join between selected runs. ``concat_stage``
and ``concat_cycle`` determine grouping; ``duration`` and ``overlap`` determine
how spectra are estimated within each resulting segment.

Customise the options
---------------------

Start from the factory functions in ``seapipe.spectrum.psa`` and change only the
values required for the analysis:

.. code-block:: python

    from seapipe.spectrum.psa import (
        default_general_opts, default_frequency_opts, default_filter_opts,
        default_epoch_opts, default_event_opts, default_norm_opts,
    )

    general_opts = default_general_opts()
    frequency_opts = default_frequency_opts()
    filter_opts = default_filter_opts()
    epoch_opts = default_epoch_opts()
    event_opts = default_event_opts()
    norm_opts = default_norm_opts()

    frequency_opts.update(duration=4, overlap=0.5, n_fft=None)
    filter_opts.update(notch_freq=50, highpass=0.25, lowpass=40)
    epoch_opts.update(reject_epoch=True, reject_artf=['Artefact', 'Arousal'])
    event_opts['buffer'] = 0

    project.power_spectrum(
        stage=['NREM3'],
        general_opts=general_opts,
        frequency_opts=frequency_opts,
        filter_opts=filter_opts,
        epoch_opts=epoch_opts,
        event_opts=event_opts,
    )

.. note::

   Nonempty option dictionaries are used as supplied; missing keys are not
   merged with defaults. Use a fresh factory dictionary for each run. The
   pipeline updates some entries, including the export frequency range,
   filename suffix and tracking-derived channel mappings.

The tables below list every key returned by the seven default-option factories
at the top of ``psa.py``. Some settings are shared with PAC or spectral-peak
analyses; these are identified explicitly.

.. _psa_general_opts:

general_opts
~~~~~~~~~~~~

Factory: ``default_general_opts()``.

.. list-table::
   :header-rows: 1
   :widths: 20 25 55

   * - Key
     - Default
     - Meaning and accepted values
   * - ``freq_full``
     - ``True``
     - Export the spectrum as frequency-bin values in a ``_freq_full_`` CSV.
   * - ``freq_band``
     - ``True``
     - Export integrated power for the bands in ``frequency_opts['bands']`` in a ``_freq_band_`` CSV.
   * - ``freq_plot``
     - ``False``
     - Save an SVG power-spectrum figure for each exported analysis group. Independent of the CSV output switches; see :ref:`psa_plotting`.
   * - ``max_freq_plot``
     - ``35``
     - Maximum displayed frequency in Hz (inclusive). Use a positive finite number, or ``None`` for no additional cap beyond the spectrum export range. Does not change CSVs or calculated spectra.
   * - ``suffix``
     - ``None``
     - Output filename suffix. When unset, ``power_spectrum()`` constructs it from the highpass and lowpass values, e.g. ``0.25-40Hz``.
   * - ``chan_grp_name``
     - ``'eeg'``
     - Legacy channel-group setting. Keep this consistent with the public ``grp_name`` argument; this key currently does not alter PSD computation.

.. _psa_frequency_opts:

frequency_opts
~~~~~~~~~~~~~~

Factory: ``default_frequency_opts()``.

.. list-table::
   :header-rows: 1
   :widths: 20 25 55

   * - Key
     - Default
     - Meaning and accepted values
   * - ``bands``
     - ``[(0.25, 1.25), (0.5, 4), (4, 7.75), (8, 11), (11.25, 16), (16.25, 19), (19.25, 35)]``
     - List of frequency-band boundaries in Hz for band-power export. See the default band definitions below.
   * - ``headers``
     - ``['SO', 'Delta', 'Theta', 'Alpha', 'Sigma', 'Low Beta', 'High Beta']``
     - Names corresponding to ``bands``, used by spectral summary routines. Keep the two lists the same length and order.
   * - ``frequency``
     - ``(11, 16)``
     - Spectrum export range in Hz. This is the factory default only: ``power_spectrum()`` overwrites it with ``(filter_opts['highpass'], filter_opts['lowpass'])``.
   * - ``output``
     - ``'spectraldensity'``
     - Wonambi transform output. Use ``'spectraldensity'`` for PSA and its CSV exporters. The underlying transform also accepts ``'complex'`` and ``'csd'`` (cross-spectral density), but these are not supported by this single-channel export workflow.
   * - ``scaling``
     - ``'power'``
     - Spectral scaling: ``'power'`` for power spectral density, or ``'energy'`` for energy spectral density. Wonambi also provides ``'fieldtrip'`` and ``'chronux'`` compatibility scaling.
   * - ``sides``
     - ``'one'``
     - ``'one'`` for nonnegative frequencies or ``'two'`` to include negative frequencies. Use ``'one'`` for ordinary power-spectrum exports.
   * - ``taper``
     - ``'hann'``
     - Window passed to Wonambi/SciPy; common choices are ``'hann'``, ``'boxcar'`` and ``'dpss'`` (multitaper).
   * - ``detrend``
     - ``'linear'``
     - ``'linear'`` removes a linear trend; ``'constant'`` removes the mean; ``None`` disables detrending.
   * - ``n_fft``
     - ``None``
     - FFT length in samples. ``None`` uses the analysis-window sample count; a larger value zero-pads and a smaller value truncates. Set this explicitly to control the frequency grid.
   * - ``fast_fft``
     - ``True``
     - When ``True``, use ``n_fft`` (including its default ``None``). When ``False``, PSA currently uses a fixed FFT length of 2048 samples; spectral-peak detection instead uses the window sample count.
   * - ``duration``
     - ``4``
     - Welch/FFT window duration in seconds. ``None`` transforms the entire selected segment without averaging subwindows. This is distinct from ``epoch_opts['epoch_dur']``.
   * - ``overlap``
     - ``0.5``
     - Fractional overlap between spectral windows, from 0 up to but excluding 1; ``0.5`` means 50%. Used when ``duration`` is not ``None``.
   * - ``step``
     - ``None``
     - Time in seconds between spectral-window starts. A positive value overrides ``overlap``.
   * - ``centend``
     - ``'mean'``
     - Average window spectra with ``'mean'`` or ``'median'``. Used when ``duration`` is not ``None``.
   * - ``log_trans``
     - ``False``
     - Natural-log transform of spectral values before window averaging. Keep ``False`` when calculating ordinary band powers or applying ratio normalization.
   * - ``halfbandwidth``
     - ``3``
     - DPSS half-bandwidth in Hz; applies only to ``taper='dpss'``.
   * - ``NW``
     - ``None``
     - DPSS time-bandwidth product. When supplied, overrides ``halfbandwidth``; the number of tapers is ``2 * NW - 1``.


Default bands and headers, in matching order:

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Header
     - Band (Hz)
   * - SO
     - 0.25 to 1.25
   * - Delta
     - 0.5 to 4
   * - Theta
     - 4 to 7.75
   * - Alpha
     - 8 to 11.25
   * - Sigma
     - 11.25 to 16
   * - Low Beta
     - 16.25 to 19
   * - High Beta
     - 19.25 to 35

The default bands contain overlaps and gaps; replace ``bands`` and ``headers``
together if a different partition is required.

Spectrum export currently looks up the lower and upper cutoffs as exact FFT
bin values. Choose numeric ``highpass`` and ``lowpass`` values that lie on the
frequency grid, even if ``bandpass=False``. For example, 4-second windows with
``n_fft=None`` produce 0.25 Hz spacing. A requested boundary that is absent from
the grid can prevent export. The upper export boundary is excluded.

.. _psa_filter_opts:

filter_opts
~~~~~~~~~~~

Factory: ``default_filter_opts()``.

.. list-table::
   :header-rows: 1
   :widths: 20 25 55

   * - Key
     - Default
     - Meaning and accepted values
   * - ``laplacian``
     - ``False``
     - Apply the spherical-spline surface Laplacian (current source density, CSD) using all loaded ``lapchan`` electrodes. The returned result contains the analysis channel only.
   * - ``oREF``
     - ``None``
     - Montage electrode name for a loaded online-reference channel labelled ``'_REF'``. Used only by Laplacian filtering; has no effect when that channel is absent.
   * - ``lapchan``
     - ``None``
     - List of recording-channel names to load for Laplacian filtering. Include the analysis channel and at least four distinct positioned electrodes for automatic sphere fitting.
   * - ``laplacian_rename``
     - ``False``
     - Enable the montage-to-recording coordinate mapping in ``renames``.
   * - ``renames``
     - ``None``
     - Dictionary mapping montage names to recording labels, e.g. ``{'FpZ': 'EEG'}``. Used only with ``laplacian_rename=True``. Montage matching ignores case and recording labels are preserved.
   * - ``montage``
     - ``'standard_alphabetic'``
     - MNE standard montage name used for Laplacian coordinates. Unknown electrodes, conflicting mappings and duplicate positions raise a logged error.
   * - ``notch``
     - ``True``
     - Apply the fixed-frequency MNE notch filter.
   * - ``notch_freq``
     - ``50``
     - Notch frequency in Hz, or a list of frequencies. Every frequency must lie strictly between zero and Nyquist, and the full filter transition band must also fit below Nyquist.
   * - ``notch_harmonics``
     - ``False``
     - Run MNE ``spectrum_fit`` sinusoid detection/removal. This detects fitted sinusoids rather than applying a fixed list of harmonics.
   * - ``bandpass``
     - ``True``
     - Apply MNE filtering with ``highpass`` and ``lowpass`` before spectral estimation.
   * - ``highpass``
     - ``0.25``
     - High-pass cutoff in Hz and lower spectrum-export boundary. Even when ``bandpass=False``, ``power_spectrum()`` uses this value to set the export range.
   * - ``lowpass``
     - ``40``
     - Low-pass cutoff in Hz and upper spectrum-export boundary. Must exceed ``highpass`` for bandpass analysis and lie below Nyquist.
   * - ``dcomplex``
     - ``'hilbert'``
     - Shared PAC option for analytic-signal extraction (``'hilbert'`` or ``'wavelet'``). Not used by ``power_spectrum()``.
   * - ``filtcycle``
     - ``(3, 6)``
     - Shared PAC filter-cycle setting. Not used by ``power_spectrum()``.
   * - ``width``
     - ``7``
     - Shared PAC wavelet-width setting. Not used by ``power_spectrum()``.


Filtering order is fixed-frequency notch, fitted-sinusoid removal, bandpass,
then optional Laplacian transformation. Temporal filters do not require a
montage or channel renaming. Laplacian filtering requires valid coordinates
for every loaded electrode and does not apply another rereferencing step to
signals already referenced when reading the recording.

For a recording label ``EEG`` that corresponds to montage electrode ``Fpz``:

.. code-block:: python

    filter_opts = default_filter_opts()
    filter_opts.update(
        laplacian=True,
        lapchan=['EEG', 'F3', 'F4', 'C3', 'C4', 'Cz', 'P3', 'P4', 'Pz', 'Oz'],
        laplacian_rename=True,
        renames={'FpZ': 'EEG'},
    )
    project.power_spectrum(chan=['EEG'], ref_chan=[], filter_opts=filter_opts)

Use electrodes that actually exist in the recording. If the configured
Laplacian electrode set cannot be read, the channel is skipped and the reason
is logged. CSD values have different units from ordinary EEG power; interpret
and compare them as Laplacian results.

.. _psa_epoch_opts:

epoch_opts
~~~~~~~~~~

Factory: ``default_epoch_opts()``.

.. list-table::
   :header-rows: 1
   :widths: 20 25 55

   * - Key
     - Default
     - Meaning and accepted values
   * - ``epoch``
     - ``None``
     - Selection mode: ``None`` returns continuous runs within each condition; ``'locked'`` uses scored epochs; ``'unlocked'`` creates windows defined by the epoch settings below.
   * - ``reject_epoch``
     - ``True``
     - Exclude epochs marked Poor quality or staged as Artefact. This setting does not reject selected annotation events.
   * - ``reject_artf``
     - ``True``
     - ``True`` removes events labelled ``'Artefact'``; ``False`` disables event-based artefact rejection. A string or list of strings selects custom artefact-event labels.
   * - ``min_dur``
     - ``0``
     - Minimum selected-segment duration in seconds, before spectral-window estimation.
   * - ``epoch_dur``
     - ``30``
     - Window duration in seconds when ``epoch='unlocked'``. This controls segment selection, not the Welch window length.
   * - ``epoch_overlap``
     - ``0.5``
     - Fractional overlap between unlocked epochs, from 0 up to but excluding 1.
   * - ``epoch_step``
     - ``None``
     - Positive step in seconds between unlocked epochs; overrides ``epoch_overlap``.
   * - ``concat_signal``
     - ``True``
     - Keep ``True`` for non-event PSA; the pipeline resets it to ``True`` otherwise. This key currently does not independently change the discontinuous-signal concatenation passed to Wonambi.

.. _psa_event_opts:

event_opts
~~~~~~~~~~

Factory: ``default_event_opts()``.

.. list-table::
   :header-rows: 1
   :widths: 20 25 55

   * - Key
     - Default
     - Meaning and accepted values
   * - ``evt_type``
     - ``None``
     - ``None`` selects scored data. A list of annotation event names, e.g. ``['spindle']``, selects those events instead.
   * - ``event_chan``
     - ``None``
     - Shared event-channel option. Currently unused by ``power_spectrum()``; this workflow supplies its analysis channel to segment selection.
   * - ``buffer``
     - ``2``
     - Seconds of signal added before and after each selected segment. Set to ``0`` to avoid padding, including when analysing events only.
   * - ``concat_events``
     - ``True``
     - Controls concatenation of event types. ``False`` also selects the per-event export branch; see the event-export limitations below.


To select an event type, set ``event_opts['evt_type']`` to its exact annotation
label and choose the input annotation directory with ``xml_dir``. Leave
``concat_events=True`` for the ordinary export path. ``concat_events=False``
uses an experimental per-event exporter whose current stage/cycle filename
handling can fail; it is not equivalent to a validated event-level table.
The ``event_chan`` key does not override the recording-channel selection.

.. _psa_norm_opts:

norm_opts
~~~~~~~~~

Factory: ``default_norm_opts()``.

.. list-table::
   :header-rows: 1
   :widths: 20 25 55

   * - Key
     - Default
     - Meaning and accepted values
   * - ``norm_cat``
     - ``(1, 1, 1, 0)``
     - Wonambi baseline-concatenation tuple: ``(cycles, stages, discontinuous signal, event types)``. Each element is ``1`` to concatenate or ``0`` to separate.
   * - ``norm_evt_type``
     - ``['norm']``
     - Baseline annotation-event names. With the default, events labelled ``'norm'`` must exist. Set to ``None`` to select scored baseline data instead.
   * - ``norm_stage``
     - ``None``
     - Baseline sleep-stage names, or ``None`` for no stage restriction.
   * - ``norm_epoch``
     - ``None``
     - Baseline epoch-selection mode: ``None``, ``'locked'`` or ``'unlocked'``.


These settings are used only with ``norm='baseline'``. Baseline events/stages
must exist in the selected annotation file. The baseline path applies the
Laplacian when enabled, but currently does not apply the analysis notch and
bandpass stack. Account for this difference when interpreting baseline ratios;
``norm_opts`` does not provide separate baseline filter controls.

.. _psa_plotting:

Plot power spectra
------------------

Set ``general_opts['freq_plot']=True`` to save power-spectrum figures alongside
the CSV outputs. Plotting is disabled by default. It requires Matplotlib,
which can be installed with:

.. code-block:: console

    python -m pip install matplotlib

.. code-block:: python

    from seapipe.spectrum.psa import default_general_opts

    general_opts = default_general_opts()
    general_opts.update(freq_plot=True, max_freq_plot=35)

    project.power_spectrum(
        stage=['NREM3'],
        general_opts=general_opts,
    )

To generate figures without CSVs, also set ``freq_full=False`` and
``freq_band=False`` in ``general_opts``. A figure is named
``<subject>_<session>_<channel>_<condition>_powerspectrum_<suffix>.svg``.
The condition identifies the selected stage(s) and, when applicable, the cycle.
The frequency suffix matches the CSV filename suffix, preventing figures from
different export ranges from sharing a filename.

Each figure contains one panel per channel. If its export group contains
multiple segments, the panel shows their arithmetic mean, with equal weight
per segment, after the requested spectral transformation and normalization.
Groups must have matching channel labels and frequency grids to be averaged.

The horizontal axis shows positive frequencies within the spectrum export
range, up to and including ``max_freq_plot``. The export range's upper boundary
remains excluded. Set ``max_freq_plot=None`` to show all positive frequencies
in that range. This display setting does not change the computed spectrum or
either CSV output.

Ordinary spectra use a logarithmic vertical axis. With
``frequency_opts['log_trans']=True``, already transformed values are shown on a
linear axis, avoiding a second logarithmic transformation. Vertical labels
identify power/energy density, Laplacian CSD, or relative/baseline-normalised
spectra as appropriate; spectral units follow the source signal and scaling.
No additional normalization is applied for plotting.

Nonfinite values are omitted. Zero and negative values are also omitted on a
logarithmic axis; omitted bins appear as gaps. If there are no usable bins, the
frequency grids differ, Matplotlib is unavailable or a file cannot be written,
the reason is logged and that figure is skipped. CSV exports and subsequent
analyses can continue. Figures are written without opening GUI windows, making
this output suitable for batch runs.

Pyparsing deprecation messages from Matplotlib's font/mathtext rendering are
suppressed only while producing these figures. Other warning categories,
analysis logs and errors remain visible, and warning settings are restored
after plotting.

.. _FOOOF analyses:

Spectral-peak detection (FOOOF)
-------------------------------

The remaining factory, ``default_fooof_opts()``, configures
``project.detect_spectral_peaks()``. It is not an argument to
``project.power_spectrum()``. This analysis fits an aperiodic component and
oscillatory peaks to linear-power spectra, then selects peaks within a target
frequency band.

.. code-block:: python

    from seapipe.spectrum.psa import default_fooof_opts

    fooof_opts = default_fooof_opts()
    frequency_opts = default_frequency_opts()
    frequency_opts['duration'] = 4

    project.detect_spectral_peaks(
        stage=['NREM2', 'NREM3'],
        frequency=(9, 16),
        frequency_opts=frequency_opts,
        fooof_opts=fooof_opts,
    )

The public method defaults to NREM2 and NREM3, and writes results under
``<root_dir>/derivatives/fooof/`` unless ``out_dir`` is supplied. Its
``frequency`` argument chooses the peak-selection band, whereas
``fooof_opts['freq_range']`` sets the model-fitting range.

fooof_opts
~~~~~~~~~~

Factory: ``default_fooof_opts()``.

.. list-table::
   :header-rows: 1
   :widths: 20 25 55

   * - Key
     - Default
     - Meaning and accepted values
   * - ``psd_dur``
     - ``2``
     - Legacy spectral-window setting. Currently unused by spectral-peak detection; set ``frequency_opts['duration']`` instead.
   * - ``peak_width_limits``
     - ``[1, 12]``
     - Lower and upper peak bandwidth limits in Hz passed to FOOOF.
   * - ``max_n_peaks``
     - ``8``
     - Maximum number of peaks fitted to a spectrum.
   * - ``min_peak_amplitude``
     - ``0.0``
     - Minimum peak height above the aperiodic fit, in log10 power units.
   * - ``peak_threshold``
     - ``2.0``
     - Relative peak-detection threshold in standard deviations above the flattened spectrum.
   * - ``freq_range``
     - ``[0.5, 35]``
     - Frequency range in Hz over which the aperiodic and peak model is fitted. The input spectrum must cover this range.
   * - ``select_highest``
     - ``True``
     - When several fitted peaks fall within a selected band, select the highest-power peak.
   * - ``thresh_param``
     - ``'PW'``
     - Peak property for the optional selection threshold: ``'PW'`` (power above the aperiodic fit) or ``'BW'`` (bandwidth).
   * - ``bands_fooof``
     - ``None``
     - Bands within which fitted peaks are selected. ``None`` lets ``Spectrum.fooof_it()`` derive bands; the public ``detect_spectral_peaks()`` method overwrites this with its ``frequency`` argument.
   * - ``thresh_select``
     - ``None``
     - Optional minimum value for the peak-selection property chosen by ``thresh_param``. ``None`` disables this additional selection threshold.

PSA output
----------

PSA results are written inside ``out_dir/sub-XXX/ses-XXX/``. Filenames include
the subject, session, channel (or tracking-sheet output alias), selected stage
or stages, and the frequency suffix. Separate-cycle analyses also include the
cycle index.

``*_freq_full_*.csv``
    Frequency-bin values with segment metadata and descriptive statistics.

``*_freq_band_*.csv``
    Integrated powers for the configured frequency bands with segment metadata.

``*_powerspectrum_*.svg``
    Optional spectrum figure, enabled with ``general_opts['freq_plot']=True``.
    See :ref:`psa_plotting` for display limits, averaging and plot-only output.

Both CSV exports are enabled by default; figure export is optional. The public
``power_spectrum()`` method writes files and logs; it does not return an
in-memory spectrum.

Check the logfile for missing annotations, absent channels, rejected segments
and filter failures. Absolute spectral units depend on the signal units read
from the recording, the selected scaling, any Laplacian transformation and
normalization.
