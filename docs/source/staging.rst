Staging
=====

.. _Overview:

Overview
------------
Sleep staging is the process of dividing a sleep recording into standard physiological states, 
typically Wake, NREM1, NREM2, NREM3, and REM, across consecutive epochs. In most PSG-style workflows, 
the recording is scored in 30-second windows using EEG together with EOG and EMG, because each stage has 
characteristic patterns such as alpha attenuation at sleep onset, spindles and K-complexes in NREM2, 
high-amplitude slow waves in NREM3, and rapid eye movements with reduced muscle tone in REM.

Traditionally, sleep staging is done manually by a trained scorer who reviews the recording epoch by epoch 
and assigns stages according to published criteria such as AASM rules. Manual scoring is still the reference 
standard because it allows contextual judgment when signals are noisy or borderline, but it is slow, 
labour-intensive, and introduces scorer variability. Automatic staging uses algorithms or machine-learning models 
to assign stages from the same signals. These methods are much faster and scale well to large datasets, but their 
quality depends on signal quality, channel availability, and how similar the new data are to the data the model was 
developed on. In practice, automated staging is often used either as a first-pass scorer or as a way to generate 
staging that is later reviewed and corrected manually.

Sleep staging is used both as a primary outcome and as a structural framework for downstream analysis. At the macro 
level, it lets you quantify sleep architecture: total sleep time, sleep efficiency, sleep latency, REM latency, 
wake after sleep onset (WASO), time spent in each stage, stage proportions, fragmentation, and sleep cycles. At the 
microstructural or dynamical level, it supports analyses of stage transitions, bout durations, transition entropy, 
stability of stages, and hypnogram similarity across nights or participants. It is also essential for event-based 
analyses, because many phenomena of interest such as spindles, slow oscillations, REMs, PAC, or power spectra are 
interpreted differently depending on the stage in which they occur. So staging is both a descriptive summary of 
sleep organization and the backbone for most stage-specific sleep neurophysiology analyses.

Seapipe allows both staging to be performed externally (manually) and can perform automatic sleep staging based on 
previously published algorithms:
  1. `Vallat & Walker (2020) <https://elifesciences.org/articles/70092>`_
  2. `Sleep ECG <https://pmc.ncbi.nlm.nih.gov/articles/PMC7355395/>`_ - TO DO 
  3. `Sleep U-Sleep (2021) <https://www.nature.com/articles/s41746-021-00440-5>`_ - TO DO 
  4. `SE-Res-UNet (2025) <https://www.nature.com/articles/s41598-025-00742-8#Abs1>`_ - TO DO 


.. _Functions:
Functions to automatically score staging
----------------

| **Detecting sleep stages will involve these functions:**

1) Detect sleep stages:  

.. code-block:: python

   project.detect_sleep_stages()
|
    This will copy the :ref:`Annotations file` from every ``/sub-XXX/ses-XXX`` in ``<xml_dir>`` to ``<root_dir>/derivatives/staging/`` and write in the detected stages. 
|

2) Summarise hypnogram dynamics (TIDE):

.. code-block:: python

   project.tide()
|
    This will read the scored hypnogram from each :ref:`Annotations file` and export
    transition matrices, bout-duration summaries, and group-level hypnogram
    similarity metrics.
|


.. _detection_staging:
Detect stages
----------------
*Command line argument:*

.. code-block:: python

    project.detect_sleep_stages(xml_dir = None, 
                                out_dir = None, 
                                subs = 'all', 
                                sessions = 'all',
                                filetype = '.edf',
                                method = 'Vallat2021',
                                qual_thresh = 0.5, 
                                eeg_chan = None,
                                ref_chan = None,
                                eog_chan = None,
                                emg_chan = None,
                                rater = None,
                                invert = False,
                                outfile = True)


*Positional arguments:*

    **xml_dir** *(str)*
        * Path to the directory with sub-directories ``/sub-XXX/ses-XXX`` containing the input :ref:`Annotations files<Annotations file>`. 

        * Default is ``None`` which will point to ``<root_dir>/derivatives/staging/`` (Annotations files with sleep stage markings and arousal/artefact events).

    **out_dir** *(str)*
        * Output path for the .xml file containing the new detected event (events will be named like the method used; e.g., ``Ray2015``)

        * Default is ``None`` which will point to ``<root_dir>/derivatives/spindle/``

    **subs** *(str, NoneType or list)*
        * Subject IDs to analyze

        * *Acceptable options:*

            * Default is ``'all'`` which will point to all the ``sub-XXX/`` directories in ``<root_dir>/rawdata/``

            * Entering ``None`` will point seapipe to the *sub* column in the :ref:`tracking file<Tracking File>`

            * Entering a list of sub IDs (e.g., ``['sub-01', 'sub-02']``) will result in detections for those subjects only

    **sessions** *(str, NoneType or list)*
        * Session IDs to analyse per subject

        * *Acceptable options:*

            * Default is ``'all'`` which will point to all the ``ses-XXX/`` directories within the ``sub-XXX/`` directories in ``<root_dir>/rawdata/``

            * Entering ``None`` will point seapipe to the *ses* column in the :ref:`tracking file<Tracking File>`

            * Entering a list of ses visits (e.g., ``['ses-V1', 'ses-V2']``) will result in detections for those session(s) within each subject only

    **filetype** *(str)*
        * Format of files containing EEG signal

        * *Acceptable options:*

            * Default is ``'.edf'`` format

            * The pipeline can also read ``.eeg``, ``.set`` formats

    **method** *(str)*
        * Method(s) of automated detection algorithm to detect staging with. 

        * *Acceptable options:*

            * Currently only ``'Vallat2021'`` is supported. `ref <https://doi.org/10.7554/eLife.70092>`_

    **qual_thresh** *(float)*
        * Quality threshold. Any stages with a confidence of prediction lower than this threshold will be set to ``'Undefined'`` for futher manual review.

    **eeg_chan** *(NoneType or str or list)*
        * EEG channel to use for sleep stage detection

        * *Acceptable options:*

            * Default is ``None`` which will point to the *chanset* columns in the :ref:`tracking file<Tracking File>`

            * Entering a list of channel names (e.g., ``['Fz', 'Cz']``) will only detect the selected channels (see NOTE in section :ref:`Channel Names<Channel Names>`)

    **ref_chan** *(NoneType or list)*
        * :ref:`Reference channel(s)<Channel Names>` for the channels of interest (e.g., mastoid A1 or A2 or joint mastoids)

        * *Acceptable options:*

            * Default is ``None`` which will point to the *refset* columns in the :ref:`tracking file<Tracking File>`. **NOTE** If the tracking file or no *refset* columns exist, then channels will not be re-referenced!

            * Entering a list of channel names (e.g., ``['A1', 'A2']``) will re-reference to these channels  

            * Entering an empty list (i.e., ``[]``) will perform no re-referencing

    **eog_chan** *(NoneType or str or list)*
        * EOG channel to use for sleep stage detection

        * *Acceptable options:*

            * Default is ``None`` which will point to the *chanset* columns in the :ref:`tracking file<Tracking File>`

            * Entering in a *str* containing a channel name (e.g., ``'EOGr'``) will use that channel (see NOTE in section :ref:`Channel Names<Channel Names>`)

            * Entering a *list* of channel names (e.g., ``['EOGl', 'EOGr']``) will use all the named channels 

    **emg_chan** *(NoneType or str or list)*
        * EMG channel to use for sleep stage detection

        * *Acceptable options:*

            * Default is ``None`` which will point to the *chanset* columns in the :ref:`tracking file<Tracking File>`

            * Entering in a *str* containing a channel name (e.g., ``'EMG1'``) will use that channel (see NOTE in section :ref:`Channel Names<Channel Names>`)

            * Entering a *list* of channel names (e.g., ``['EMG1', 'EMG2']``) will use all the named channels 

    **rater** *(NoneType or list)*
        * Name of the rater in the :ref:`Annotations file` to save the detections under

        * *Acceptable options:*

            * Default is ``None`` which will discard the name of the rater. 

            .. note::
                This assumes there is only one rater per Annotations file (``.xml``) 
                !! make sure you don't have multiple raters!!
    
            * Entering a list of rater names (e.g., ``['Rater1', 'Rater2']``) will only save detected events on this rater in the Annotations file

    **invert** *(NoneType or logical)*
        * Option to invert polarity

        * *Acceptable options:*

            * Default is ``None`` which will point to the *chanset_invert* columns in the :ref:`tracking file<Tracking File>`. However, if the *tracking* file does not specify *chanset_invert* 
            columns, the detection will default to ``False``

            * Entering ``False`` will keep the polarity of the recording as it is

            * Entering ``True`` will reverse (flip) the polarity of the recording 

    **outfile** *(str or logical)*
        * Logging of detection

        * *Acceptable options:*

            * Default is ``True`` which will create a logfile *detect_spindles_{method}_{datetime}_log.txt* in ``<root_dir>/derivatives/audit/logs/``

            * Entering a string ``<custom_outfile_name.txt>`` will save the logfile under that custom name
            
            * Entering ``False`` won't save a logfile


.. _tide_overview:
TIDE: Transitions, Intervals, and Dynamics of Epochs
----------------------------------------------------

TIDE reads sleep stages from the staging ``.xml`` annotations and exports three
families of hypnogram-derived metrics:

1. **Transition matrices**
    * Row-normalised stage-to-stage transition probabilities.
    * ``resolution='complete'`` exports a 5 x 5 matrix over
      ``['Wake', 'NREM1', 'NREM2', 'NREM3', 'REM']`` by default.
    * ``resolution='reduced'`` collapses ``NREM1``/``NREM2``/``NREM3`` into a
      single ``NREM`` class and exports a 3 x 3 matrix over ``Wake``, ``NREM``,
      and ``REM``.

2. **Stage bout duration distributions**
    * Bout durations are calculated from contiguous runs of the same stage.
    * Duration-valued columns are reported in **minutes** and include:

      * ``stage_mean_bout_dur_min_*``
      * ``stage_median_bout_dur_min_*``
      * ``stage_p75_bout_dur_min_*``

    * Additional per-stage bout metrics include:

      * ``stage_skew_bout_dur_*``
      * ``stage_num_bouts_*``
      * ``stage_prop_short_bouts_*`` (proportion of bouts shorter than 2 min)

    * The same six bout metrics are also exported for
      ``*_all_stages`` (all selected stages pooled together).

3. **Hypnogram similarity**
    * ``hyp_sim_epoch``: epoch-by-epoch agreement after aligning hypnograms from
      sleep onset.
    * ``hyp_sim_kappa``: Cohen's kappa on the aligned hypnograms.
    * ``hyp_sim_transition_corr``: correlation between subject-level transition
      matrices.


.. _tide_outputs:
TIDE Outputs
------------

By default TIDE separates subject/session outputs from group summaries:

* **Subject/session outputs** are written to
  ``<root_dir>/derivatives/hypnogram/sub-XXX/ses-XXX/``

  * ``*_tide_transition_matrix_complete.csv``
  * ``*_tide_transition_counts_complete.csv``
  * ``*_tide_transition_matrix_reduced.csv``
  * ``*_tide_transition_counts_reduced.csv``
  * ``*_tide_stage_duration_distributions.csv``

* **Group-level outputs** are written to
  ``<root_dir>/derivatives/datasets/hypnogram/``

  * ``tide_transition_matrix_complete_summary.csv``
  * ``tide_transition_matrix_reduced_summary.csv``
  * ``tide_stage_duration_distributions_summary.csv``
  * ``hyp_sim_epoch.csv``
  * ``hyp_sim_kappa.csv``
  * ``hyp_sim_transition_corr.csv``
  * ``hypnogram_similarity_manifest.csv``


.. _tide_pipeline:
Run TIDE from the pipeline
--------------------------

*Command line argument:*

.. code-block:: python

    project.tide(xml_dir = None,
                 out_dir = None,
                 subject_out_dir = None,
                 subs = 'all',
                 sessions = 'all',
                 stage = None,
                 rater = None,
                 resolution = 'complete',
                 analyses = 'all',
                 keyword = None,
                 outfile = True)


*Positional arguments:*

    **xml_dir** *(str)*
        * Path to the directory with sub-directories ``/sub-XXX/ses-XXX``
          containing the staging :ref:`Annotations files<Annotations file>`.

        * Default is ``None`` which will point to the staging derivatives
          directory resolved by seapipe (usually
          ``<root_dir>/derivatives/staging/`` or
          ``<root_dir>/derivatives/staging_auto/`` depending on what exists).

    **out_dir** *(str)*
        * Output path for group-level ``.csv`` summaries.

        * Default is ``None`` which will point to
          ``<root_dir>/derivatives/datasets/hypnogram/``

    **subject_out_dir** *(str)*
        * Output path for subject/session-level files.

        * Default is ``None`` which will point to
          ``<root_dir>/derivatives/hypnogram/``

    **subs** *(str, NoneType or list)*
        * Subject IDs to analyse

        * *Acceptable options:*

            * Default is ``'all'`` which will analyse all ``sub-XXX/``
              directories in ``xml_dir``

            * Entering ``None`` will point seapipe to the *sub* column in the
              :ref:`tracking file<Tracking File>`

            * Entering a list of subject IDs (e.g., ``['sub-01', 'sub-02']``)
              will analyse those subjects only

    **sessions** *(str, NoneType or list)*
        * Session IDs to analyse per subject

        * *Acceptable options:*

            * Default is ``'all'`` which will analyse all ``ses-XXX/``
              directories within each ``sub-XXX/`` directory in ``xml_dir``

            * Entering ``None`` will point seapipe to the *ses* column in the
              :ref:`tracking file<Tracking File>`

            * Entering a list of session IDs (e.g., ``['ses-V0', 'ses-V1']``)
              will analyse those sessions only

    **stage** *(NoneType, str or list)*
        * Stages to include in the calculations.

        * *Acceptable options:*

            * Default is ``['Wake', 'NREM1', 'NREM2', 'NREM3', 'REM']``

            * Entering a single stage name (e.g., ``'REM'``) or a list of stage
              names restricts the calculations to those stages

    **rater** *(NoneType or str)*
        * Name of the rater in the :ref:`Annotations file` to read staging from

        * *Acceptable options:*

            * Default is ``None`` which will select the first rater found in the
              file

            * Entering a rater name (e.g., ``'Vallat2021'``) will read staging
              only from that rater

    **resolution** *(str)*
        * Resolution of the transition matrix export.

        * *Acceptable options:*

            * ``'complete'`` exports a 5 x 5 matrix

            * ``'reduced'`` exports a 3 x 3 Wake/NREM/REM matrix

    **analyses** *(str or list)*
        * Which TIDE analyses to run.

        * *Acceptable options:*

            * Default is ``'all'`` which runs:

              * ``'transition_matrix'``
              * ``'stage_duration_distributions'``
              * ``'hypnogram_similarity'``

            * Entering one string (e.g., ``'hypnogram_similarity'``) runs only
              that analysis

            * Entering a list of names runs the selected subset

    **keyword** *(NoneType or str)*
        * Optional substring used to select the correct ``.xml`` file when more
          than one annotations file exists in a subject/session folder.

    **outfile** *(str or logical)*
        * Logging of the analysis

        * *Acceptable options:*

            * Default is ``True`` which will create a logfile in
              ``<root_dir>/derivatives/audit/logs/``

            * Entering a string ``<custom_outfile_name.txt>`` will save the
              logfile under that custom name

            * Entering ``False`` won't save a logfile


.. _tide_examples:
Examples
--------

Run the full TIDE workflow:

.. code-block:: python

    project.tide()


Run only the reduced transition matrices:

.. code-block:: python

    project.tide(analyses = 'transition_matrix',
                 resolution = 'reduced')


Run TIDE on a subset of subjects/sessions and read a specific rater:

.. code-block:: python

    project.tide(subs = ['sub-01', 'sub-02'],
                 sessions = ['ses-V0'],
                 rater = 'Vallat2021')


.. _tide_metric_details:
TIDE metric details
-------------------

**Transition matrices**
    * ``*_transition_matrix_complete.csv`` and
      ``*_transition_matrix_reduced.csv`` contain row-normalised probabilities.
    * Each row sums to 1 within the selected state space.
    * ``*_transition_counts_complete.csv`` and
      ``*_transition_counts_reduced.csv`` contain the underlying transition
      counts used to calculate those probabilities.

**Bout duration metrics**
    * A bout is a contiguous run of epochs with the same stage label.
    * ``stage_mean_bout_dur_min_*``:
      arithmetic mean bout duration in minutes.
    * ``stage_median_bout_dur_min_*``:
      median bout duration in minutes.
    * ``stage_p75_bout_dur_min_*``:
      75th percentile of the bout duration distribution in minutes.
    * ``stage_skew_bout_dur_*``:
      skewness of the bout duration distribution.
    * ``stage_num_bouts_*``:
      number of bouts for that stage.
    * ``stage_prop_short_bouts_*``:
      proportion of bouts shorter than 2 minutes.

**Whole-hypnogram metrics**
    * ``p_stay_same_stage``:
      proportion of epoch-to-epoch transitions that stay in the same stage.
    * ``transition_entropy``:
      entropy of the transition probability structure across the selected stages.
    * ``num_sleep_cycles``:
      approximate number of NREM-to-REM cycles, estimated from bout structure.
    * ``rem_first_half_prop``:
      proportion of epochs scored as REM in the first half of the sleep period.
    * ``rem_second_half_prop``:
      proportion of epochs scored as REM in the second half of the sleep period.
    * ``delta_n3_early_late_ratio``:
      ratio of N3 proportion in the first half versus the second half of the
      sleep period.

**Hypnogram similarity**
    * ``hyp_sim_epoch.csv``:
      pairwise epoch-by-epoch agreement between hypnograms, aligned from sleep
      onset.
    * ``hyp_sim_kappa.csv``:
      pairwise Cohen's kappa between aligned hypnograms.
    * ``hyp_sim_transition_corr.csv``:
      pairwise correlation between subject-level transition matrices.
    * ``hypnogram_similarity_manifest.csv``:
      manifest listing the subject/session IDs included in the similarity
      matrices and the number of epochs used from sleep onset.


.. _tide_internal_api:
Lower-level TIDE class
----------------------

Advanced users can call the lower-level class directly:

.. code-block:: python

    from seapipe.stats.tide import tide

    T = tide(xml_dir,
             out_dir = None,
             stage = None,
             rater = None,
             subs = 'all',
             sessions = 'all',
             keyword = None,
             subject_out_dir = None)

    T.transition_matrix(stage = None, resolution = 'complete')
    T.stage_duration_distributions(stage = None)
    T.hypnogram_similarity(stage = None)


Here:

    **out_dir** controls the group-level outputs, while
    **subject_out_dir** controls the per-subject/per-session files.

    The method-level ``stage`` argument can be used to override the stage list
    stored on the class instance.










