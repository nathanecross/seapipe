Event synchrony
===============

.. _overview_sync:

Overview
--------

Event synchrony in seapipe is handled by CORAL: *Co-Occurrence and Recall
Analysis of Labels*. CORAL compares two event types in the same annotations
file:

* a **target** event that you want to split
* a **probe** event that defines co-occurrence

The target events are matched to probe events using an
intersection-over-union threshold (``iu_thresh``). Matched target events can
be written back to a copied annotations file under a new event name, and the
same comparison can also be summarised into a dataset of recall, precision,
F1-score, true positives, false positives, and false negatives.

Typical uses include:

* separating spindles that co-occur with slow oscillations from those that do not
* marking slow oscillations that contain a spindle
* creating cohort-level synchrony summary tables for downstream statistics


.. _functions_sync:

Functions
---------

| **Event synchrony can be run in two stages:**

1) Create synchronised events in copied annotations files:

.. code-block:: python

   project.event_synchrony()

|
    This copies the input :ref:`Annotations file` from each
    ``/sub-XXX/ses-XXX`` folder into ``<root_dir>/derivatives/sync/`` and writes
    new event labels for target events that do or do not co-occur with the
    probe event.
|

2) Export a summary dataset:

.. code-block:: python

   project.event_synchrony_dataset()

|
    This creates a ``.csv`` dataset under
    ``<root_dir>/derivatives/datasets/sync/`` summarising synchrony metrics for
    each subject and session.
|


.. _detection_sync:

Run Event Synchrony
-------------------

*Command line argument:*

.. code-block:: python

    project.event_synchrony(xml_dir = None,
                            out_dir = None,
                            subs = 'all',
                            sessions = 'all',
                            stage = None,
                            chan = None,
                            ref_chan = None,
                            grp_name = 'eeg',
                            rater = None,
                            evttype_target = None,
                            evttype_probe = None,
                            evttype_tp_target = None,
                            evttype_fn = None,
                            iu_thresh = 0.5,
                            concat_stage = True,
                            concat_cycle = True,
                            reject_artf = None,
                            filetype = ('.edf', '.rec', '.eeg'),
                            outfile = True)


*Positional arguments:*

    **xml_dir** *(str)*
        * Path to the directory with sub-directories ``/sub-XXX/ses-XXX``
          containing the input :ref:`Annotations files<Annotations file>`.

        * Default is ``None`` which resolves from ``evttype_target`` using
          seapipe's derivatives conventions via ``select_input_dirs()``.

    **out_dir** *(str)*
        * Output path for copied annotations containing the synchrony labels.

        * Default is ``None`` which points to ``<root_dir>/derivatives/sync/``

    **subs** *(str, NoneType or list)*
        * Subject IDs to analyze.

        * *Acceptable options:*

            * Default is ``'all'`` which uses all ``sub-XXX/`` directories in
              ``<root_dir>/rawdata/``

            * Entering ``None`` uses the *sub* column in the
              :ref:`tracking file<Tracking File>`

            * Entering a list of subject IDs (e.g.,
              ``['sub-01', 'sub-02']``) restricts analysis to those subjects

    **sessions** *(str, NoneType or list)*
        * Session IDs to analyse per subject.

        * *Acceptable options:*

            * Default is ``'all'`` which uses all ``ses-XXX/`` directories for
              each subject

            * Entering ``None`` uses the *ses* column in the
              :ref:`tracking file<Tracking File>`

            * Entering a list of sessions (e.g., ``['ses-V1', 'ses-V2']``)
              restricts analysis to those visits

    **stage** *(NoneType or list)*
        * Sleep stage(s) to include during event fetching.

        * *Acceptable options:*

            * Default is ``None``, which means the value already configured in
              the pipeline call will be passed through unchanged

            * Entering a list of stages (e.g., ``['NREM2', 'NREM3']``) limits
              synchrony analysis to those stages

    **chan** *(NoneType or list)*
        * Channel(s) of interest.

        * *Acceptable options:*

            * Default is ``None`` which loads channels from the *chanset*
              columns in the :ref:`tracking file<Tracking File>`

            * Entering a list of channel names (e.g., ``['Fz', 'Cz']``)
              restricts analysis to those channels

    **ref_chan** *(NoneType or list)*
        * Reference channel(s) corresponding to the selected channels.

        * *Acceptable options:*

            * Default is ``None`` which loads references from the *refset*
              columns in the :ref:`tracking file<Tracking File>`

            * Entering a list of channel names (e.g., ``['A1', 'A2']``)
              supplies the reference set directly

    **grp_name** *(str)*
        * Name of the annotation group used for event fetching.

        * Default is ``'eeg'``

    **rater** *(NoneType or str)*
        * Name of the rater in the :ref:`Annotations file`.

        * Default is ``None``

    **evttype_target** *(str)*
        * Name of the target event to be split by co-occurrence.

    **evttype_probe** *(str)*
        * Name of the probe event used to test co-occurrence.

    **evttype_tp_target** *(str)*
        * Output event name for target events that co-occur with the probe.

    **evttype_fn** *(NoneType or str)*
        * Optional output event name for target events that do not co-occur
          with the probe.

    **iu_thresh** *(float)*
        * Intersection-over-union threshold used by Wonambi's
          ``match_events()``.

        * Default is ``0.5``

    **concat_stage** *(logical)*
        * Stage concatenation setting used during fetch.

        * Default is ``True``

    **concat_cycle** *(logical)*
        * Cycle concatenation setting used during fetch.

        * Default is ``True``

    **reject_artf** *(NoneType or list)*
        * Event labels to reject during fetching.

        * Default is ``None`` which falls back to CORAL's internal defaults:
          ``['Artefact', 'Arou', 'Arousal']``

    **filetype** *(tuple or list)*
        * Recording file extensions to search for.

        * Default is ``('.edf', '.rec', '.eeg')``

    **outfile** *(bool or str)*
        * Controls logging output.

        * Default is ``True`` which creates a timestamped log file


Example
-------

.. code-block:: python

    project.event_synchrony(
        xml_dir = None,
        subs = 'all',
        sessions = 'all',
        chan = ['Cz'],
        ref_chan = ['A1', 'A2'],
        stage = ['NREM2', 'NREM3'],
        evttype_target = 'spindle',
        evttype_probe = 'slowwave',
        evttype_tp_target = 'spindle_so',
        evttype_fn = 'spindle_noso',
        iu_thresh = 0.5,
    )


.. _dataset_sync:

Export Event Synchrony Dataset
------------------------------

*Command line argument:*

.. code-block:: python

    project.event_synchrony_dataset(xml_dir = None,
                                    out_dir = None,
                                    subs = 'all',
                                    sessions = 'all',
                                    chan = None,
                                    stage = None,
                                    grp_name = 'eeg',
                                    rater = None,
                                    evttype_target = None,
                                    evttype_probe = None,
                                    evttype_tp_target = None,
                                    evttype_fn = None,
                                    iu_thresh = 0.5,
                                    concat_stage = True,
                                    concat_cycle = True,
                                    outfile_suffix = None,
                                    reject_artf = None,
                                    filetype = ('.edf', '.rec', '.eeg'),
                                    outfile = True)


*Positional arguments:*

    **xml_dir** *(str)*
        * Path to the input annotations directory.

        * Default is ``None`` which resolves from ``evttype_target`` using
          ``select_input_dirs()``

    **out_dir** *(str)*
        * Output path for the synchrony summary dataset.

        * Default is ``None`` which points to
          ``<root_dir>/derivatives/datasets/sync/``

    **subs** *(str, NoneType or list)*
        * Subject IDs to include in the dataset.

    **sessions** *(str, NoneType or list)*
        * Sessions to include in the dataset.

    **chan** *(NoneType or list)*
        * Channel(s) to summarise.

    **stage** *(NoneType or list)*
        * Sleep stage(s) to summarise.

    **grp_name** *(str)*
        * Annotation group name used for fetching.

        * Default is ``'eeg'``

    **rater** *(NoneType or str)*
        * Rater name in the annotations file.

    **evttype_target** *(str)*
        * Target event name.

    **evttype_probe** *(str)*
        * Probe event name.

    **evttype_tp_target** *(str)*
        * Label used for matched target events.

    **evttype_fn** *(NoneType or str)*
        * Optional label used for unmatched target events.

    **iu_thresh** *(float)*
        * Intersection-over-union threshold.

        * Default is ``0.5``

    **concat_stage** *(logical)*
        * Stage concatenation setting used during fetch.

    **concat_cycle** *(logical)*
        * Cycle concatenation setting used during fetch.

    **outfile_suffix** *(str)*
        * Filename for the summary dataset.

        * Default is
          ``<evttype_target>_x_<evttype_probe>_sync_stats.csv``

    **reject_artf** *(NoneType or list)*
        * Event labels to reject during fetching.

        * Default is ``None``

    **filetype** *(tuple or list)*
        * Recording file extensions to search for.

        * Default is ``('.edf', '.rec', '.eeg')``

    **outfile** *(bool or str)*
        * Controls logging output.

        * Default is ``True``


Outputs
-------

``project.event_synchrony()`` writes copied XML annotations with new event
labels under ``<root_dir>/derivatives/sync/sub-XXX/ses-XXX/``.

``project.event_synchrony_dataset()`` writes one or more ``.csv`` files under
``<root_dir>/derivatives/datasets/sync/`` containing:

* ``Recall``
* ``Precision``
* ``F1 score``
* count of matched target events (``evttype_tp_target``)
* count of unmatched probe events
* count of unmatched target events (if ``evttype_fn`` is requested)


Notes
-----

* Both wrappers build the fetch concatenation tuple as
  ``(int(concat_cycle), int(concat_stage), 0, 0)`` so event types are not
  concatenated during synchrony analysis.
