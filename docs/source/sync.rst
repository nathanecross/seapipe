Event synchrony
=====

.. _overview_sync:

Overview
------------

Event synchrony (CORAL: *Co-Occurrence and Recall Analysis of Labels*) splits a
**target** event into two classes based on whether a **probe** event overlaps in time.
This is useful when you want to separate events that co-occur with another event
(e.g., spindles that occur during slow oscillations) from those that do not.

| CORAL writes new event labels into a copied :ref:`Annotations file` and can also
  export a summary dataset with recall/precision/F1 statistics.


.. _detection_sync:
Run Event Synchrony
----------------
|
    This will copy the :ref:`Annotations file` from every ``/sub-XXX/ses-XXX`` in
    ``<xml_dir>`` to ``<root_dir>/derivatives/sync/`` and add two new event types:
    ``evttype_tp_target`` (co-occurring) and ``evttype_fn`` (non co-occurring, optional).
|

*Command line argument:*

.. code-block:: python

    project.event_synchrony(xml_dir = None, out_dir = None,
                            subs = 'all', sessions = 'all',
                            chan = None, stage = None,
                            grp_name = 'eeg', rater = None,
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
        * Path to the directory with sub-directories ``/sub-XXX/ses-XXX`` containing
          the input :ref:`Annotations files<Annotations file>`.

        * Default is ``None`` which will point to ``<root_dir>/derivatives/<evttype_target>/``
          (or a known derivatives folder for common event names).

    **out_dir** *(str)*
        * Output path for the new annotations with synchronized events.

        * Default is ``None`` which will point to ``<root_dir>/derivatives/sync/``

    **subs** *(str, NoneType or list)*
        * Subject IDs to analyze

        * *Acceptable options:*

            * Default is ``'all'`` which will point to all the ``sub-XXX/`` directories
              in ``<root_dir>/rawdata/``

            * Entering ``None`` will point seapipe to the *sub* column in the
              :ref:`tracking file<Tracking File>`

            * Entering a list of sub IDs (e.g., ``['sub-01', 'sub-02']``) will result
              in detections for those subjects only

    **sessions** *(str, NoneType or list)*
        * Session IDs to analyse per subject

        * *Acceptable options:*

            * Default is ``'all'`` which will point to all the ``ses-XXX/`` directories
              within the ``sub-XXX/`` directories in ``<root_dir>/rawdata/``

            * Entering ``None`` will point seapipe to the *ses* column in the
              :ref:`tracking file<Tracking File>`

            * Entering a list of ses visits (e.g., ``['ses-V1', 'ses-V2']``) will result
              in detections for those session(s) within each subject only

    **chan** *(NoneType or list)*
        * Channel(s) of interest

        * *Acceptable options:*

            * Default is ``None`` which will point to the *chanset* columns in the
              :ref:`tracking file<Tracking File>`

            * Entering a list of channel names (e.g., ``['Fz', 'Cz']``) will only
              detect the selected channels

    **stage** *(list)*
        * Stages of interest

        * *Acceptable options:*

            * Default is ``['NREM2', 'NREM3']`` if staging is available

            * Entering a list of stages (e.g., ``['NREM3']``) will only analyze those
              stages

    **evttype_target** *(str)*
        * Name of the target event to split based on co-occurrence

    **evttype_probe** *(str)*
        * Name of the probe event used to test co-occurrence

    **evttype_tp_target** *(str)*
        * Name of the event written for target events that co-occur with probe events

    **evttype_fn** *(str, NoneType)*
        * Name of the event written for target events that do *not* co-occur with probe
          events (optional)

    **iu_thresh** *(float)*
        * Intersection-over-union threshold to count a co-occurrence (default ``0.5``)


.. _dataset_sync:
Export Event Synchrony Dataset
----------------
|
    This will export a summary ``.csv`` with recall, precision, F1 score and counts
    of true/false positives/negatives for each subject and session.
|

*Command line argument:*

.. code-block:: python

    project.event_synchrony_dataset(xml_dir = None, out_dir = None,
                                    subs = 'all', sessions = 'all',
                                    chan = None, stage = None,
                                    grp_name = 'eeg', rater = None,
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

    **out_dir** *(str)*
        * Output path for the dataset ``.csv``

        * Default is ``None`` which will point to ``<root_dir>/derivatives/datasets/sync/``

    **outfile_suffix** *(str)*
        * Filename for the dataset ``.csv``

        * Default is ``<evttype_target>_x_<evttype_probe>_sync_stats.csv``
