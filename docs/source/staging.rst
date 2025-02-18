Staging
=====

.. _Overview:

Overview
------------

Seapipe provides some automatic sleep staging algorithms:

1. `Vallat & Walker (2020) <https://elifesciences.org/articles/70092>`_
2.  'Sleep ECG' - TO DO



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

            * Default is ``'all'`` which will point to all the ``sub-XXX/`` directories in ``<root_dir>/DATA/``

            * Entering ``None`` will point seapipe to the *sub* column in the :ref:`tracking file<Tracking File>`

            * Entering a list of sub IDs (e.g., ``['sub-01', 'sub-02']``) will result in detections for those subjects only

    **sessions** *(str, NoneType or list)*
        * Session IDs to analyse per subject

        * *Acceptable options:*

            * Default is ``'all'`` which will point to all the ``ses-XXX/`` directories within the ``sub-XXX/`` directories in ``<root_dir>/DATA/``

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











