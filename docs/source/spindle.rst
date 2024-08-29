Spindles
=====

.. _overview:

Overview
------------
TODO
Slow oscillations (SOs) are coherent waves corresponding to the alternation between biphasic membrane potential levels (UP states = depolarization 
and DOWN states = hyperpolarization). Oscillating below ~1 Hz, SOs are generated during sleep stages NREM2 and NREM3.

| Slow oscillations can be detected as events and their characteristics (see definitions in section :ref:`Output<Outputs of Slow Oscillations>`) can be extracted across NREM (NREM2+NREM3), per stage and/or per cycle.

| Seapipe provides 4 published methods to automatically detect SOs:

    * `Staresina et al. (2015) <https://doi.org/10.1038/nn.4119>`_: SO detection (<1.25Hz) with an adapted amplitude criteria per individual and channel.
    
       Method in brief:

        1. Filter the signal (FIR bandpass filter, 0.5–1.25 Hz, order = 3); 

        2. Identify events with a positive-to-negative zero crossing and a subsequent negative-to-positive zero crossing separated by 0.8-2 sec;

        3. Keep the top 25% of events with the largest trough-to-peak amplitudes. 




.. _Functions:
Functions to detect spindles
----------------
| **Detecting spindles and extracting their parameters will involve three functions:**

1) Detect spindles events:  

.. code-block:: python

   project.detect_spindles()
|
    This will copy the :ref:`Annotations file` from every ``/sub-XXX/ses-XXX`` in ``<xml_dir>`` to ``<root_dir>/OUT/spindle/`` and write in the detected events. 
|
2) Export event characteristics: 

.. code-block:: python

   project_name.export_eventparams()
|   
    This will extract a ``.csv`` file for every channel and/or stage and/or cycle into every ``/sub-XXX/ses-XXX`` directory in ``<root_dir>/OUT/slowwave/`` 
|
3) Create datasets combining all the subjects: 

.. code-block:: python

   project_name.event_dataset()
|
    This will combine all of the ``.csv`` files from the previous step into a single dataset (one row per subject) ``<root_dir>/OUT/datasets/``
|

.. _detection_spindle:
Detect spindles
----------------
*Command line argument:*

.. code-block:: python

    project.detect_spindles(xml_dir = None, 
                            out_dir = None, 
                            subs = 'all', 
                            sessions = 'all', 
                            filetype = '.edf', 
                            method = ['Moelle2011'], 
                            chan = None, 
                            ref_chan = None, 
                            rater = None, 
                            stage = ['NREM2','NREM3'], 
                            grp_name = 'eeg', 
                            cycle_idx = None, 
                            concat_cycle = True, ##should be False
                            frequency = None, 
                            adap_bands = 'Fixed', 
                            adap_bw = 4, 
                            duration =( 0.5, 3),
                            reject_artf = ['Artefact', 'Arou', 'Arousal'], 
                            outfile = True)


*Positional arguments:*

    **xml_dir** *(str)*
        * Path to the directory with sub-directories ``/sub-XXX/ses-XXX`` containing the input :ref:`Annotations files<Annotations file>`. 

        * Default is ``None`` which will point to ``<root_dir>/OUT/staging/`` (Annotations files with sleep stage markings and arousal/artefact events).

    **out_dir** *(str)*
        * Output path for the .xml file containing the new detected event (events will be named like the method used; e.g., ``Ray2015``)

        * Default is ``None`` which will point to ``<root_dir>/OUT/spindle/``

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

    **method** *(list)*
        * Method(s) of interest for spindles detection

        * *Acceptable options:*

            * Default is ``['Moelle2011']`` method  
            
            * All methods can be run simultaneously (e.g., ``['Ray2015', 'Lacourse2018']``)

    **chan** *(NoneType or list)*
        * Channel(s) of interest

        * *Acceptable options:*

            * Default is ``None`` which will point to the *chanset* columns in the :ref:`tracking file<Tracking File>`

            * Entering a list of channel names (e.g., ``['Fz', 'Cz']``) will only detect the selected channels (see NOTE in section :ref:`Channel Names<Channel Names>`)

    **ref_chan** *(NoneType or list)*
        * :ref:`Reference channel(s)<Channel Names>` for the channels of interest (e.g., mastoid A1 or A2 or joint mastoids)

        * *Acceptable options:*

            * Default is ``None`` which will point to the *refset* columns in the :ref:`tracking file<Tracking File>`. **NOTE** If the tracking file or no *refset* columns exist, then channels will not be re-referenced!

            * Entering a list of channel names (e.g., ``['A1', 'A2']``) will re-reference to these channels  

            * Entering an empty list (i.e., ``[]``) will perform no re-referencing

    **rater** *(NoneType or list)*
        * Name of the rater in the :ref:`Annotations file` to save the detections under

        * *Acceptable options:*

            * Default is ``None`` which will discard the name of the rater. 

            .. note::
                This assumes there is only one rater per Annotations file (``.xml``) 
                !! make sure you don't have multiple raters!!
    
            * Entering a list of rater names (e.g., ``[<Rater1>, <Rater2>]``) will only save detected events on this rater in the Annotations file

    **stage** *(list)*
        * Stages of interest

        * *Acceptable options:*

            * Default is ``['NREM2', 'NREM3']`` 

            * Entering a list of stages (e.g., ``['NREM3']``), it will only detect the events for this specific stage. **It is recommended that you leave the default option**


    **grp_name** *(str)*
        * Name of the tab in the :ref:`Annotations file` to save the detections to. This is for visualization in Wonambi only, however it will impact the `exporting <Export slow oscillations characteristics>` of events also

        * *Acceptable options:*

            * Default is ``eeg`` which is the recommended naming convention
           
            * Entering a list of group names (e.g., ``['eeg_hemiR']``) will save the events to a tab of this name in the Annotations file. The events can only be visualised in :ref:`Wonambi` with a montage that includes a tab with this name

    **cycle_idx** *(NoneType or tuple)*
        * Sleep cycle numbers

        * *Acceptable options:*

            * Default is ``None`` which will infer no cycles 

            * Entering a `tuple <https://docs.python.org/3/tutorial/datastructures.html#tuples-and-sequences>`_ of integers corresponding to sleep cycle numbers (e.g., ``(1,2,3,4,5)``), it will only detect the events for these specific 
            cycles' numbers. If a ``sub`` has less than the number of cycles entered, then the maximum number of cycles possible will be used for that subject.

    **concat_cycle** *(logical)*
        * Concatenation options for sleep cycles

        * *Acceptable options:*

            * Default is ``False`` which means that detection will be performed per stage

            * Entering ``True`` which means that all cycles will be concatenated (i.e., merged) before detection **It is recommended that you leave the default option**

    **frequency** *(tuple)*
        * Frequency range of interest 

        * *Acceptable options:*

            * Default is ``None`` which will depend to the options selected for **adap_bands**. If ``adap_band = 'Fixed'``, frequency will be (11,16) while ``adap_band = 'Auto'``
            will be (9,16) for the peak frequency detection

            * Enter a `tuple <https://docs.python.org/3/tutorial/datastructures.html#tuples-and-sequences>`_ containing the frequency range of interest that 
            will be used if selecting ``adap_bands = 'Fixed'`` or ``adap_bands = 'Auto'`
 

    **adap_bands** *(str)*
        * Options to set an adapted sigma band of spindle detection tailored to each individual based on their peak in sigma per channel, stage and/or session

        * *Acceptable options:*

            * Default is ``'Fixed'`` which will point to the frequency range set up in **frequency**

            * Entering ``'Auto'`` will perform :ref:`FOOOF analyses<FOOOF analyses>` which will detect the peak in sigma characterized in terms of their specific
            center frequency, power and bandwidth within the frequency range set up in **frequency** and controlling for the aperiodic component. By default, if left 
            ``frequency = None``, the range set-up for fooof peak detection is 9-16Hz. It will add *_adap* at the end of the event name (e.g., Moelle2011_adap).

            * Entering ``Manual`` will point to the *chanset_peaks* columns in the :ref:`tracking file<Tracking File>`. It will add *_adap* at the end of the event name (e.g., Moelle2011_adap).

    **adap_bw** *(str)*
        * Size of the frequency range around sigma peak frequency when entering ``Auto``or ``Manual`` to **adap_bands**

        * *Acceptable options:*

            * Default is ``4``meaning 2Hz on both side of the sigma peak frequency

            * Enter a number 

    **duration** *(tuple)*
        * Minimum and maximum duration of events that will be detected. Any events with durations that are outside these limits will be discarded

        * *Acceptable options:*

            * Default is ``(0.5, 3)`` (in seconds)

            * Entering a `tuple <https://docs.python.org/3/tutorial/datastructures.html#tuples-and-sequences>`_ of float with length 2 (e.g., ``(0.5, 2)``)  will limit the detection to events with a duration within this range

    **reject_artf** *(list)*
        * Options to discard detection within specific events such as Artefact events

        * *Acceptable options:*

            * Default is ``['Artefact', 'Arou', 'Arousal']``which will discard detection during events with these specific names

            * Entering a list of events will discard detection within those events

    **outfile** *(str or logical)*
        * Logging of detection

        * *Acceptable options:*

            * Default is ``True`` which will create a logfile *detect_spindles_{method}_{datetime}_log.txt* in ``<root_dir>/OUT/audit/logs/``

            * Entering a string ``<custom_outfile_name.txt>`` will save the logfile under that custom name
            
            * Entering ``False`` won't save a logfile


.. _export_spindle:
Export spindle characteristics
----------------
*Command line argument:*

.. code-block:: python


    project.export_eventparams(evt_name,
                               frequency = None,
                               xml_dir = None, 
                               out_dir = None, 
                               subs = 'all', 
                               sessions = 'all', 
                               chan = None, 
                               ref_chan = None, 
                               stage = ['NREM2','NREM3'], 
                               grp_name = 'eeg',
                               rater = None, 
                               cycle_idx = None, 
                               concat_cycle = True, 
                               concat_stage = False, 
                               keyword = None, 
                               segs = None,
                               adap_bands = 'Fixed',
                               adap_bw = 4,
                               params = 'all',  
                               epoch_dur = 30, 
                               outfile = True)

*Required arguments:*

    **evt_name** *(list or str)*
        * Name of the event of interest to export from the :ref:`Annotations file` 

        * Enter a string (e.g ``Ray2015``) which refers to the event as it is named in the Annotations file. **NOTE** This will be the method name used in the :ref:`detection<Detect spindle>`

        * Entering a list of event names (e.g ``['Ray2015', 'Lacourse2018']``) will export the parameters for each event *separately*

*Positional arguments:*

    **frequency** (tuple)
        * Frequency range of interest in which to export event parameters (e.g. *frequency*, *power*)

        * *Acceptable options:*

            * Default is ``None`` which will depend to the options selected for **adap_bands**. If ``adap_band = 'Fixed'``, frequency will be (11,16) while ``adap_band = 'Auto'``
            will be (9,16) for the peak frequency detection

            * Enter a `tuple <https://docs.python.org/3/tutorial/datastructures.html#tuples-and-sequences>`_ containing the frequency range of interest that 
            will be used if selecting ``adap_bands = 'Fixed'`` or ``adap_bands = 'Auto'`

    **xml_dir** *(str)*
        * Path to the directory with sub-directories ``/sub-XXX/ses-XXX`` containing the :ref:`Annotations files<Annotations file>` where the :ref:`detections<Detect spindle>` were saved. 

            * Default is ``None`` which will point to ``<root_dir>/OUT/spindle/``

    **out_dir** *(str)*
        * Output path for the where to save the ``.csv`` file containing the parameters of the spindle events per subject, session, and/or stage, and/or channel.

            * Default is ``None`` which will point to ``root_dir/OUT/spindle/``

    **subs** *(str, NoneType or list)*
        * Subject IDs to analyze

        * *Acceptable options:*

            * Default is ``'all'`` which will point to all the ``sub-XXX/`` directories in ``<root_dir>/DATA/``

            * Entering ``None`` will point seapipe to the *sub* column in the :ref:`tracking file<Tracking File>`

            * Entering a list of sub IDs (e.g., ``['sub-01', 'sub-02']``) will the parameters of the spindle events for those subjects only

    **sessions** *(str, NoneType or list)*
        * Session IDs to analyze

        * *Acceptable options:*

            * Default is ``'all'`` which will point to all the ``ses-XXX/`` directories within the ``sub-XXX/`` directories in ``<root_dir>/DATA/``

            * Entering ``None`` will point seapipe to the *ses* column in the :ref:`tracking file<Tracking File>`

            * Entering a list of session IDs (e.g., ``['ses-V1', 'ses-V2']``) will result in detections for those session(s) within each subject only

    **chan** *(NoneType or list)*
        * Channel(s) of interest

        * *Acceptable options:*

            * Default is ``None`` which will point to the *chanset* columns in the :ref:`tracking file<Tracking File>`

            * Entering a list of channel names (e.g., ``['Fz', 'Cz']``) will only export parameters for the events in the selected channels (see NOTE in section :ref:`Channel Names<Channel Names>`)

    **ref_chan** *(NoneType or list)*
        * :ref:`Reference channel(s)<Channel Names>` for the channels of interest (e.g., mastoid A1 or A2 or joint mastoids)

        * *Acceptable options:*

            * Default is ``None`` which will point to the *refset* columns in the :ref:`tracking file<Tracking File>`. **NOTE** If the tracking file or no *refset* columns exist, then channels will not be re-referenced!

            * Entering a list of channel names (e.g., ``['A1', 'A2']``) will re-reference to these channels  

            * Entering an empty list (i.e., ``[]``) will perform no re-referencing

        .. note::
            If the reference channels are not the same as were entered in the :ref:`detection spindle<Detect spindles>`, the event parameters will still be stored,
            however the parameters (e.g. frequency, amplitude, power) might be affected. Be careful to remain consistent across these steps!           

    **stage** *(list)*
        * Stages of interest

        * *Acceptable options:*

            * Default is ``['NREM2', 'NREM3']`` 

            * Entering a list of stages (e.g., ``['NREM3']``), it will only export parameters for the events in this specific stage

    **grp_name** *(str)*
        * Name of the tab in the :ref:`Annotations file` where the detected events are saved 

        * *Acceptable options:*

            * Default is ``eeg`` which is the recommended naming convention
           
            * If entering a list of group names (e.g., ``['eeg_hemiR']``), ensure that this matches ``grp_name`` used in the :ref:`detection spindle<Detect spindle>`

    **rater** *(NoneType or list)*
        * Name of the rater in the :ref:`Annotations file` where the detected events are saved

        * *Acceptable options:*

            * Default is ``None`` which will discard the name of the rater

            * Entering a list of raters names (e.g., ``[<Rater1>, <Rater2>]``) will only export the the parameters for events saved under this rater. **NOTE** An create an empty extraction ``.csv`` file will be created if the rater is absent

    **cycle_idx** *(NoneType or tuple)*
        * Sleep cycle numbers

        * *Acceptable options:*

            * Default is ``None`` which will infer no cycles 

            * Entering a `tuple <https://docs.python.org/3/tutorial/datastructures.html#tuples-and-sequences>`_ of integers corresponding to sleep cycle numbers (e.g., ``(1,2,3,4,5)``) will only detect the events for these specific 
            cycles. If a ``sub`` has less than the number of cycles entered, then the maximum number of cycles possible will be used for that subject

    **concat_cycle** *(logical)*
        * Concatenation options for sleep cycles

        * *Acceptable options:*

            * Default is ``True`` which means that all cycles will be concatenated (i.e., merged) before exporting the parameters of the spindle events

            * Entering ``False`` will export the spindle parameters per sleep cycle (saving each cycle as a separate ``.csv`` output file)

    **concat_stage** *(logical)*
        * Concatenation options for stages

        * *Acceptable options:*

            * Default is ``False`` which means the parameters of spindle events will be exported per stage (e.g. NREM2, NREM3) separately (saving each stage as a separate ``.csv`` output file)

            * Entering ``True`` will concatenate (i.e., merge) all stages before exporting the parameters of the spindle events

    **keyword** *(str)*
        * Allow seapipe to search for a Annotations filename containing a specific wildcard (keyword)

        * *Acceptable options:*

            * Default is ``None`` which will infer no keyword to search for.

            * Entering a string (e.g. ``Moelle_adapted_custom``) will only export event parameters from this specific Annotations file

    **seg** *(NoneType or list of tuples)*
        * Option to extract parameters of SOs that only occur in between certain markers. These markers need to be events saved in the :ref:`Annotations file`

        * *Acceptable options:*

            * Default is ``None`` which will infer no segmentation prior to exporting event parameters

            * Entering a list of `tuples <https://docs.python.org/3/tutorial/datastructures.html#tuples-and-sequences>`_, with both start and end tags named (e.g. ``[('N2_ON', 'N2_OFF'), ('N3_ON', 'N3_OFF')]``) will export event parameters that only occur between these event markers

    **adap_bands** *(str)*
        * Options to set an adapted sigma band for the extraction of spindle parameters. Enter the options used in the :ref:`detection spindle<Detect spindle>`

        * *Acceptable options:*

            * Default is ``'Fixed'`` which will point to the frequency range set up in **frequency**. By default, if left ``frequency = None``, the range set-up is (11-16).

            * Entering ``'Auto'`` will extract parameters of the events set up in **evt_name** with *_adap* at the end of the name. By default, if left ``frequency = None``, the range set-up is (9-16).

            * Entering ``Manual`` will point to the *chanset_peaks* columns in the :ref:`tracking file<Tracking File>`. It will extract parameters of the events set up in **evt_name** with *_adap* at the end of the name.

    **adap_bw** *(str)*
        * Size of the frequency range around sigma peak frequency when entering ``Auto``or ``Manual`` to **adap_bands**

        * *Acceptable options:*

            * Default is ``4``meaning 2Hz on both side of the sigma peak frequency

            * Enter a even number 

    **params** *(str or dict)*
        * The names of specific parameters to export

        * *Acceptable options:*

            * Default is ``all`` which will export all characteristics (see :ref:`Output`) -  *Recommended* 

            * To specify only specific parameters to export, enter a `dictionary <https://docs.python.org/3/tutorial/datastructures.html#dictionaries>`_ with ``True`` or ``False`` for each parameter (e.g., ``params = ['dur':True, 'minamp':False, 'maxamp':False, 'ptp':True, 'rms':False, 'power':True, 'peakpf':False, 'energy':False, 'peakef':False]``)

    **epoch_dur** *(int)*
        * Options to change the denominator (duration) for the *spindle density* index 

        * *Acceptable options:*

            * Default is ``30`` (this infers 30-second epochs)

            * Entering a number (e.g., ``60``) will imply that the SO density value equals the number of events per this time period (e.g. per *60 seconds*)

    **outfile** *(str or logical)*
        * Logging of event parameter export

        * *Acceptable options:*

            * Default is ``True`` which will create a logfile *export_params_{method}_{datetime}_log.txt* in ``<root_dir>/OUT/audit/logs/``

            * Entering a string ``<custom_outfile_name.txt>`` will save the logfile under that custom name
            
            * Entering ``False`` won't save a logfile

     .. note::
        By default

        * *export_eventparams* cannot extract spindle characteristics without required arguments for ``evt_name``.

        * *export_eventparams* will extract characteristics per stage (NREM2 vs NREM3). If you want the extraction for NREM2+NREM3 combined as well, re-run *export_eventparams* with ``concat_stage = True``.

        * *export_eventparams* will extract characteristics for the whole-night. If you want the extraction per cycle, re-run *export_eventparams* with ``concat_cycle = False``.



.. _create_datasets:
Create datasets
----------------
*Command line argument:* 

.. code-block:: python

   project.event_dataset(chan, 
                         evt_name, 
                         xml_dir = None, 
                         out_dir = None, 
                         subs = 'all', 
                         sessions = 'all', 
                         stage = None,
                         concat_stage = False, 
                         concat_cycle = True,
                         cycle_idx = None, 
                         grp_name = 'eeg',
                         adap_bands = 'Fixed',
                         params = 'all', 
                         outfile=True))


*Required arguments:*
    **chan** *(str or list)*
        * Channel(s) of interest

        * *Acceptable options:*

            * Entering a string (e.g ``Fz``) will create separate datasets for that channel only.

            * Entering a list of channels' names (e.g., ``['Fz', 'Cz', 'Pz']``) will create separate datasets for each channel. The the names will be taken from the *chanset_rename* columns in the :ref:`tracking file<Tracking File>`

    **evt_name** *(str or list)*
        * Name of the events of interest 

        * *Acceptable options:*
        
            * Enter a string (e.g ``Ray2015``) which refers to the event as it was used in the :ref:`export event parameters step<Export spindle characteristics>`

            * Entering a list of event names (e.g ``['Ray2015', 'Lacourse2018']``) will create a dataset for each event *separately*

*Positional arguments:*
    **xml_dir** *(str)*
        * Path to the directory with sub-directories ``/sub-XXX/ses-XXX`` containing the ``.csv`` files from the :ref:`export event parameters step<Export spindle characteristics>`

        * *Acceptable options:*

            * Default is ``None`` which will point to ``<root_dir>/OUT/spindle/``

    **out_dir** *(str)*
        * Output path for the created datasets

        * *Acceptable options:*

            * Default is ``None`` which will point to ``<root_dir>/OUT/datasets/``

    **subs** *(str, NoneType or list)*
        * Subject IDs to export into the dataset

        * *Acceptable options:*

            * Default is ``'all'`` which will point to all the ``sub-XXX/`` directories in ``<xml_dir>/``

            * Entering ``None`` will point seapipe to the *sub* column in the :ref:`tracking file<Tracking File>`

            * Entering a list of sub IDs (e.g., ``['sub-01', 'sub-02']``) will export those subjects only into the dataset 

    **sessions** *(str, NoneType or list)*
        * Session IDs to export into the dataset

        * *Acceptable options:*

            * Default is ``'all'`` which will point to all the ``ses-XXX/`` directories within the ``sub-XXX/`` directories in ``<xml_dir>/``

            * Entering ``None`` will point seapipe to the *ses* column in the :ref:`tracking file<Tracking File>`

            * Entering a list of session IDs (e.g., ``['ses-V1', 'ses-V2']``) will export those sessions only into the dataset
    
    **stage** *(list)*
        * Stages of interest

        * *Acceptable options:*

            * Default is ``None`` which will create datasets for all stages extracted in the :ref:`export event parameters step<Export spindle characteristics>`

            * Entering a list of stages (e.g., ``['NREM3']``) will only export parameters for the events in this specific stage

                .. note::
                    This will only work if the :ref:`export event parameters step<Export spindle characteristics>` was run with ``concat_stage = False``

    **concat_stage** *(logical)*
        * Concatenation options for stages

        * *Acceptable options:*

            * Default is ``False`` which will create datasets per stage (NREM2 vs NREM3) *separately* (saving each stage as a separate ``.csv`` output file)

            * Entering ``True`` will concatenate (i.e., merge) all stages before exporting the parameters of the spindle events

            .. note::
                Pay caution to how the argument ``concat_stage`` was defined in the :ref:`export event parameters step<Export spindle characteristics>` .
                If in this step (**Create datasets**) the argument is set to: ``concat_stage = False``, but in the :ref:`export event parameters step<Export spindle characteristics>` this was set to ``concat_stage = True`` , then this will fail as the spindle events have not been exported for stages combined. The previous step will need to be re-run with ``concat_stage = False``

    **concat_cycle** *(logical)*
        * Concatenation options for sleep cycles

        * *Acceptable options:*

            * Default is ``True`` will create datasets for all cycles concatenated (i.e., merged) in one ``.csv`` dataset file.

            * Entering ``False`` create datasets per sleep cycle *separately* (saving each cycle as a separate ``.csv`` output file)

            .. note::
                Similar to ``concat_stage`` - pay caution to how the argument ``concat_cycle`` was defined in the :ref:`export event parameters step<Export spindle characteristics>` .
                If in this step (**Create datasets**) the argument is set to: ``concat_cycle = False``, but in the :ref:`export event parameters step<Export spindle characteristics>` this was set to ``concat_cycle = True`` , then this will fail as the spindle events have not been exported for cycle combined. The previous step will need to be re-run with ``concat_cycle = False``

    **cycle_idx** *(NoneType or tuple)*
        * Cycles of interest

        * *Acceptable options:*

            * Default is ``None`` which will infer to not take into consideration the cycle and either extract cycle for the whole night if ``concat_cycle = True`` 
            or for all the cycles if ``concat_cycle = False``

            * Entering list of cycle numbers (e.g., ``[1,2,3]``) will extract the spindle parameters for those cycles only. It requires that you have 
            defined ``cycle_idx`` during :ref:`export event parameters<Export spindle characteristics>` and have also set up ``concat_cycle = False``.

    **grp_name** *(str)*
        * Name of the tab in the :ref:`Annotations file` where the detected events are saved 

        * *Acceptable options:*

            * Default is ``eeg`` which is the recommended naming convention
           
            * If entering a list of group names (e.g., ``['eeg_hemiR']``), ensure that this matches ``grp_name`` used in the :ref:`export event parameters step<Export spindle characteristics>`

    **params** *(str or dict)*
        * The names of specific parameters to export into the dataset 

        * *Acceptable options:*

            * Default is ``all`` which will export all characteristics (see :ref:`Output`) -  *Recommended* 

            * To specify only specific parameters to export, enter a `dictionary <https://docs.python.org/3/tutorial/datastructures.html#dictionaries>`_ with ``True`` or ``False`` for each parameter (e.g., ``params = ['dur':True, 'minamp':False, 'maxamp':False, 'ptp':True, 'rms':False, 'power':True, 'peakpf':False, 'energy':False, 'peakef':False]``)

    **outfile** *(str or logical)*
        * Logging of event parameter export

        * *Acceptable options:*

            * Default is ``True`` which will create a logfile *event_dataset_{method}_{datetime}_log.txt* in ``<root_dir>/OUT/audit/logs/``

            * Entering a string ``<custom_outfile_name.txt>`` will save the logfile under that custom name
            
            * Entering ``False`` won't save a logfile

.. hint::
    To combine datasets, use the *trawl* function (see XXXX)


.. _output_spindle:
Outputs of Spindle
----------------

*Parameters of spindle characteristics:*

    **Count** : Number of spindles detected 

    **Density** :  Mean number of spindles detected per period (e.g., 30s, 60s - depend on ``epoch_dur`` argument in *export_eventparams*)

    **Duration_mean** : Mean spindles duration (s)

    **Duration_stdv** : Standard deviation of spindles duration (s)

    **Min_amplitude_mean** : Mean amplitude of the spindles trough (uV)

    **Min_amplitude_stdv** : Standard deviation of the amplitude of the spindles trough (uV)

    **Max_amplitude_mean** : Mean amplitude of the spindles peak (uV)

    **Max_amplitude_stdv** : Standard deviation of the amplitude of the spindles peak (uV)

    **Ptp_amplitude_mean** : Mean peak-to-peak spindles amplitude (uV)

    **Ptp_amplitude_stdv** : Standard deviation of the peak-to-peak spindles amplitude (uV)

    **Power_mean** : Mean absolute spectral power within the ``frequency`` range set in *export_eventparams* (uV2)

    **Power_stdv** : Standard deviation of the absolute spectral power within the ``frequency`` range set in *export_eventparams* (uV2)

    **Peak_power_frequency_mean** : Mean peak power frequency of the spindles events (Hz)

    **Peak_power_frequency_stdv** : Standard deviation of the peak power frequency of the spindles events (Hz)











