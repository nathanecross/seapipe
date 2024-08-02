Slow Oscillations
=====

.. _overview:

Overview
------------
Slow oscillations (SOs) are coherent waves corresponding to the alternation between biphasic membrane potential levels (UP states = depolarization 
and DOWN states = hyperpolarization). Oscillating below ~1 Hz, SOs are generated during sleep stages NREM2 and NREM3.

| Slow oscillations can be detected as events and their characteristics (see definitions in section :ref:`Output<Outputs of Slow Oscillations>`) can be extracted across NREM (NREM2+NREM3), per stage and/or per cycle.

| Seapipe provides 4 published methods to automatically detect SOs:

    * `Staresina et al. (2015) <https://doi.org/10.1038/nn.4119>`_: SO detection (<1.25Hz) with an adapted amplitude criteria per individual and channel.
    
       Method in brief:

        1. Filter the signal (FIR bandpass filter, 0.5–1.25 Hz, order = 3); 

        2. Identify events with a positive-to-negative zero crossing and a subsequent negative-to-positive zero crossing separated by 0.8-2 sec;

        3. Keep the top 25% of events with the largest trough-to-peak amplitudes. 

    * `Ngo et al. (2015) <https://doi.org/10.1016/j.neuron.2013.03.006>`_: Slow wave detection (<3.5Hz) with an adapted amplitude criteria per individual averaged across several channels. 
    
        Method in brief: 

        1. Filter the signal (FIR lowpass filter, 3.5 Hz); 

        2. Identify events with a positive-to-negative zero crossing and a subsequent negative-to-positive zero crossing separated by 0.833-2 sec; 

        3. Keep the events with a negative peak amplitude lower than 1.25 times the mean negative peak amplitude per subject.

    * `Massimini et al. (2004) <https://doi.org/10.1523/JNEUROSCI.1318-04.2004>`_: SO detection with a rigid amplitude criteria.
    
        Method in brief: 

        1. Filter the signal (FIR bandpass filter, 0.1-4 Hz); 

        2. Identify events with a positive-to-negative zero crossing and a subsequent negative-to-positive zero crossing separated by 0.3-1 sec; 

        3. Keep the events with a negative peak less than -80 uV,; 

        4. Keep the events with a negative-to-positive peak-to-peak amplitude >140 uV.

    * *Adapted Massimini et al*: adapted SO detection with rigid amplitude criteria (lowered based on  `AASM criteria for slow wave activity <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5406946/>`_).
   
        Method in brief:

        1. Same as Massimini et al. 2004, except to keep the events with a negative-to-positive peak-to-peak amplitude >75 uV.

.. admonition:: polarity of recordings
    
    The shape and orientation of slow oscillations depends on the manner of recording, specifically whether they are detected from inside (intracranial EEG, iEEG) or from outside (scalp EEG) the brain.
    
    .. image:: images/polarity.png
        :width: 700
    |    
    When working with scalp EEG, it is also common that recordings are 'inverted' before they are exported.
    The importance of keeping track of the polarity of the EEG data is related to which direction corresponds to the (‘UP’ or 'DOWN'), 
    as this will determine the underlying physiological and biological interpretation! This is especially relevant when running Phase-Amplitude Coupling. 
    Therefore it is recommended that you confirm the polarity of your recordings prior to commencing any analyses.

.. _Functions:
Functions to detect Slow Oscillations
----------------
| **Detecting Slow Oscillations and extracting their parameters will involve three functions:**

1) Detect slow oscillation events:  

.. code-block:: python

   project.detect_slow_oscillations()
|
    This will copy the :ref:`Annotations file` from every ``/sub-XXX/ses-XXX`` in ``<xml_dir>`` to ``<root_dir>/OUT/slowwave/`` and write in the detected events. 
|
2) Export event characteristics: 

.. code-block:: python

   project_name.export_eventparams()
|   
    This will extract a ``.csv`` file for every channel and/or stage and/or cycle into the ``/sub-XXX/ses-XXX`` directory in ``<root_dir>/OUT/slowwave/`` 
|
3) Create datasets combining all the subjects: 

.. code-block:: python

   project_name.event_dataset()
|
    This will combine all of the ``.csv`` files from the previous step into a single dataset (one row per subject) ``<root_dir>/OUT/datasets/``
|

.. _detection_SO:
Detect slow oscillations
----------------
*Command line argument:*

.. code-block:: python

    project.detect_slow_oscillations(xml_dir = None, 
                                     out_dir = None, 
                                     subs = 'all', 
                                     sessions = 'all', 
                                     filetype = '.edf', 
                                     method = ['Staresina2015'], 
                                     chan = None,
                                     ref_chan = None, 
                                     rater = None, 
                                     grp_name = 'eeg', 
                                     stage = ['NREM2','NREM3'], 
                                     cycle_idx = None, 
                                     duration = (0.2, 2), 
                                     invert = None,
                                     average_channels = False, 
                                     outfile = True)


*Positional arguments:*

    **xml_dir** *(str)*
        * Path to the directory with sub-directories ``/sub-XXX/ses-XXX`` containing the input :ref:`Annotations files<Annotations file>`. 

        * Default is ``None`` which will point to ``<root_dir>/OUT/staging/`` (Annotations files with sleep stage markings and arousal/artefact events).

    **out_dir** *(str)*
        * Output path for the .xml file containing the new detected event (events will be named like the method used; e.g., ``Staresina2015``)

        * Default is ``None`` which will point to ``<root_dir>/OUT/slowwave/``

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
        * Method for SO detection

        * *Acceptable options:*

            * Default is ``['Staresina2015']`` method  
            
            * Only ``['Staresina2015', 'Massimini2004', 'AASM/Massimini2004']`` methods can be run simultaneously. ``['Ngo2015']`` can only be runned separately with ``average_channels = True``

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

.. admonition:: NOTE

    This assumes there is only one rater per Annotations file (``.xml``) 
    !! make sure you don't have multiple raters!!
    
            * Entering a list of rater names (e.g., ``[<Rater1>, <Rater2>]``) will only save detected events on this rater in the Annotations file
|

    **grp_name** *(str)*
        * Name of the tab in the :ref:`Annotations file` to save the detections to. This is for visualization in Wonambi only, however it will impact the `exporting <Export slow oscillations characteristics>` of events also

        * *Acceptable options:*

            * Default is ``eeg`` which is the recommended naming convention
           
            * Entering a list of group names (e.g., ``['eeg_hemiR']``) will save the events to a tab of this name in the Annotations file. The events can only be visualised in :ref:`Wonambi` with a montage that includes a tab with this name

    **stage** *(list)*
        * Stages of interest

        * *Acceptable options:*

            * Default is ``['NREM2', 'NREM3']`` 

            * Entering a list of stages (e.g., ``['NREM3']``), it will only detect the events for this specific stage. **It is recommended that you leave the default option**

    **cycle_idx** *(NoneType or tuple)*
        * Sleep cycle numbers

        * *Acceptable options:*

            * Default is ``None`` which will infer no cycles 

            * Entering a `tuple <https://docs.python.org/3/tutorial/datastructures.html#tuples-and-sequences>`_ of integers corresponding to sleep cycle numbers (e.g., ``(1,2,3,4,5)``), it will only detect the events for these specific 
            cycles' numbers. If a ``sub`` has less than the number of cycles entered, then the maximum number of cycles possible will be used for that subject.

    **duration** *(tuple)*
        * Minimum and maximum duration of events that will be detected. Any events with durations that are outside these limits will be discarded

        * *Acceptable options:*

            * Default is ``(0.2, 2)`` 

            * Entering a `tuple <https://docs.python.org/3/tutorial/datastructures.html#tuples-and-sequences>`_ of float with length 2 (e.g., ``(0.5, 1)``)  will limit the detection to events with a duration within this range

    **invert** *(NoneType or logical)*
        * Option to invert polarity

        * *Acceptable options:*

            * Default is ``None`` which will point to the *chanset_invert* columns in the :ref:`tracking file<Tracking File>`. However, if the *tracking* file does not specify *chanset_invert* 
            columns, the detection will default to ``False``

            * Entering ``False`` will keep the polarity of the recording as it is

            * Entering ``True`` will reverse (flip) the polarity of the recording 

    **average_channels** *(logical)*
        * Options to average multiple channels before the detection 

        * Default is ``False``
        
        * Only pass ``True`` if using the ``['Ngo2015']`` method

    **outfile** *(str or logical)*
        * Logging of detection

        * *Acceptable options:*

            * Default is ``True`` which will create a logfile *detect_slowosc_{method}_{datetime}_log.txt* in ``<root_dir>/OUT/audit/logs/``

            * Entering a string ``<custom_outfile_name.txt>`` will save the logfile under that custom name
            
            * Entering ``False`` won't save a logfile


.. _export_SO:
Export slow oscillations characteristics
----------------
*Command line argument:*
To run per method if usin multiple detection methods

.. code-block:: python

    project.export_eventparams(evt_name,
                               frequency,
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
                               params = 'all',  
                               epoch_dur = 30, 
                               average_channels = False,
                               outfile = True)

*Required arguments:*

    **evt_name** *(list or str)*
        * Name of the event of interest to export from the :ref:`Annotations file` 

        * Enter a string (e.g ``Staresina2015``) which refers to the event as it is named in the Annotations file. **NOTE** This will be the method name used in the :ref:`detection<Detect slow oscillations>`

        * Entering a list of event names (e.g ``['Staresina2015', 'Massimini2004']``) will export the parameters for each event *separately*

    **frequency** (tuple)
        * Frequency range of interest in which to export certain event parameters (e.g. *frequency*, *power*)

        * Enter a `tuple <https://docs.python.org/3/tutorial/datastructures.html#tuples-and-sequences>`_ containing the frequency range depending on the method used in the :ref:`detection<Detect slow oscillations>`: 
            - Staresina2015 requires ``(0.5, 1.25)``
            - Ngo2015 requires ``(0, 3.5)``
            - Massimini2004 and AASM/Massimini2004 require ``(0.1, 4)``

*Positional arguments:*

    **xml_dir** *(str)*
        * Path to the directory with sub-directories ``/sub-XXX/ses-XXX`` containing the :ref:`Annotations files<Annotations file>` where the :ref:`detections<Detect slow oscillations>` were saved. 

            * Default is ``None`` which will point to ``<root_dir>/OUT/slowwave/``

    **out_dir** *(str)*
        * Output path for the where to save the ``.csv`` file containing the parameters of the slow oscillation events per subject, session, and/or stage, and/or channel.

            * Default is ``None`` which will point to ``root_dir/OUT/slowwave/``

    **subs** *(str, NoneType or list)*
        * Subject IDs to analyze

        * *Acceptable options:*

            * Default is ``'all'`` which will point to all the ``sub-XXX/`` directories in ``<root_dir>/DATA/``

            * Entering ``None`` will point seapipe to the *sub* column in the :ref:`tracking file<Tracking File>`

            * Entering a list of sub IDs (e.g., ``['sub-01', 'sub-02']``) will the parameters of the slow oscillation events for those subjects only

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

.. admonition:: NOTE2

    If the reference channels are not the same as were entered in the :ref:`detection stage<Detect slow oscillations>`, the event parameters will still be stored,
    however the parameters (e.g. frequency, amplitude, power) might be affected. Be careful to remain consistent across these stages!           
|

    **stage** *(list)*
        * Stages of interest

        * *Acceptable options:*

            * Default is ``['NREM2', 'NREM3']`` 

            * Entering a list of stages (e.g., ``['NREM3']``), it will only export parameters for the events in this specific stage

    **grp_name** *(str)*
        * Name of the tab in the :ref:`Annotations file` where the detected events are saved 

        * *Acceptable options:*

            * Default is ``eeg`` which is the recommended naming convention
           
            * If entering a list of group names (e.g., ``['eeg_hemiR']``), ensure that this matches ``grp_name`` used in the :ref:`detection stage<Detect slow oscillations>`

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

            * Default is ``True`` which means that all cycles will be concatenated (i.e., merged) before exporting the parameters of the SO events

            * Entering ``False`` will export the SO parameters per sleep cycle (saving each cycle as a separate ``.csv`` output file)

    **concat_stage** *(logical)*
        * Concatenation options for stages

        * *Acceptable options:*

            * Default is ``False`` which means the parameters of SO events will be exported per stage (e.g. NREM2, NREM3) separately (saving each stage as a separate ``.csv`` output file)

            * Entering ``True`` will concatenate (i.e., merge) all stages before exporting the parameters of the SO events

    **keyword** *(str)*
        * Allow seapipe to search for a Annotations filename containing a specific wildcard (keyword)

        * *Acceptable options:*

            * Default is ``None`` which will infer no keyword to search for.

            * Entering a string (e.g. ``Staresina_adapted_custom``) will only export event parameters from this specific Annotations file

    **seg** *(NoneType or list of tuples)*
        * Option to extract parameters of SOs that only occur in between certain markers. These markers need to be events saved in the :ref:`Annotations file`

        * *Acceptable options:*

            * Default is ``None`` which will infer no segmentation prior to exporting event parameters

            * Entering a list of `tuples <https://docs.python.org/3/tutorial/datastructures.html#tuples-and-sequences>`_, with both start and end tags named (e.g. ``[('N2_ON', 'N2_OFF'), ('N3_ON', 'N3_OFF')]``) will export event parameters that only occur between these event markers

    **params** *(str or dict)*
        * The names of specific parameters to export

        * *Acceptable options:*

            * Default is ``all`` which will export all characteristics (see :ref:`Output`) -  *Recommended*

            * To specify only specific parameters to export, enter a `dictionary <https://docs.python.org/3/tutorial/datastructures.html#dictionaries>`_ with ``True`` or ``False`` for each parameter (e.g., ``params = ['dur':True, 'minamp':False, 'maxamp':False, 'ptp':True, 'rms':False, 'power':True, 'peakpf':False, 
                         'energy':False, 'peakef':False]``)

    **epoch_dur** *(int)*
        * Options to change the denominator (duration) for the *SO density* index 

        * *Acceptable options:*

            * Default is ``30`` (this infers 30-second epochs)

            * Entering a number (e.g., ``60``) will imply that the SO density value equals the number of events per this time period (e.g. per *60 seconds*)
    
    **average_channels** *(logical)*
        * Options to average multiple channels before the detection (Refer to this option in :ref:`Detect slow oscillations`)

        * Default is ``False``
        
        * Only pass ``True`` if using the ``['Ngo2015']`` method

    **outfile** *(str or logical)*
        * Logging of event parameter export

        * *Acceptable options:*

            * Default is ``True`` which will create a logfile *export_params_{method}_{datetime}_log.txt* in ``<root_dir>/OUT/audit/logs/``

            * Entering a string ``<custom_outfile_name.txt>`` will save the logfile under that custom name
            
            * Entering ``False`` won't save a logfile

     .. note::
        By default
        * - *export_eventparams* cannot extract SOs characteristics without required arguments for ``evt_name`` and ``frequency``. 

        * - *export_eventparams* will extract characteristics per stage (NREM2 vs NREM3). If you want the extraction for NREM2+NREM3 combined as well, re-run *export_eventparams* 
        with ``concat_stage = True``.

        * - *export_eventparams* will extract characteristics for the whole-night. If you want the extraction per cycle and per stage as well, re-run *export_eventparams* 
        with ``concat_cycle = False`` and ``concat_stage = False``.



.. _create_datasets:
Create datasets
----------------
*Command line argument:*


.. code-block:: python

   project.event_dataset(chan, #input required
                        xml_dir = None, 
                        out_dir = None, 
                        subs = 'all', 
                        sessions = 'all',  
                        stage = None, 
                        concat_stage = False, 
                        concat_cycle = True, 
                        cycle_idx = None, 
                        grp_name = 'eeg', 
                        evt_name = 'spindle', #input required 
                        params = 'all', 
                        outfile=True)


*Required arguments:*
    **chan**
        * Channel(s) of interest

        * *Input Required:*
        
            * Write a string of channels' names (e.g., *['Fz','Cz', 'Pz']*). Use the names written in the *chanset_rename* columns in the :ref:`tracking file<Tracking File>`

*Positional arguments:*
    **xml_dir**
        * Path to folder containing the .csv extracted with the *export_eventparams* function

        * Default is ``None`` which will point to ``root_dir/OUT/slowwave/``

    **out_dir**
        * Output path for the created datasets

        * Default is ``None`` which will point to ``root_dir/OUT/datasets/``

    **subs**
        * Subject to export in the datasets

        * Default is ``'all'`` which will point to all the *sub* folders in ``root_dir/OUT/staging``

            * If put ``None``, it will point to the *sub* column in the :ref:`tracking file<Tracking File>`

            * If put list of sub ID (e.g., *['sub-01', 'sub-02']*), it will only detect those sub folders

    **sessions**
        * Sessions/Visits to extract per subject

        * Default is ``'all'`` which will point to all the *ses* folders within the sub folder in ``root_dir/OUT/staging``

            * If put ``None``, it will point to the *ses* column in the :ref:`tracking file<Tracking File>`

            * If put string of ses visit (e.g., *['ses-V1']*), it will only detect that/these session(s) within each subject

    **stage**
        * Stages of interest

        * *Acceptable options:*

            * Default is ``None`` which will create datasets for all stages extracted with the *export_eventparams* function

            * If you put string of stage (e.g., *['NREM3']*), it will only export the SOs' characteristics from this specific stage (if you 
            runmed *export_eventparams* with ``concat_stage = False``)

    **concat_stage**
        * Concatenation options for stages

        * *Acceptable options:*

            * Default is ``False`` which means that it will create datasets per stage (NREM2 vs NREM3). It requires that you have runned *export_eventparams* 
            with ``concat_stage = False``.

            * If you put ``True``, it will create datasets "whole_night" combining NREM2+NREM3. It requires that you have runned *export_eventparams* 
            with ``concat_stage = True``.

    **concat_cycle**
        * Concatenation options for sleep cycle

        * *Acceptable options:*

            * Default is ``True`` which means  that it will create datasets "whole_night" combining all cycles. It requires that you have runned *export_eventparams* 
            with ``concat_cycle = True``.

            * If you put ``False``, it will create datasets per cycle. It requires that you have runned *export_eventparams* with ``concat_cycle = False``.

    **cycle_idx**
        * Cycles of interest

        * *Acceptable options:*

            * Default is ``None`` which will infer to not take into consideration the cycle and either extract cycle for the whole night if ``concat_cycle = True`` 
            or for all the cycles if ``concat_cycle = False``

            * If put a list of cycle number (e.g., [1,2,3]), it will extract the SOs' characteristics for those cycles only. It requires that you have 
            define ``cycle_idx`` during *export_eventparams* and have also set up ``concat_cycle = False``.

    **grp_name**
        * Name of the tab in the montage which includes the channels of interest. 

        * *Acceptable options:*

            * Default is ``eeg`` which is the name we recommend
           
            * Need to match whatever was written in *detect_slowocillation* and *export_eventparams*

    **evt_name**
        * Name of the events of interest 

        * *Input Required for SO extraction:*

            * Default is ``spindle`` which refer to the Whale spindle detection (will lead to an ERROR argument)

            * Put the name of the method used for *detect_slow_oscillations* and *export_eventparams* (e.g., ``['Staresina2015']``) !! One method per extraction !!

    **params**
        * Options to create dataset with specific characteristics only

        * *Acceptable options:*

            * Default is ``all`` which will export all characteristics (see :ref:`Output`) -  *Recommended*

            * You can specify characteristics of interest using ``True/False`` arguments (e.g., ``params = ['dur':True, 'minamp':False, 'maxamp':False, 'ptp':True, 'rms':False, 'power':True, 'peakpf':False, 
                         'energy':False, 'peakef':False]``)

    **outfile**
        * Extraction of output file

        * Default is ``True`` which will create a .csv dataset file combining all subjects in ``root_dir/OUT/datasets/evt_name`` per session and per channel
    
            * If put ``False``, it won't extract .csv file 


.. hint::
    To combine datasets, use the *trawl* function (see XXXX)


.. _output_so:
Outputs of Slow Oscillations
----------------

*Parameters of SOs characteristics:*

    **Count** : Number of SOs detected 

    **Density** :  Mean number of SOs detected per period (e.g., 30s, 60s - depend on ``epoch_dur`` argument in *export_eventparams*)

    **Duration_mean** : Mean SOs duration (s)

    **Duration_stdv** : Standard deviation of SOs duration (s)

    **Min_amplitude_mean** : Mean amplitude of the SOs trough (uV)

    **Min_amplitude_stdv** : Standard deviation of the amplitude of the SOs trough (uV)

    **Max_amplitude_mean** : Mean amplitude of the SOs peak (uV)

    **Max_amplitude_stdv** : Standard deviation of the amplitude of the SOs peak (uV)

    **Ptp_amplitude_mean** : Mean peak-to-peak SOs amplitude (uV)

    **Ptp_amplitude_stdv** : Standard deviation of the peak-to-peak SOs amplitude (uV)

    **Power_mean** : Mean absolute spectral power within the ``frequency`` range set in *export_eventparams* (uV2)

    **Power_stdv** : Standard deviation of the absolute spectral power within the ``frequency`` range set in *export_eventparams* (uV2)

    **Peak_power_frequency_mean** : Mean peak power frequency of the SO events (Hz)

    **Peak_power_frequency_stdv** : Standard deviation of the peak power frequency of the SO events (Hz)











