Slow Oscillations
=====

.. _overview:

Overview
------------
Slow oscillations (SOs) are coherent waves corresponding to the alternation between biphasic membrane potential levels (UP states = depolarization and DOWN states = hyperpolarization). Oscillating below ~1 Hz, SOs are generated during sleep stages NREM2 and NREM3.

| Slow oscillations can be detected as events and their characteristics (see definitions in section :ref:`Output <Outputs of Slow Oscillations>`) can be extracted across NREM (NREM2+NREM3), per stage and/or per cycle.

| Seapipe provides 4 published methods to automatically detect SOs:

    * `Staresina et al. (2015) <https://doi.org/10.1038/nn.4119>`_: SO detection (<1.25Hz) with an adapted amplitude criteria per individual and channel
    
       Method in brief:

        1. Filter the signal (FIR bandpass filter, 0.5–1.25 Hz, order = 3); 

        2. Identify events with a positive-to-negative zero crossing and a subsequent negative-to-positive zero crossing separated by 0.8-2 sec;

        3. Keep the top 25% of events with the largest trough-to-peak amplitudes. 

    * `Ngo et al. (2015) <https://doi.org/10.1016/j.neuron.2013.03.006>`_: Slow wave detection (<3.5Hz) with an adapted amplitude criteria per individual averaged across several channels 
    
        Method in brief: 
        1. Filter the signal (FIR lowpass filter, 3.5 Hz); 
        2. Identify events with a positive-to-negative zero crossing and a subsequent negative-to-positive zero crossing separated by 0.833-2 sec; 
        3. Keep the events with a negative peak amplitude lower than 1.25 times the mean negative peak amplitude per subject.

    * `Massimini et al. (2004) <https://doi.org/10.1523/JNEUROSCI.1318-04.2004>`_: SO detection with a rigid amplitude criteria
    
        Method in brief: 
        1. Filter the signal (FIR bandpass filter, 0.1-4 Hz); 
        2. Identify events with a positive-to-negative zero crossing and a subsequent negative-to-positive zero crossing separated by 0.3-1 sec; 
        3. Keep the events with a negative peak less than -80 uV,; 
        4. Keep the events with a negative-to-positive peak-to-peak amplitude >140 uV.

    * *Adapted Massimini et al*: adapted SO detection with rigid amplitude criteria (lowered based on `AASM criteria for slow wave activity <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5406946/>`_)
   
        Method in brief: 
        1. Same as Massimini et al. 2004, except to keep the events with a negative-to-positive peak-to-peak amplitude >75 uV.

.. admonition:: polarity of recordings
    
    The shape and orientation of slow oscillations depends on the manner of recording, specifically whether they are detected from inside (intracranial EEG, iEEG) or from outside (scalp EEG) the brain.
    
    .. image:: images/polarity.png
        :width: 600
    |    
    When working with scalp EEG, it is also common that recordings are 'inverted' before they are exported.
    The importance of keeping track of the polarity of the EEG data is related to which direction corresponds to the (‘UP’ or 'DOWN'), 
    as this will determine the underlying physiological and biological interpretation! This is especially relevant when running Phase-Amplitude Coupling. 
    Therefore it is recommended that you confirm the polarity of your recordings prior to commencing any analyses.

.. _Functions:
Functions
----------------
| **You will need to run three functions:**

1) Detect SOs events: it will copy the .xml from ``root_dir/OUT/staging/`` to ``root_dir/OUT/slowwave/`` and write events detected 

.. code-block:: python

   project.detect_slow_oscillations()


2) Export event characteristics per method: it will extract a .csv file per channel and/or stage in the subject and session folders in ``root_dir/OUT/slowwave/`` 

.. code-block:: python

   project_name.export_eventparams()
 
3) Create datasets combining all the subjects: it will combine all .csv into a single dataset (one row per subject) per session, stage and channel in ``root_dir/OUT/datasets/``

.. code-block:: python

   project_name.event_dataset()
 

.. _detection_SO:
Detect slow oscillations
----------------
*Command line argument:*

.. code-block:: python

    project.detect_slow_oscillations(xml_dir=None, 
                                    out_dir=None, 
                                    subs='all', 
                                    sessions='all', 
                                    filetype='.edf', 
                                    method = ['Staresina2015'], 
                                    chan=None,
                                    ref_chan=None, 
                                    rater=None, 
                                    grp_name='eeg', 
                                    stage = ['NREM2','NREM3'], 
                                    cycle_idx=None, 
                                    duration=(0.2, 2), 
                                    invert = None,
                                    average_channels = False, 
                                    outfile=True)


*Positional arguments:*

    **xml_dir**
        * Path to folder with the .xml file containing sleep stages and arousal/artefact events. 

        * Default is ``None`` which will point to ``root_dir/OUT/staging``

    **out_dir**
        * Output path for the .xml file containing the new detected event (events will be named like the method used; e.g., Staresina2015)

        * Default is ``None`` which will point to ``root_dir/OUT/slowwave``

    **subs**
        * Subject to analyze

        * *Acceptable options:*

            * Default is ``'all'`` which will point to all the *sub* folders in ``root_dir/DATA``

            * If you put ``None``, it will point to the *sub* column in *tracking* file

            * If you put a string of sub IDs (e.g., *['sub-01', 'sub-02']*), it will only detect those sub folders

    **sessions**
        * Sessions/Visits to analyse per subject

        * *Acceptable options:*

            * Default is ``'all'`` which will point to all the *ses* folders within the sub folder in ``root_dir/DATA``

            * If you put ``None``, it will point to the *ses* column in *tracking* file

            * If you put a string of ses visits (e.g., *['ses-V1']*), it will only detect the selected session(s) within each subject

    **filetype**
        * Format of files containing EEG signal

        * *Acceptable options:*

            * Default is ``'.edf'`` format

            * The pipeline can also read .eeg, .set formats

    **method**
        * Method of SOs detection (i.e., Staresina2015, Ngo2015, Massimini2004,AASM/Massimini2004) 

        * *Acceptable options:*

            * Default is ``['Staresina2015']`` method  
            
            * Only ``['Staresina2015', 'Massimini2004', 'AASM/Massimini2004']`` methods can be run simultaneously. ``['Ngo2015']`` can only be runned separately with ``average_channels = True``

    **chan**
        * Channel(s) of interest

        * *Acceptable options:*

            * Default is ``None`` which will point to the *chanset* columns in *tracking* file

            * If you put string of channels' names (e.g., *['Cz']*), it will only detect the selected channels  

    **ref_chan**
        * Reference channel(s) for the channels of interest (e.g., mastoid A1 or A2 or joint mastoids)

        * *Acceptable options:*

            * Default is ``None`` which will point to the *refset* columns in *tracking* file

            * If you put string of channels' names (e.g., *['A1', 'A2']*), it will only re-reference to the channels written 

    **rater**
        * Name of the rater to analyze

        * *Acceptable options:*

            * Default is ``None`` which will discard the name of the rater and expect only one rater per .xml (!! make sure you don't have multiple raters!!)
    
            * If put string of rater's name (e.g., *[Rater1]*), it will only detects events from this rater per .xml (and create an empty extraction file if the 
            rater is absent)

    **grp_name**
        * Name of the tab in the montage which includes the channels of interest !! It is for visualization in Wonambi only !!

        * *Acceptable options:*

            * Default is ``eeg`` which is the name we recommend
           
            * If you put string of channels' names (e.g., *['eeg_hemiR']*), events can only be seen in Wonambi with a montage that includes a tab with this name

    **stage**
        * Stages of interest

        * *Acceptable options:*

            * Default is ``['NREM2', 'NREM3']`` 

            * If you put string of stage (e.g., *['NREM3']*), it will only detect the events for this specific stage

    **cycle_idx**
        * Sleep cycle numbers

        * *Acceptable options:*

            * Default is ``None`` which will infer no cycles 

            * If you put a list of indices corresponding to sleep cycle numbers (e.g., *(1,2,3,4,5,6,7)*), it will only detect the events for these specific 
            cycles' numbers

    **duration**
        * Minimum and maximum duration of events

        * *Acceptable options:*

            * Default is ``(0.2, 2)`` 

            * If you put a list of 2 indices (e.g., *(0.2,1)*), it will only detect the events with a duration within this range

    **invert**
        * Option to invert polarity

        * *Acceptable options:*

            * Default is ``None`` which will point to the *chanset_invert* columns in *tracking* file. However, if the *tracking* file does not specify *chanset_invert* 
            columns, it will keep the polarity of the recording as it is 

            * If you put ``False``, it will keep the polarity of the recording as it is

            * If you put ``True``, it will reverse the polarity of the recording 

    **average_channels**
        * Options to average channels before the detection 

        * Default is ``False``: only pass ``True`` if using the ['Ngo2015'] method

    **outfile**
        * Extraction of output file

        * *Acceptable options:*

            * Default is ``True`` which will create a .xml file per subject and per session in ``root_dir/OUT/slowwave/``
            
            * If put ``False``, it won't extract the .xml file with the events detection


.. _export_SO:
Export slow oscillations characteristics
----------------
*Command line argument:*
To run per method if usin multiple detection methods

.. code-block:: python

    project.export_eventparams(xml_dir = None, 
                        out_dir = None, 
                        subs = 'all', 
                        sessions = 'all', 
                        chan = None, 
                        ref_chan = None, 
                        stage = ['NREM2','NREM3'], 
                        grp_name = 'eeg',
                        rater=None, 
                        cycle_idx = None, 
                        concat_cycle = True, 
                        concat_stage = False, 
                        keyword = None, 
                        segs = None,
                        evt_name = 'spindle', #input required
                        frequency = None,  #input required
                        params = 'all',  
                        epoch_dur = 30, 
                        average_channels = False,
                        outfile = True)


*Positional arguments:*

    **xml_dir**
        * Path to folder with the .xml file containing sleep stages, arousal/artefact events and newly detected slow oscillations events.

        * Default is ``None`` which will point to ``root_dir/OUT/slowwave``

    **out_dir**
        * Output path for the created .csv file containing the characteristics of the slow oscillation events per subject, session, stage, channel

        * Default is ``None`` which will point to ``root_dir/OUT/slowwave``

    **subs**
        * Subject to analyze

        * *Acceptable options:*

            * Default is ``'all'`` which will point to all the *sub* folders in ``root_dir/DATA``

            * If you put ``None``, it will point to the *sub* column in *tracking* file

            * If you put a string of sub IDs (e.g., *['sub-01', 'sub-02']*), it will only export the SOs' characteristics from those sub folders

    **sessions**
        * Sessions/Visits to analyse per subject

        * *Acceptable options:*

            * Default is ``'all'`` which will point to all the *ses* folders within the sub folder in ``root_dir/DATA``

            * If you put ``None``, it will point to the *ses* column in *tracking* file

            * If you put a string of ses visits (e.g., *['ses-V1']*), it will only export the SOs' characteristics from the selected session(s) within each subject

    **chan**
        * Channel(s) of interest

        * *Acceptable options:*

            * Default is ``None`` which will point to the *chanset* columns in *tracking* file - *Recommended*

            * If you put string of channels' names (e.g., *['Cz']*), it will only export the SOs' characteristics from the selected channels  

    **ref_chan**
        * Reference channel(s) for the channels of interest (e.g., mastoid A1 or A2 or joint mastoids)

        * *Acceptable options:*

            * Default is ``None`` which will point to the *refset* columns in *tracking* file - *Recommended*

            * If you put string of channels' names (e.g., *['A1', 'A2']*), it will only export the SOs' characteristics from the selected channels and reference written

    **stage**
        * Stages of interest

        * *Acceptable options:*

            * Default is ``['NREM2', 'NREM3']`` 

            * If you put string of stage (e.g., *['NREM3']*), it will only export the SOs' characteristics from this specific stage

    **grp_name**
        * Name of the tab in the montage which includes the channels of interest. 

        * *Acceptable options:*

            * Default is ``eeg`` which is the name we recommend
           
            * Need to match ``grp_name`` used in *detect_slowocillation*

    **rater**
        * Name of the rater to analyze

        * *Acceptable options:*

            * Default is ``None`` which will discard the name of the rater and expect only one rater per .xml (!! make sure you don't have multiple raters!!)
    
            * If put string of rater's name (e.g., *[Rater1]*), it will only export the the event's characteristics from this rater (and create an empty extraction file if the 
            rater is absent)

    **cycle_idx**
        * Sleep cycle numbers

        * *Acceptable options:*

            * Default is ``None`` which will infer no cycle

            * If you put a list of indices corresponding to sleep cycle numbers (e.g., *(1,2)*), it will only export the SOs' characteristics from these 
            specific cycles. Also requires ``concat_cycle = False``

    **concat_cycle**
        * Concatenation options for sleep cycle

        * *Acceptable options:*

            * Default is ``True`` which means that cycles will be concatenated (i.e., merge) before the exportation of the SOs' characteristics

            * If you put ``False``, it will export SOs' characteristics per cycle

    **concat_stage**
        * Concatenation options for stages

        * *Acceptable options:*

            * Default is ``False`` which means that it will export SOs' characteristics per stage (NREM2 vs NREM3)

            * If you put ``True``, stages will be concatenated (i.e., merge) before the exportation of SOs' characteristics

    **keyword**
        * Allow search for a filename with a specific wildcard (keyword)

        * *Acceptable options:*

            * Default is ``None`` which will infer no keyword to search for

            * If you put string of keywords, it will only export the event's characteristics from this specific .xml

    **seg**
        * Option to extract parameters between certain markers, which need to be defined in the .xml file in ``root_dir/OUT/staging``

        * *Acceptable options:*

            * Default is ``None`` which will infer no segmentation

            * If you put a list of tuples, with both tags named (e.g. *[('N2_ON','N2_OFF'), ('N3_ON','N3_OFF')]*), it will only export the event's characteristics within the events markers (segments)

    **evt_name**
        * Name of the event of interest to export from the .xml 

        * *Input Required for SO extraction:*

            * Default is ``spindle`` which refer to the Whale spindle detection (will lead to an ERROR argument)

            * Put the name of the method used for *detect_slow_oscillations* (e.g., ``['Staresina2015']``) !! One method per extraction !!

    **frequency**
        * Frequency range of interest

        * *Input Required:*

            * Put the frequency range depending on the method used for *detect_slow_oscillations*: Staresina2015 requires ``(0.5,1.25)``; Ngo2015 requires
            ``(0,3.5)``; Massimini2004 and AASM/Massimini2004 requires ``(0.1,4)``

    **params**
        * Options to export specific characteristics only

        * *Acceptable options:*

            * Default is ``all`` which will export all characteristics (see :ref:`Output`) -  *Recommended*

            * You can specify characteristics of interest using ``True/False`` arguments (e.g., ``params = ['dur':True, 'minamp':False, 'maxamp':False, 'ptp':True, 'rms':False, 'power':True, 'peakpf':False, 
                         'energy':False, 'peakef':False]``)

    **epoch_dur**
        * Options to change the denominator (duration for index density)

        * *Acceptable options:*

            * Default is ``30`` infers 30-seconds epoch

            * If you put a number (e.g., *60*), it will use that number as denominator for the computation of SO density

    **average_channels**
        * Refer to the options to average channels before the detection - only relevant if you used the ``['Ngo2015']`` method in *detect_slow_oscillations*

        * Default is ``False``: only pass ``True`` if used the ``['Ngo2015']`` method to detect SOs

    **outfile**
        * Extraction of output file

        * *Acceptable options:*

            * Default is ``True`` which will create a .csv file per subject, session, channel, stage in ``root_dir/OUT/slowwave/``
            
            * If put ``False``, it won't extract the .csv file with the events' characteristics


     .. note::
        By default
        * - *export_eventparams* cannot extract SOs characteristics without required arguments for ``evt_name`` and ``frequency``. 

        * - it will extract characteristics per stage (NREM2 vs NREM3). If you want the extraction for NREM2+NREM3 combined as well, re-run *export_eventparams* 
        with ``concat_stage = True``.

        * - it will extract characteristics for the whole-night. If you want the extraction per cycle and per stage as well, re-run *export_eventparams* 
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


*Positional arguments:*
    **chan**
        * Channel(s) of interest

        * *Input Required:*
        
            * Write a string of channels' names (e.g., *['Fz','Cz', 'Pz']*). Use the names written in the *chanset_rename* columns in *tracking* file

    **xml_dir**
        * Path to folder containing the .csv extracted with the *export_eventparams* function

        * Default is ``None`` which will point to ``root_dir/OUT/slowwave/``

    **out_dir**
        * Output path for the created datasets

        * Default is ``None`` which will point to ``root_dir/OUT/datasets/``

    **subs**
        * Subject to export in the datasets

        * Default is ``'all'`` which will point to all the *sub* folders in ``root_dir/OUT/staging``

            * If put ``None``, it will point to the *sub* column in *tracking* file

            * If put list of sub ID (e.g., *['sub-01', 'sub-02']*), it will only detect those sub folders

    **sessions**
        * Sessions/Visits to extract per subject

        * Default is ``'all'`` which will point to all the *ses* folders within the sub folder in ``root_dir/OUT/staging``

            * If put ``None``, it will point to the *ses* column in *tracking* file

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


.. _output:
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











