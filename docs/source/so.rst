Slow Oscillations
=====

.. _overview:

Overview
------------
Slow oscillations (SOs) are biphasic waves corresponding to the alternation between two stable membranes potential levels (UP states = depolarization and 
DOWN states = hyperpolarization).Oscillating below 1.25 Hz, SOs are generated during NREM2 and NREM3.

| Slow oscillations can be detected as events and their characteristics (see definitions in section :ref:`Output`) can be extracted across NREM (N2+N3), 
per stage and/or per cycle.

| We propose 4 standardized published methods to automatically detect SOs :
    * *Staresina et al. (2015)*: recommended for event detection (<1.25Hz) with amplitude adapted per individual per channel
    Method in brief: 1. Filter the signal (two-pass FIR bandpass filter, 0.5–1.25 Hz, order = 3); 2. A positive-to-negative zero crossing and a subsequent 
    negative-to-positive zero crossing separated by 0.8-2 sec; 3. Top 25% of events with the largest amplitudes for trough-to-peak amplitude between two 
    positive-to-negative zero crossings. `see reference`_.
.. _see reference: https://doi.org/10.1038/nn.4119

    * *Ngo et al. (2015)*: recommended for event detection (<3.5Hz) with amplitude adapted per individual across several channels (average)
    Method in brief: 1. Filter the signal (lowpass filter, 3.5 Hz); 2. A positive-to-negative zero crossing and a subsequent negative-to-positive zero crossing 
    separated by 0.833-2 sec ; 3. A negative peak amplitude lower than 1.25 times the mean negative peak amplitude per subject. `see reference`_.
.. _see reference: https://doi.org/10.1016/j.neuron.2013.03.006

    * *Massimini et al. (2004)*: recommended for detection of SOs with rigid criteria
    Method in brief: 1. Filter the signal (bandpass, 0.1-4 Hz); 2. A positive-to-negative zero crossing and a subsequent negative-to-positive zero crossing 
    separated by 0.3-1 sec; 3. A negative peak between the two zero crossings with voltage less than -80 uV,; 4. A negative-to-positive peak-to-peak 
    amplitude >140 uV. `see reference`_.
.. _see reference: https://doi.org/10.1523/JNEUROSCI.1318-04.2004

    * *Adapted Massimini et al*: recommended for detection of SOs with rigid criteria (based on AASM)
    Method in brief: 1. Filter the signal (bandpass, 0.1-4 Hz); 2. A positive-to-negative zero crossing and a subsequent negative-to-positive zero crossing 
    separated by 0.25-1 sec; 3. A negative peak between the two zero crossings with voltage less than -40 uV,; 4. A negative-to-positive peak-to-peak 
    amplitude >75 uV.


**You will need to run three functions:**

1) Detect SOs events: it will copy the .xml from ``root_dir/OUT/staging/`` to ``root_dir/OUT/slowwave/`` and write events detected 

.. code-block:: python

   project.detect_slow_oscillations()


2) Export event characteristics per method: it will extract a .csv file per channel and/or stage in the subject and session folders in ``root_dir/OUT/slowwave/`` 

.. code-block:: python

   project_name.export_eventparams()
 
3) Create datasets combining all the subjects: it will combine all .csv into a single dataset (one row per subject) per session, stage and channel in ``root_dir/OUT/datasets/``

.. code-block:: python

   project_name.event_dataset()
 

.. _extraction_SO:
Extract slow oscillations
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
        * Method of SOs detection (i.e., Staresina2015, Ngo2015, Massimini2004) 

        * Default is ``['Staresina2015']`` method  

     .. note::
    Only ``['Staresina2015', 'Massimini2004', 'AASM/Massimini2004']`` methods can be run simultaneously. ``['Ngo2015']`` can only be runned separately with ``average_channels = True``

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
                        rater=None, 
                        stage = None, 
                        grp_name = 'eeg', 
                        concat_cycle = True, 
                        concat_stage = False, 
                        cycle_idx = None, 
                        keyword = None, 
                        evt_name = 'spindle', #input required for SO extraction
                        frequency = (11,16),  #input required for SO extraction
                        segs = None, 
                        adap_bands = False, 
                        param_keys = 'all',  
                        epoch_dur = 30, 
                        n_fft_sec = 4, 
                        Ngo = {'run':False},
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

            * If you put a string of sub IDs (e.g., *['sub-01', 'sub-02']*), it will only export the event's characteristics from those sub folders

    **sessions**
        * Sessions/Visits to analyse per subject

        * *Acceptable options:*

            * Default is ``'all'`` which will point to all the *ses* folders within the sub folder in ``root_dir/DATA``

            * If you put ``None``, it will point to the *ses* column in *tracking* file

            * If you put a string of ses visits (e.g., *['ses-V1']*), it will only export the event's characteristics from the selected session(s) within each subject

    **chan**
        * Channel(s) of interest

        * *Acceptable options:*

            * Default is ``None`` which will point to the *chanset* columns in *tracking* file

            * If you put string of channels' names (e.g., *['Cz']*), it will only export the event's characteristics from the selected channels  

    **ref_chan**
        * Reference channel(s) for the channels of interest (e.g., mastoid A1 or A2 or joint mastoids)

        * *Acceptable options:*

            * Default is ``None`` which will point to the *refset* columns in *tracking* file

            * If you put string of channels' names (e.g., *['A1', 'A2']*), it will only re-reference to the channels written 

    **rater**
        * Name of the rater to analyze

        * *Acceptable options:*

            * Default is ``None`` which will discard the name of the rater and expect only one rater per .xml (!! make sure you don't have multiple raters!!)
    
            * If put string of rater's name (e.g., *[Rater1]*), it will only export the the event's characteristics from this rater (and create an empty extraction file if the 
            rater is absent)

    **stage**
        * Stages of interest

        * *Acceptable options:*

            * Default is ``['NREM2', 'NREM3']`` 

            * If you put string of stage (e.g., *['NREM3']*), it will only export the event's characteristics from this specific stage

    **grp_name**
        * Name of the tab in the montage which includes the channels of interest. 

        * *Acceptable options:*

            * Default is ``eeg`` which is the name we recommend
           
            * Need to match whatever was written in *detect_slowocillation*

    **concat_cycle**
        * Concatenation options for sleep cycle

        * *Acceptable options:*

            * Default is ``True`` which means that cycles will be concatenated (i.e., merge) before the exportation of the event's characteristics

            * If you put ``False``, it will export the event's characteristics per cycle

    **concat_stage**
        * Concatenation options for stages

        * *Acceptable options:*

            * Default is ``False`` which means that it will export the event's characteristics per stage (NREM2 and NREM3)

            * If you put ``True``, stages will be concatenated (i.e., merge) before the exportation of the event's characteristics

    **cycle_idx**
        * Sleep cycle numbers

        * *Acceptable options:*

            * Default is ``None`` which will infer no cycles 

            * If you put a list of indices corresponding to sleep cycle numbers (e.g., *(1,2,3,4,5,6,7)*), it will only detect the events for these specific 
            cycles' numbers

    **keyword**
        * Allow search for a filename with a specific wildcard (keyword)

        * *Acceptable options:*

            * Default is ``None`` which will infer no keyword to search for

            * If you put string of keywords, it will only export the event's characteristics from this specific .xml

    **evt_name**
        * Name of the events to export in the .xml (location specified in ``xml_dir``)

        * *Input Required for SO extraction:*

            * Default is ``spindle`` which refer to the Whale spindle detection (will lead to an ERROR argument)

            * Put the name of the method used for *detect_slow_oscillations* (e.g., ``['Staresina2015']``) !! One method per extraction !!

    **frequency**
        * Frequency range of interest

        * *Input Required for SO extraction:*

            * Default is ``(11,16)`` which refer to the default frequency range for spindle extraction

            * Put the frequency range depending on the method used for *detect_slow_oscillations*: Staresina2015 requires ``(0.5,1.25)``; Ngo2015 requires
            ``(0,3.5)``; Massimini2004 and AASM/Massimini2004 requires ``(0.1,4)``



                        segs = None, 
                        adap_bands = False, 
                        param_keys = 'all',  
                        epoch_dur = 30, 
                        n_fft_sec = 4, 
                        Ngo = {'run':False},
                        outfile = True)


    **average_channels**
        * Options to average channels before the detection 

        * Default is ``False``: only pass ``True`` if using the ['Ngo2015'] method

    **outfile**
        * Extraction of output file

        * *Acceptable options:*

            * Default is ``True`` which will create a .xml file per subject and per session in ``root_dir/OUT/slowwave/``
            
            * If put ``False``, it won't extract the .xml file with the events detection








.. _create_datasets:
Create datasets
----------------
*Command line argument:*


.. code-block:: python

   project.event_dataset(chan,
                        xml_dir = None, 
                        out_dir = None, 
                        subs = 'all', 
                        sessions = 'all',  
                        ref_chan = None, 
                        stage = None, 
                        concat_stage = False, 
                        concat_cycle = True, 
                        cycle_idx = None, 
                        grp_name = 'eeg', 
                        adap_bands = 'Fixed', 
                        evt_name = 'spindle', 
                        params = 'all', 
                        outfile=True)


*Positional arguments:*
    **chan**
        * Channel(s) of interest 

        * Required arguments: write a string of channels' names (e.g., *['Fz','Cz']*). Use the names written in the *chanset_rename* column in *tracking* file

    **xml_dir**
        * Path to folder with the .xml file which also contains the .csv extracted with the *detect_slow_oscillations* function

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

            * Default is ``['NREM2', 'NREM3']`` 

            * If you put string of stage (e.g., *['NREM3']*), it will only detect the events for this specific stage
    

    **cycle_idx**
        * Extract sleep macro-architecture per cycle

        * Default is ``None`` which will create a .csv extracting macro-architecture for whole-night only (from light off to light on)
    
            * If put a list of cycle number (e.g., [1,2,3]), it will extract macro-architecture per cycle 
            .. note::
            Make sure the cycles are marked on the .xml in ``root_dir/OUT/staging/``

    **outfile**
        * Extraction of output file

        * Default is ``True`` which will create a .csv dataset file combining all subjects in ``root_dir/OUT/datasets/macro/`` per session
    
            * If put ``False``, it won't extract .csv file 


.. note::
    To combine datasets, use the *trawl* function (see XXXX)


.. _output:
Output
----------------

*Markers of SOs characteristics:*

    **Count** : Number of events detected 




