Slow Oscillations
=====

.. _overview:

Overview
------------
Slow oscillations (SOs) are biphasic waves corresponding to the alternation between two stable membranes potential levels (UP states = depolarization and 
DOWN states = hyperpolarization).Oscillating below 1.25 Hz, SOs are generated during NREM2 and NREM3.

| Slow oscillations can be detected as events and their characteristics (see definitions in section :ref:`Output`) can be extracted across NREM (N2+N3), 
per stage and/or per cycle.

| We propose 3 standardized published methods to automatically detect SOs :
    * *Staresina et al. (2015)*: recommended for event detection (<1.25Hz) with amplitude adapted per individual per channel
    Method in brief: 1. Filter the signal (two-pass FIR bandpass filter, 0.5–1.25 Hz, order = 3); 2. A positive-to-negative zero crossing and a subsequent 
    negative-to-positive zero crossing separated by 0.8-2.0 sec; 3. Top 25% of events with the largest amplitudes for trough-to-peak amplitude between two 
    positive-to-negative zero crossings. `see reference`_.
.. _see reference: https://doi.org/10.1038/nn.4119

    * *Ngo et al. (2013)*: recommended for event detection (<3.5Hz) with amplitude adapted per individual across several channels (average)
    Method in brief: 1. Filter the signal (lowpass filter, 3.5 Hz); 2. A positive-to-negative zero crossing and a subsequent negative-to-positive zero crossing 
    separated by 0.833-2.0 sec ; 3. A negative peak amplitude lower than 1.25 times the mean negative peak amplitude per subject. `see reference`_.
.. _see reference: https://doi.org/10.1016/j.neuron.2013.03.006

    * *Massimini et al. (2004)*: recommended for detection of textbook SOs (75uV peak-to-peak, 1-4Hz)
    Method in brief: 1. Filter the signal (bandpass, 0.1-4 Hz); 2. A positive-to-negative zero crossing and a subsequent negative-to-positive zero crossing 
    separated by 0.3-1.0 sec; 3. A negative peak between the two zero crossings with voltage less than -80 uV,; 4. A negative-to-positive peak-to-peak 
    amplitude >140 uV. `see reference`_.
.. _see reference: https://doi.org/10.1523/JNEUROSCI.1318-04.2004


**You will need to run three functions:**

1) Detect SOs events
    * It will copy the .xml from ``root_dir/OUT/staging/`` to ``root_dir/OUT/slowwave/`` and write events detected 

.. code-block:: python

   project.detect_slow_oscillations()


2) Export event characteristics 
    * It will extract a .csv file per channel and/or stage in the subject and session folders in ``root_dir/OUT/slowwave/`` 

.. code-block:: python

   project_name.export_eventparams()
 
3) Create datasets combining all the subjects
   * It will combine all .csv into a single dataset (one row per subject) per session, stage and channel in ``root_dir/OUT/datasets/``

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
        * Path to folder with the .xml file containing sleep stages and arousal events. 

        * Default is ``None`` which will point to ``root_dir/OUT/staging``

    **out_dir**
        * Output path for the .xml file containing the new detected event (named like the method used; e.g., Staresina2015)

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
    Only ``['Staresina2015', 'Massimini2004']`` methods can be run simultaneously. ``['Ngo2015']`` can only be runned separately with ``average_channels = True``

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


.. _create_datasets:
Create datasets
----------------
*Command line argument:*

.. code-block:: python

   project_name.macro_dataset(xml_dir = None, 
                              out_dir = None, 
                              subs = 'all', 
                              sessions = 'all', 
                              cycle_idx = None,
                              outfile = True)


*Positional arguments:*

    **xml_dir**
        * Path to folder with the .xml file which also contains the .csv extracted with the *export_macro_stats* function

        * Default is ``None`` which will point to ``root_dir/OUT/staging``

    **out_dir**
        * Output path for the created datasets

        * Default is ``None`` which will point to ``root_dir/OUT/datasets/macro/``

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

    **cycle_idx**
        * Extract sleep macro-architecture per cycle

        * Default is ``None`` which will create a .csv extracting macro-architecture for whole-night only (from light off to light on)
    
            * If put a list of cycle number (e.g., [1,2,3]), it will extract macro-architecture per cycle *!!! Make sure you marked the cycles on the .xml in staging (see wonambi)!!!*

    **outfile**
        * Extraction of output file

        * Default is ``True`` which will create a .csv dataset file combining all subjects in ``root_dir/OUT/datasets/macro/`` per session
    
            * If put ``False``, it won't extract .csv file 


.. note::
    To combine datasets, use the *trawl* function (see XXXX)


.. _output:
Output
----------------

*Markers of macro-architecture:*

    **TIB_min** : time in bed from light off to light on - in minutes

    **TotalWake_min** : total wake duration between light off and light on (including SL, WASO, Wmor) - in minutes

    **SL_min** : sleep onset latency from light off to first epoch of sleep - in minutes

    **WASOintra_min** : wake after sleep onset (wake duration from SOL to last epoch of sleep) - in minutes

    **Wmor_min** : wake duration from last epoch of sleep to light on - in minutes

    **TSP_min** : total sleep period (duration from SOL to last epoch of sleep, includes epochs of N1, N2, N3, REM and Wake) - in minutes

    **TST_min** : total sleep time (only includes epochs of N1, N2, N3, REM) - in minutes

    **SE_%** : sleep efficiency (TST/TiB*100) - in percentage

    **N1_min** : time spent in stage N1 - in minutes

    **N2_min** : time spent in stage N2 - in minutes

    **N3_min** : time spent in stage N3 - in minutes

    **REM_min** : time spent in stage REM - in minutes

    **W_%tsp** : proportion of time spent in wake relative to TSP (WASO_intra/TSP*100) - in percentage

    **N1_%tsp** : proportion of time spent in N1 relative to TSP (N1/TSP*100) - in percentage

    **N2_%tsp** : proportion of time spent in N2 relative to TSP (N2/TSP*100) - in percentage

    **N3_%tsp** : proportion of time spent in N3 relative to TSP (N3/TSP*100) - in percentage

    **REM_%tsp** : proportion of time spent in REM relative to TSP (REM/TSP*100) - in percentage

    **SSI** : stage switching index (number of change from one stage to another) - number per hour (TSP)

    **SFI** : sleep fragmentation index (number of change from one stage to a lighter stage) - number per hour (TSP)

    **SL_toN2_min** : sleep latency to reach first epoch of N2 - in minutes

    **SL_toN3_min** : sleep latency to reach first epoch of N3 - in minutes

    **SL_toREM_min** : sleep latency to reach first epoch of REM - in minutes

    **SL_toNREM_5m_min** : sleep latency to reach 5 minutes of consolidated NREM (N2+N3) - in minutes

    **SL_toNREM_10m_min** : sleep latency to reach 10 minutes of consolidated NREM (N2+N3) - in minutes

    **SL_toN3_5m_min** : sleep latency to reach 5 minutes of consolidated N3 - in minutes

    **SL_toN3_10m_min** : sleep latency to reach 10 minutes of consolidated N3 - in minutes
