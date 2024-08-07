Macro-architecture
=====

.. _overview:

Overview
------------
You can extract the different markers of macro-architecture (see definitions in section :ref:`Output`) whole night or per cycle

    **adap_bands**
        * Precise whether you performed a detection using adapted frequency range (i.e., tailored to the individual's peak frequency of a specific frequency band) 
        for *detect_slow_oscillations*

        * *Acceptable options:*

            * Default is ``False`` which infer a traditional detection of events based on fixed frequency range which will point to frequency 
            input - *recommended for SO detection*

            * If you put ``Auto``, it will point to the ``root_dir/OUT/fooof`` folder, which comprises the band of interest's peak frequency for each subject, session and channel

            * If you put ``Manual``, it will point to the *chanset_peaks* column in *tracking* file

    **adap_bw**
        * Bandwidth for frequency range (Hz) - only relevant if you passed ``adap_bands=Auto`` or ``adap_bands=Manual``

        * *Acceptable options:*

            * Default is ``4`` which infers 2hz before and after the individual's peak frequency (fooof) for each channel

            * You can put another even number to expand the bandwidth for frequency range


*You will need to run two functions:*

1. Extract macro sleep characteristics for each subject.
    It will extract a .csv file including macro-architecture variables wholenight and per cycle for each subject and each session in ``root_dir/OUT/staging/``

.. code-block:: python

   project_name.export_macro_stats()


2. Create datasets combining all the subjects
    It will combine all .csv into a single dataset per session (one row per subject)

.. code-block:: python

   project_name.macro_dataset()
 

.. _extraction_macro:
Extract macro-architecture
----------------
*Command line argument:*

.. code-block:: python

   project_name.export_macro_stats(xml_dir = None, 
                                   out_dir = None, 
                                   subs = 'all', 
                                   sessions = 'all', 
                                   times = None, 
                                   rater = None, 
                                   outfile = True)


*Positional arguments:*

    **xml_dir**
        * Path to folder with the .xml file containing sleep stages and arousal events. 

        * Default is ``None`` which will point to ``root_dir/OUT/staging``

    **out_dir**
        - Output path for the outcomes of charactertistics extraction per subject.

        - Default is ``None`` which will point to ``root_dir/OUT/staging``

    **subs**
        * Subject to analyze

        * Default is ``'all'`` which will point to all the *sub* folders in ``root_dir/DATA``

            * If put ``None``, it will point to the *sub* column in *tracking* file
            * If put string of sub ID (e.g., *['sub-01', 'sub-02']*), it will only detect those sub folders

    **sessions**
        * Sessions/Visits to analyse per subject

        * Default is ``'all'`` which will point to all the *ses* folders within the sub folder in ``root_dir/DATA``

            * If put ``None``, it will point to the *ses* column in *tracking* file

            * If put string of ses visit (e.g., *['ses-V1']*), it will only detect the selected session(s) within each subject

    **times**
        * Light off and light on in seconds from beginning of recording

        * Default is ``None`` which will point to the *loff* and *lon* columns in *tracking* file

    **rater**
        * Name of the rater to analyze

        * Default is ``None`` which will discard the name of the rater and expect only one rater per .xml (!! make sure you don't have multiple raters!!)
    
            * If put string of rater's name (e.g., *[Rater1]*), it will only extract sleep architecture from this rater per .xml (and create an empty extraction file if the rater is absent)

    **outfile**
        * Extraction of output file

        * Default is ``True`` which will create a .csv file per subject and per session in ``root_dir/OUT/staging/``
            
            * If put ``False``, it won't extract .csv file of macro-sleep characteristics which will impact creation of datasets


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
