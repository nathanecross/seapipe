Signal preprocessing
=====

.. _overview:

Overview
------------
Signal quality is a critical determinant of reliable electrophysiological analysis, particularly in sleep EEG where artefacts, 
channel dropouts, and physiological contamination can vary substantially across recordings. Automated and standardised quality 
control (QC) procedures are therefore essential to ensure that downstream analyses are not biased by poor signal integrity.

Seapipe provides this functionality via S.Q.U.I.D. (Signal Quality, Uncertainty, Inspection and Denoising) - an automated framework 
to quantify channel-level signal quality using a range of complementary metrics. These metrics capture temporal stability, amplitude 
plausibility, spectral composition, signal distribution, and physiological contamination, enabling the identification and flagging of 
low-quality channels. QC outputs can be generated for individual recordings and aggregated across datasets, facilitating both data-driven 
channel selection and transparent reporting of signal quality.

These cases are different from :ref:`Artefacts` in that they contaminate a signal for the majority of the recording.
Some examples of poor signal quality include:

1. Hardware noise

    .. image:: images/noise.png
        :width: 700

2. ECG contamination

    .. image:: images/ecg_contamination.png
        :width: 700


| S.Q.U.I.D provides a summary of signal quality using the following metrics:

    * time_not_flatline: Proportion of time the signal is not flat (based on a moving average).

    * time_above_10:  Proportion of time the signal amplitude is above 10 µV. Flags channels where the signal is of abnormally low amplitude.

    * time_below_200:  Proportion of time the signal amplitude is above 200 µV. Flags channels where the signal is of abnormally high amplitude.
    
    * `gini_coefficient <https://doi.org/10.1155/2022/7731309>`_: Gini coefficient of signal amplitude distribution. The Gini coefficient, 
        traditionally used to measure economic inequality, is applied in signal processing to measure sparsity, impulsiveness, or concentration. 
        A Gini coefficient closer to 1 indicates that the signal energy is concentrated in a few components (high sparsity, low noise), while a 
        value closer to 0 indicates a uniform distribution (low sparsity, high noise/spread).
    
    * inverse_power_ratio:  Inverse of the ratio of high-frequency to low-frequency power. 

    * ecg_artefact:  Measure of ECG contamination based on average correlation between EEG signal and ECG peak . 

    * ecg_artefact_perc:  Proportion of EEG signal time with high (r>0.5) correlation with ECG (i.e. contamination).


.. admonition:: When to use S.Q.U.I.D
    
    
    S.Q.U.I.D. is designed for use at the preprocessing or quality control stage of electrophysiological analyses, particularly 
    when working with large datasets, multi-channel recordings, or heterogeneous data sources. It is especially useful when 
    channel quality is uncertain, variable across recordings, or when automated channel selection is required prior to downstream 
    analyses (e.g., event detection, spectral analysis, or connectivity measures). By providing standardised, quantitative QC metrics, 
    S.Q.U.I.D. facilitates reproducible data cleaning, objective channel inclusion criteria, and transparent reporting of data quality 
    across subjects and sessions.
    

.. _Functions:
Functions to implement S.Q.U.I.D
----------------
| **Running S.Q.U.I.D on individual recordings:**

1) Calculate QC metrics:  

.. code-block:: python

   project.QC_channels(()
|
    This will run the function on every ``/sub-XXX/ses-XXX`` in ``<rec_dir>`` and write QC metrics in to ``<root_dir>/derivatives/QC/``. 
|
2) Export QC dataset: 

.. code-block:: python

   project_name.QC_summary()
|   
    This will gather all QC reports, including every channel, for a given modality into a single ``.csv`` DataFrame file in ``<root_dir>/derivatives/datasets/QC/`` 
|


.. _qc_channels:
Calculate QC metrics
----------------
*Command line argument:*

.. code-block:: python

    project.QC_channels(subs = 'all', 
                        sessions = 'all', 
                        filetype = '.edf',
                        filt = None, 
                        chantype = ['eeg', 'eog', 'emg', 'ecg'],
                        outfile=True)


*Positional arguments:*

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

    **filt** *(str or NoneType)*
        * Explicit list of channels to process for EEG.

        * *Acceptable options:*

            * Default is ``None``  Standard 10-20 EEG channel names will be used. 
            
            * If entering a ``str``, then this will work to filter only EEG channel names containing this string in the name. 

    **chantype** *(list)*
        * Channel type(s) to process. Accepted values: .

        * *Acceptable options:*

            * Default is ``['eeg', 'eog', 'emg', 'ecg']`` which will process each of the channel types.

            * Entering a subset of this list is accepted, e.g. ``['eeg']``

    **outfile** *(bool or str)*
    
        * Wheter to log the process. 

        * *Acceptable options:*

            * ``True``, writes log to auto-generated timestamped file in `self.log_dir`.

            * ``str``, writes log to the specified file path.

            * ``False``, logs only to console.

    


.. _export_QC:
Export QC metrics to dataset
----------------
*Command line argument:*

.. code-block:: python

    project.QC_summary(qc_dir = None, 
                       chantype = ['eeg', 'eog', 'emg', 'ecg'])


*Positional arguments:*

    **qc_dir** *(str)*
        * Path to the directory with sub-directories ``/sub-XXX/ses-XXX`` containing the saved QC reports. 

            * Default is ``None`` which will point to ``<root_dir>/derivatives/QC/``

    **chantype** *(list)*
        * Channel type(s) to process. Accepted values: .

        * *Acceptable options:*

            * Default is ``['eeg', 'eog', 'emg', 'ecg']`` which will process each of the channel types.

            * Entering a subset of this list is accepted, e.g. ``['eeg']``












