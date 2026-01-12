Phase Amplitude Coupling
=====

.. _overview:

Overview
------------

Phase amplitude coupling (PAC) is a phenomenon where the amplitude of a high-frequency oscillation statistically changes in relation to the phase of a low-frequency oscillation. 
It is a well-studied interaction between brain oscillations in different frequency bands (e.g. slow oscillations and the spindle *sigma* frequency band).
There is no convention yet of how to calculate phase-amplitude coupling.

| Seapipe (via `Tensorpac <https://github.com/EtienneCmb/tensorpac>`_) provides 6 published methods to automatically calculate PAC:

    * `'Mean Vector Length (MVL) [Canolty et al. 2006]' <https://www.science.org/doi/10.1126/science.1128115>`_:
    
       **Method in brief**: The phase-amplitude coupling measure MVL quantifies the coupling between low-frequency phase and high-frequency amplitude by averaging the vectors of the analytic signal, 
       where the length of the mean vector represents the strength of coupling and its direction indicates the phase at which amplitude is strongest. In the absence of coupling, the vectors cancel out, 
       resulting in a short mean vector with no meaningful phase direction.

    * `'Modulation Index (MI) [Tort 2010]' <https://journals.physiology.org/doi/full/10.1152/jn.00106.2010>`_:
    
       **Method in brief**: This metric is based on the Kullback-Leibler divergence, a model-free, information theoretic measure of the distance between two distributions. By calculating Shannon entropy, 
       MI measures the level of irregularity in the distribution of the mean amplitudes of a high frequency signal binned by the concomitant phase of a low frequency signal. The MI ranges from 0 to 1, 
       where 0 indicates a uniform distribution of amplitudes across all phase bins, and 1 signifies a total amplitude concentration inside a single-phase bin.

    * `'Heights Ratio (HR) [Lakatos et al. 2005]' <https://journals.physiology.org/doi/full/10.1152/jn.00263.2005>`_:
    
       **Method in brief**: Instantaneous power and phase are extracted by wavelet decomposition (Morlet wavelet). Oscillatory amplitude is aligned to oscillatory phase by sorting phase values 
       from −π to π radians and applying the corresponding permutation vector to rearrange the amplitude time series. Pooled amplitude and frequency values are evaluated statistically by ANOVA. 

    * `'ndPAC [Ozkurt 2012]', <https://ieeexplore.ieee.org/document/6184293>`_:
    
       **Method in brief**: The normalized direct PAC estimate (ndPAC) uses a statistical threshold to ensure reliable coupling estimates while eliminating distortions from DC components in the data.
       The amplitude vector is normalized by removing its mean and scaling it to have unit variance. For each time point 𝑛, the normalized amplitude is multiplied by the complex phase, creating a set of complex 
       vectors. These vectors are summed over all time points to produce a single complex value. The magnitude (length) of this summed vector measures the consistency of the relationship between amplitude and 
       phase across time.
 
    * `'Phase-Locking Value (PLV) [Penny et al. 2008]' <https://www.sciencedirect.com/science/article/pii/S0165027008003816>`_:
    
       **Method in brief**: To calculate the Phase Locking Value (PLV), phase angles are extracted from both a low-frequency filtered signal and the Hilbert-transformed high-frequency amplitude signal, 
       and the phase angle differences are computed for each time point. The coupling strength is determined by the length of the mean vector obtained from averaging these phase difference vectors in the polar plane, 
       where a constant phase lag results in a longer mean vector and stronger coupling.

    * `Gaussian Copula PAC (GCPAC) [Ince 2017 et al. 2013] <https://onlinelibrary.wiley.com/doi/10.1002/hbm.23471>`_:
    
       **Method in brief**: This is an approach to estimate MI with continuous variables (rather than binning data) by leveraging mutual information and Gaussian copula (a statistical structure that expresses a 
       multivariate distribution as the combination of univariate marginal distributions). Amplitude is calculated as the absolute value of the complex number and optionally squared for power. Phase is isolated by 
       normalizing the complex number by its amplitude, resulting in a 2D unit circle representation where Gaussian Copula Mutual Information (GCMI) is applied, ensuring a valid estimate of mutual information (MI) 
       carried by phase without adding extraneous information.


There are also various implemented methods in order to generate the distribution of surrogates and then to correct the PAC for spurious coupling:

    * No surrogates
    * Swap phase / amplitude across trials `[Tort 2010] <https://journals.physiology.org/doi/full/10.1152/jn.00106.2010>`_
    * Swap amplitude time blocks `[Bahramisharif 2013] <https://www.jneurosci.org/content/33/48/18849>`_
    * Introducing a time lag on phase series `[Canolty et al. 2006] <https://www.science.org/doi/10.1126/science.1128115>`_


Finally, the true PAC is corrected by the calculated surrogate distribution. Tensorpac includes four types of normalization that should give similar results:

    * Substract the mean of surrogates (1)
    * Divide by the mean of surrogates (2)
    * Substract then divide by the mean of surrogates (3)
    * Substract the mean then divide by the deviation of surrogates (z-score, 4)


To read more about thse approaches, `see here <https://etiennecmb.github.io/tensorpac/auto_examples/index.html#tutorials>`_



.. _detection_pac:
Run Phase Amplitude Coupling
----------------

|
    This will copy the :ref:`Annotations file` from every ``/sub-XXX/ses-XXX`` in ``<xml_dir>`` to ``<root_dir>/OUT/pac/`` and calculate PAC. Output parameters will be stored in ``_pac_parameters.csv``
|

*Command line argument:*

.. code-block:: python

    seapipe.pac(xml_dir = None, out_dir = None, 
                subs = 'all', sessions = 'all', 
                filetype = '.edf',
                chan = None, ref_chan = None, rater = None, grp_name = 'eeg', 
                stage = ['NREM2','NREM3'], concat_stage = True, 
                cycle_idx = None, concat_cycle = True,  
                method = 'MI', surrogate = 'Time lag', correction = 'Z-score',
                evt_name = None, min_dur = 1, nbins = 18, invert = None,
                adap_bands_phase = 'Fixed', frequency_phase = (0.5, 1.25), 
                adap_bands_amplitude = 'Fixed', frequency_amplitude = (11, 16),
                adap_bw = 4,
                frequency_opts = None, 
                filter_opts = None, 
                epoch_opts = None, 
                event_opts = None, 
                reject_artf = ['Artefact', 'Arou', 'Arousal'], 
                progress = True, 
                outfile = True)


*Positional arguments:*

    **xml_dir** *(str)*
        * Path to the directory with sub-directories ``/sub-XXX/ses-XXX`` containing the input :ref:`Annotations files<Annotations file>`. 

        * Default is ``None`` which will point to ``<root_dir>/OUT/staging/`` (Annotations files with sleep stage markings and arousal/artefact events).

    **out_dir** *(str)*
        * Output path for the .xml file containing the new detected event (events will be named like the method used; e.g., ``Ray2015``)

        * Default is ``None`` which will point to ``<root_dir>/OUT/pac/``

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

    **concat_stage** *(logical)*
        * Concatenation options for sleep stages

        * *Acceptable options:*

            * Default is ``True`` which means that detection will be performed per stage

            * Entering ``False`` which means that all stages will be concatenated (i.e., merged) before detection **It is recommended that you leave the default option**

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

    **method** *(str)*
        * Method of calculating phase amplitude coupling

        * *Acceptable options:*

            * Default is ``'MI'``  which is the `'Modulation Index (MI) [Tort 2010]' <https://journals.physiology.org/doi/full/10.1152/jn.00106.2010>` method  
            
            * Other available methods include: ``'MVL', 'HR', 'ndPAC', 'PLV', 'GCPAC'`` (see :ref:`Overview<overview>`)

    **surrogate** *(str)*
        * Method of calculating surrogate (artificial) distribution for correcting spurious coupling.

        * *Acceptable options:*

            * Default is ``'Time lag'`` which involves introducing a time lag on phase series `[Canolty et al. 2006] <https://www.science.org/doi/10.1126/science.1128115>`_

            * Other available methods include: ``'No surrogates', 'Swap phase', 'Swap amplitude'`` (see :ref:`Overview<overview>`)
        

    **correction** *(str)*
        * Method of correcting correcting spurious coupling using the surrogates.

        * *Acceptable options:*

            * Default is ``'Z-score'`` which involves subtracting the mean then dividing by the standard deviation of surrogates.

            * Other available methods include: ``'No normalization', 'Subtract', 'Divide', 'Subtract then divide'`` (see :ref:`Overview<overview>`)

    **evt_name** *(NoneType or str)*
        * Event name to run PAC across (e.g. slow oscillations). Events will be isolated before phase and amplitude filtering is applied.

        * * *Acceptable options:*

            * Default is ``None`` which will perform PAC across continuous signal (e.g. all NREM or all REM or all NREM3) depending on how the cycles and stages are concatenated.

            * If entering an event name (e.g. ``'SO'``) that event will need to be already detected in the :ref:`Annotations file` and named exactly as entered here.

    **min_dur** *(float)*
        * Minimum duration (seconds) of events that will be detected. Any events with durations that are outside these limits will be discarded

        * *Acceptable options:*

            * Default is ``1`` (in seconds)

            * Entering a float (e.g., ``1.4``)  will limit the detection to events with a duration above this threshold

    **nbins** *(int)*
        * Number of phase bins to discretize the signal for calculation of preferred phase and coupling strength.

        * *Acceptable options:*
            * Default is ``18`` which will provde phase bins of 20˚ (ie. 360˚/18)

            * Any integer is allowed, but it is recommended that it be a factor of 360.

    **invert** *(NoneType or logical)*
        * Option to invert polarity

        * *Acceptable options:*

            * Default is ``None`` which will point to the *chanset_invert* columns in the :ref:`tracking file<Tracking File>`. However, if the *tracking* file does not specify *chanset_invert* 
            columns, the detection will default to ``False``

            * Entering ``False`` will keep the polarity of the recording as it is

            * Entering ``True`` will reverse (flip) the polarity of the recording 

    **adap_bands_phase** *(str)*
        * Options to set an adapted frequency band tailored to each individual for the phase portion of the PAC

        * *Acceptable options:*

            * Default is ``'Fixed'`` which will point to the frequency range set up in **frequency_phase**

            * Entering ``'Auto'`` will perform :ref:`FOOOF analyses<FOOOF analyses>` which will detect the peak in sigma characterized in terms of their specific
            center frequency, power and bandwidth within the frequency range set up in **frequency_phase** and controlling for the aperiodic component. By default, if left 
            ``frequency = None``, the range set-up for fooof peak detection is 0.5-1.25Hz. **THIS IS NOT RECOMMENDED FOR LOW FREQUENCY RANGES**. It will add *_adap_phase* at 
            the end of ``out_dir`` (e.g., *pac_adap_phase*).

            * Entering ``Manual`` will point to the *chanset_peaks* columns in the :ref:`tracking file<Tracking File>`. It will add *_adap_phase* at the end of ``out_dir`` (e.g., *pac_adap_phase*).

    **frequency_phase** *(tuple)*
        * Frequency range of interest 

        * *Acceptable options:*

            * Default is ``None`` which will depend to the options selected for **adap_bands**. If ``adap_band = 'Fixed'``, frequency will be (11,16) while ``adap_band = 'Auto'``
            will be (9,16) for the peak frequency detection

            * Enter a `tuple <https://docs.python.org/3/tutorial/datastructures.html#tuples-and-sequences>`_ containing the frequency range of interest that 
            will be used if selecting ``adap_bands = 'Fixed'`` or ``adap_bands = 'Auto'`
 
    **adap_bands_amplitude** *(str)*
        * Options to set an adapted frequency band tailored to each individual for the phase portion of the PAC

        * *Acceptable options:*

            * Default is ``'Fixed'`` which will point to the frequency range set up in **frequency_amplitude**

            * Entering ``'Auto'`` will perform :ref:`FOOOF analyses<FOOOF analyses>` which will detect the peak in sigma characterized in terms of their specific
            center frequency, power and bandwidth within the frequency range set up in **frequency_amplitude** and controlling for the aperiodic component. By default, if left 
            ``frequency = None``, the range set-up for fooof peak detection is 0.5-1.25Hz. It will add *_adap_amplitude* at 
            the end of ``out_dir`` (e.g., *pac_adap_amplitude*).

            * Entering ``Manual`` will point to the *chanset_peaks* columns in the :ref:`tracking file<Tracking File>`. It will add *_adap_amplitude* at the end of ``out_dir`` 
            (e.g., *pac_adap_amplitude*).

    **frequency_amplitude** *(tuple)*
        * Frequency range of interest 

        * *Acceptable options:*

            * Default is ``None`` which will depend to the options selected for **adap_bands**. If ``adap_band = 'Fixed'``, frequency will be (11,16) while ``adap_band = 'Auto'``
            will be (9,16) for the peak frequency detection

            * Enter a `tuple <https://docs.python.org/3/tutorial/datastructures.html#tuples-and-sequences>`_ containing the frequency range of interest that 
            will be used if selecting ``adap_bands = 'Fixed'`` or ``adap_bands = 'Auto'`

    **adap_bw** *(str or float)*
        * Size of the frequency range around sigma peak frequency when entering ``Auto``or ``Manual`` to **adap_bands**

        * *Acceptable options:*

            * Default is ``4``meaning 2Hz on both side of the sigma peak frequency

            * Any `float <https://docs.python.org/3/tutorial/floatingpoint.html>`_ is allowed

    **reject_artf** *(list)*
        * Options to discard detection within specific events such as Artefact events

        * *Acceptable options:*

            * Default is ``['Artefact', 'Arou', 'Arousal']``which will discard detection during events with these specific names

            * Entering a list of events will discard detection within those events

    **frequency_opts** *(NoneType or dict)*
        * Options for parameters for power spectral analyses 

        * *Acceptable options:* 
            * For formatting the dictionary, see :ref:`Power spectrum<Power_spectrum>`

            * Entering ``None`` will use default parameters for power spectral analyses.

    **filter_opts** *(NoneType or dict)*
        * Options for parameters for filtering 

        * *Acceptable options:* 
            * For formatting the dictionary, see :ref:`Power spectrum<Power_spectrum>`
            
            * Entering ``None`` will use default parameters for power spectral analyses.

    **epoch_opts** *(NoneType or dict)*
        * Options for parameters for epoch analyses 

        * *Acceptable options:* 
            * For formatting the dictionary, see :ref:`Power spectrum<Power_spectrum>`
            
            * Entering ``None`` will use default parameters for power spectral analyses.

    **event_opts** *(NoneType or dict)*
        * Options for parameters for event analyses 

        * *Acceptable options:* 
            * For formatting the dictionary, see :ref:`Power spectrum<Power_spectrum>`
            
            * Entering ``None`` will use default parameters for power spectral analyses.

    **progress** *(logical)*
        * Show the progress bar for each ``sub`` and ``ses``. 

        * Default is ``True`` - set to ``False`` if running on HPC clusters.

    **outfile** *(str or logical)*
        * Logging of PAC

        * *Acceptable options:*

            * Default is ``True`` which will create a logfile *detect_pac_{method}_{datetime}_log.txt* in ``<root_dir>/OUT/audit/logs/``

            * Entering a string ``<custom_outfile_name.txt>`` will save the logfile under that custom name
            
            * Entering ``False`` won't save a logfile











