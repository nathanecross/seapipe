Usage
=====

.. _installation:

Installation
------------

To use seapipe, first install it using pip:

.. code-block:: console

   (.venv) $ pip install seapipe

.. _data_preparation_and_setup:
Data Preparation and Setup
----------------

Seapipe is a `bids-standard <https://bids-specification.readthedocs.io/en/stable/>`_ data processing pipeline, 
and as such the data will need to be organised according to `bids-specification <https://bids-specification.readthedocs.io/en/stable/common-principles.html#source-vs-raw-vs-derived-data>`_ 
for seapipe to run properly. 
By running seapipe, raw data will be transformed and partial as well as final results will be saved. Therefore, the original, raw data must be separated from the outputs (derivatives).
This is done by placing the raw dataset inside a directory labelled ``sourcedata`` inside the root directory.

For example, a eeg datafile should be in the structure ``~/rootdir/sourcedata/sub-01/ses-01/eeg/sub-01_ses-01_task-sleep_acq-PSG_eeg.edf``

An example of the datastructure would look like this:
::
   └─ my_project-1/
      ├─ sourcedata/
      │  ├─ sub-01/
      │  │  ├─ ses-01
      │  │  │  └─ eeg
      │  │  │     ├─ sub-01_ses-01_task-sleep_acq-PSG_eeg.edf
      │  │  │     ├─ sub-01_ses-01_task-sleep_acq-PSG_eeg.json
      │  │  │     ├─ sub-01_ses-01_task-sleep_acq-PSG_events.tsv     *optional for seapipe*
      │  │  │     └─ sub-01_ses-01_task-sleep_acq-PSG_channels.tsv   *optional for seapipe*
      │  │  ├─ ses-02/
      │  │  └─ ...
      │  ├─ sub-02/
      │  ├─ ... 
      │  ├─ dataset_description.json 
      │  └─ participants.tsv
      └─ derivatives/
         ├─ seapipe/
         └─ ...


.. admonition:: NOTE - BIDS🧠
   1 - You can look at example bids-standard EEG datasets `here <https://openneuro.org/search/modality/eeg>`_.
   2 - It is recommended that you use a `bids-validator <https://bids-standard.github.io/bids-validator/>`_ on the ``/sourcedata`` folder *prior to using seapipe*.

.. _creating_a_pipeline:
Creating a pipeline
----------------
To begin, open python and load seapipe

.. code-block:: console

   (.venv) $ python
>>> from seapipe import pipeline

Then you can initiate the pipeline by specifying the path to your dataset.

>>> project_name = pipeline('~/my_project-1/') 

.. _checking_your_dataset:
Checking your dataset
----------------

Before running any analyses, it is important to check your data.
For seapipe to run properly, the data needs to be organised in the **Brain Imaging Data Structure (BIDS)** (see :ref:`Data Preparation and Setup`).

However, seapipe also works almost symbiotically with the `Wonambi <https://wonambi-python.github.io/>`_ package.
Therefore, any annotations (sleep scoring, artefact markings etc.) need to be inside a wonambi annotations (.xml) file. 
For more information, see :ref:`Annotations file`.

To receive an overview of your dataset, including whether the each participant's directory is BIDS compatible, as well as 
how many sessions, recording (e.g. edfs) and annotation files they contain, you can call the ``pipeline.audit`` property 
of every dataset:
 
>>> pipeline.audit()
 
::

                      Summary:
                      2 files, 2.80 GB
                      Subjects: 2
                      Sessions: 2

   2024-12-02 18:35:54 - Audit - The dataset appears compatible for SEAPIPE analysis.

             BIDS?  #sessions  #recordings    
   sub-001  False          2            1  !!
   sub-003  False          1            1    


This will be automatically saved to a file *dataset_audit.csv*

To retrieve a list of all the files inside the root directory, along with the
directories 1 and 2 levels preceding the files,
you can use the ``pipeline.list_dataset()`` function:

>>> project_name.list_dataset()

:: 

   Directory: project/bids
   Files = ['dataset_description.json', 'participants.tsv']
   ----------
   Directory: ses-01/eeg
   Files = ['sub-001_ses-01_eeg.edf']
   ----------
   Directory: ses-02/eeg
   Files = ['sub-001_ses-02_eeg.edf']
   ----------
   Directory: ses-01/eeg
   Files = ['sub-002_ses-01_eeg.edf']
   ----------
   etc.


To retrieve a table of all the analyses that have been run (and are located in ``<root_dir>/OUT/``), run the following command:

.. code-block:: python

    project.track(subs = 'tracking.tsv',
                  step=['staging','spindles', 'so', 'fooof'], 
                  chan = ['Fz (eeg)'],
                  outfile='progress.csv')

This will output a table of each stage provided for the subs, sessions and channels specified:
::
   2024-12-02 18:42:41 - Tracking - Slow oscillation detection has NOT been run. 

                           ses      staging      spindle slow_osc        fooof
   sub-001  [ses-V1, ses-V2]  [ses-V1, -]  [ses-V1, -]      [-]  [ses-V1, -]
   sub-003          [ses-V1]     [ses-V1]     [ses-V1]      [-]          [-]
   sub-004               [-]          [-]          [-]      [-]          [-]

|
.. _tracking file:
Tracking File
----------------

Uniformity in EEG electrode placement is crucial for ensuring consistent signal capture, minimizing artifacts, 
and improving comparability across recordings. EEG measures scalp electrical activity, meaning even slight variations 
in electrode positioning can alter recorded signals, affecting amplitude and frequency analyses, and source localization 
accuracy. Unlike MRI, which provides high-resolution brain images and allows for spatial normalization to a common 
template, EEG lacks a direct post hoc standardization method, making uniform electrode placement essential.

This has led to systems of EEG application, most notably the `10-20 system. <https://en.wikipedia.org/wiki/10%E2%80%9320_system_(EEG)>`_

Therefore, when working with EEG data, each timeseries is affiliated to a source electrode. And because EEG is a measure 
of `electrical potentials <https://doi.org/10.1016/j.cub.2018.11.052>`_ there is the need for reference channels.

However, despite uniformity in spatial placement of recording electrode sites, not all recording software use the same 
EEG configurations (e.g. channel names, online references, sampling_frequencies etc). This can cause headaches when
trying to conduct pipeline analyses across datasets with inconsistences in these certain parameters.

One way that **seapipe** gets around this is with the use of a tracking file. This file can be in .tsv or .xlsx format.
However it *must* be named: **tracking.csv**

It's structure should look like this:
::
   sub        ses         loff      lon       format     chanset1            chanset1_rename    refset1
   sub-01     ses-1       330       31500     .edf       F3, C3              F3, C3             M1, M2
   sub-01     ses-2       4320      32390     .edf      	F3, C3              F3, C3             M1, M2
   sub-02     ses-1       1900      29945     .edf       F3 (A2), C3 (A2)    F3, C3             A1, A2
   sub-02     ses-2       670       31010     .edf       F3 (A2), C3 (A2)    F3, C3             A1, A2
   ...
|
As you can see with this dataset, there are some inconsistences in the channel naming: 

sub-01 has channels named 'F3' and 'C3' <-> sub-01 has channels named 'F3 (A2)' and 'C3 (A2)'

sub-01 has references named 'M1' and 'M2' <-> sub-01 has channels named 'A2' and 'A2'

All subjects and sessions have different lights out (loff) and lights on (lon) times, corresponding to the `time in bed. <https://www.sleepfoundation.org/how-sleep-works/sleep-dictionary#:~:text=Time%20in%20bed%3A%20The%20total,studies%20to%20calculate%20sleep%20efficiency.>`_

If you create this tracking file, then you can read parameters such as channel names by setting this:
.. code-block:: python

   chan = None
|
 in the calls to functions.

** Coming soon ** The function to read from a channels.tsv file in a BIDS dataset

.. The ``kind`` parameter should be either ``"meat"``, ``"fish"``,
.. or ``"veggies"``. Otherwise, :py:func:`lumache.get_random_ingredients`
.. will raise an exception.

.. .. autoexception:: lumache.InvalidKindError

.. For example:

.. >>> import lumache
.. >>> lumache.get_random_ingredients()
.. ['shells', 'gorgonzola', 'parsley']

