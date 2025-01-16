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
and as such the data will need to be organised according to `bids-specification <https://bids-specification.readthedocs.io/en/stable/common-principles.html#source-vs-raw-vs-derived-data>` 
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
      │  │  │     ├─ sub-01_ses-01_task-sleep_acq-PSG_events.tsv   *optional for seapipe*
      │  │  │     └─ sub-01_ses-01_task-sleep_acq-PSG_channels.tsv *optional for seapipe*
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
   1. You can look at example bids-standard EEG datasets `here <https://openneuro.org/search/modality/eeg>`_.
   2. It is recommended that you use a `bids-validator <https://bids-standard.github.io/bids-validator/>`_ on the ``/sourcedata`` folder *prior to using seapipe*.

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
For seapipe to run properly, the data needs to be organised in the **Brain Imaging Data Structure (BIDS)**.
The compatibility of the dataset with BIDS can be validated `online <https://bids-standard.github.io/bids-validator/>`_.

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

*[To be completed]*

.. The ``kind`` parameter should be either ``"meat"``, ``"fish"``,
.. or ``"veggies"``. Otherwise, :py:func:`lumache.get_random_ingredients`
.. will raise an exception.

.. .. autoexception:: lumache.InvalidKindError

.. For example:

.. >>> import lumache
.. >>> lumache.get_random_ingredients()
.. ['shells', 'gorgonzola', 'parsley']

