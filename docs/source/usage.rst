Usage
=====

.. _installation:

Installation
------------

To use seapipe, first install it using pip:

.. code-block:: console

   (.venv) $ pip install seapipe

.. _creating_a_pipeline:
Creating a pipeline
----------------
To begin, open python and load seapipe

.. code-block:: console

   (.venv) $ python
>>> from seapipe import pipeline

Then you can initiate the pipeline by specifying the path to your dataset.

>>> project_name = pipeline('/home/username/project/') 

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
      
::
   2024-12-02 18:42:41 - Tracking - Slow oscillation detection has NOT been run. 

                           ses      staging      spindle slow_osc        fooof
   sub-IN001  [ses-V1, ses-V2]  [ses-V1, -]  [ses-V1, -]      [-]  [ses-V1, -]
   sub-IN003          [ses-V1]     [ses-V1]     [ses-V1]      [-]          [-]
   sub-IN076               [-]          [-]          [-]      [-]          [-]

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

