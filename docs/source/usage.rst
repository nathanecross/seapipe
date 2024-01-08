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

   Participants   BIDS?	   #sessions #recordings   #annotations
   sub-002        TRUE	      2	         2	            2
   sub-004        TRUE	      3	         3	            3
   sub-006	  FALSE	      2	         2	            1
   sub-007	  TRUE	      2	         2	            2
   sub-008	  FALSE	      2	         1	            1
   sub-009	  TRUE	      2	         2	            2
   sub-011	  FALSE	      0	         0	            0
   sub-013	  TRUE	      2	         2	            2
   sub-014	  FALSE       0          2                  2
   sub-015	  TRUE	      2	         2	            2
   sub-016	  FALSE	      2	         2	            0


This will be automatically saved to a file *dataset_audit.csv*

To retrieve a list of all the files inside the root directory, along with the
directories 1 and 2 levels preceding the files,
you can use the ``pipeline.list_dataset()`` function:

>>> pipeline.list_dataset()

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

.. The ``kind`` parameter should be either ``"meat"``, ``"fish"``,
.. or ``"veggies"``. Otherwise, :py:func:`lumache.get_random_ingredients`
.. will raise an exception.

.. .. autoexception:: lumache.InvalidKindError

.. For example:

.. >>> import lumache
.. >>> lumache.get_random_ingredients()
.. ['shells', 'gorgonzola', 'parsley']

