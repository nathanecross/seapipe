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

>>> pipeline = pipeline('/home/username/project/') 

.. _checking_your_dataset:
Checking your dataset
----------------

To retrieve a list of all the files inside the root directory, along with the
directories 1 and 2 levels preceding the files,
you can use the ``pipeline.list_dataset()`` function:

>>> pipeline.list_dataset()

.. code-block:: console
   Directory: projects/bids
   Files = ['dataset_description.json', 'participants.tsv']
   ----------
   Directory: ses-m/eeg
   Files = ['sub-001_ses-01_eeg.edf']
   ----------
   Directory: ses-c/eeg
   Files = ['sub-001_ses-02_eeg.edf']
   ----------
   Directory: ses-m/eeg
   Files = ['sub-002_ses-01_eeg.edf']
   ----------
   etc.
s
The ``kind`` parameter should be either ``"meat"``, ``"fish"``,
or ``"veggies"``. Otherwise, :py:func:`lumache.get_random_ingredients`
will raise an exception.

.. autoexception:: lumache.InvalidKindError

For example:

>>> import lumache
>>> lumache.get_random_ingredients()
['shells', 'gorgonzola', 'parsley']

