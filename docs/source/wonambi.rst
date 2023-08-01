Wonambi
=====

.. _Overview:

Overview
------------

`Wonambi <https://wonambi-python.github.io/>`_ is an open-source package for the analysis of EEG, ECoG 
and other electrophysiology modalities. It allows for the visualization of data and sleep stage scoring 
in a GUI, as well as providing automatic detectors for sleep spindles and slow waves. 

It is a dependency of **seapipe** and therefore it is included when you install seapipe.

.. _reading_data:
Reading data
----------------
To begin, open python and load wonambi:

.. code-block:: console

   (.venv) $ python
>>> from wonambi import dataset


Then you can load a recording by specifying the path to your data.

>>> data = Dataset('/home/username/project/participant/session/eeg/sub-001_ses-01_eeg.edf') 

.. _Annotations file:
Annotations file
----------------

When working with EEG, it is often necessary to tag, or annotate, certain events or periods of
the recording that will be analysed at a later stage. In the context of **sleep analyses,** the 
most common example of this would be `sleep stage scoring <https://aasm.org/clinical-resources/scoring-manual/>`_,
but other events can also be annotated (e.g. spindles, slow oscillations, arousals etc.)

Wonambi utilises an Annotations file, which is in the XML format

.. admonition:: BIDS 🧠

   When working with BIDS, `it has been proposed <https://www.nature.com/articles/s41597-019-0104-8>`_ that an *events.tsv* file 
   should exist that specifies all events which have been recorded during the session (e.g. which reference presented stimuli 
   with the stimuli directory of the dataset). To learn more about converting between Wonambi's annotations format and BIDS 
   compatible *events.tsv*, see :ref:`Converting between Formats`

First load in the function Annotations

>>> from wonambi.attr import Annotations

Then,


.. _Converting between Formats:


.. Before running any analyses, it is important to check your data.
.. For seapipe to run properly, the data needs to be organised in the **Brain Imaging Data Structure (BIDS)**.
.. The compatibility of the dataset with BIDS can be validated `online <https://bids-standard.github.io/bids-validator/>`_.

.. However, seapipe also works almost symbiotically with the `Wonambi <https://wonambi-python.github.io/>`_ package.
.. Therefore, any annotations (sleep scoring, artefact markings etc.) need to be inside a wonambi annotations file. 
.. For more information, see :doc:`Wonambi`

.. To receive an overview of your dataset, including whether the each participant's directory is BIDS compatible, as well as 
.. how many sessions, recording (e.g. edfs) and annoation files they contain, you can call the property of every dataset:
.. ``pipeline.audit`` 
.. ::
..    Participants   BIDS?	   #sessions #recordings   #annotations
..    sub-002        TRUE	      2	         2	            2
..    sub-004        TRUE	      3	         3	            3
..    sub-006	  TRUE	      2	         2	            1
..    sub-007	  TRUE	      2	         2	            2
..    sub-008	  TRUE	      2	         1	            1
..    sub-009	  TRUE	      2	         2	            2
..    sub-011	  TRUE	      0	         0	            0
..    sub-013	  TRUE	      2	         2	            2
..    sub-014	  FALSE       0          2                  2
..    sub-015	  TRUE	      2	         2	            2
..    sub-016	  TRUE	      2	         2	            0


.. This will be automatically saved to a file *dataset_audit.csv*

.. To retrieve a list of all the files inside the root directory, along with the
.. directories 1 and 2 levels preceding the files,
.. you can use the ``pipeline.list_dataset()`` function:

.. >>> pipeline.list_dataset()

.. :: 

..    Directory: project/bids
..    Files = ['dataset_description.json', 'participants.tsv']
..    ----------
..    Directory: ses-01/eeg
..    Files = ['sub-001_ses-01_eeg.edf']
..    ----------
..    Directory: ses-02/eeg
..    Files = ['sub-001_ses-02_eeg.edf']
..    ----------
..    Directory: ses-01/eeg
..    Files = ['sub-002_ses-01_eeg.edf']
..    ----------
..    etc.

.. The ``kind`` parameter should be either ``"meat"``, ``"fish"``,
.. or ``"veggies"``. Otherwise, :py:func:`lumache.get_random_ingredients`
.. will raise an exception.

.. .. autoexception:: lumache.InvalidKindError

.. For example:

.. >>> import lumache
.. >>> lumache.get_random_ingredients()
.. ['shells', 'gorgonzola', 'parsley']

