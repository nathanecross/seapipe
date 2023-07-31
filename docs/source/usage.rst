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

To retrieve a list of random ingredients,
you can use the ``lumache.get_random_ingredients()`` function:

 lumache.get_random_ingredients

The ``kind`` parameter should be either ``"meat"``, ``"fish"``,
or ``"veggies"``. Otherwise, :py:func:`lumache.get_random_ingredients`
will raise an exception.

.. autoexception:: lumache.InvalidKindError

For example:

>>> import lumache
>>> lumache.get_random_ingredients()
['shells', 'gorgonzola', 'parsley']

