Welcome to the seapipe documentation!
===================================

.. image:: https://img.shields.io/github/tag/nathanecross/seapipe?include_prereleases=&sort=semver&color=informational
  :target: https://github.com/nathanecross/seapipe/releases/
  :alt: version

.. image:: https://img.shields.io/github/issues/nathanecross/seapipe
  :target: https://github.com/nathanecross/seapipe/issues
  :alt: GitHub issues

.. image:: https://img.shields.io/badge/License-MIT-informational
  :target: https://opensource.org/license/mit/
  :alt: License: MIT

.. image:: https://app.codacy.com/project/badge/Grade/d9dd7acc4601421e8d52aab5015404c9
  :target: https://app.codacy.com/gh/nathanecross/seapipe/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_grade
  :alt: Codacy Badge

**seapipe** is a Python *pipeline* for researchers intending to investigate sleep neurophysiology.
The package reads EEG recordings with the purpose of detecting and analysing sleep events, such as *sleep architecture statistics*, *sleep spindles*, *slow oscillations*, *phase amplitude coupling*, and *power spectral analyses*. 

.. role:: raw-html(raw)
    :format: html 

It pulls functions from a range of other open source packages, including: `Wonambi <https://wonambi-python.github.io/>`_, the `MNE <https://mne.tools/stable/index.html>`_ and `Tensorpac <https://github.com/EtienneCmb/tensorpac>`_ 
  

It offers a *simple* and *intuitive*, yet *fully customisable* API. The major advantage of **seapipe** is that it can be deployed on full, diverse datasets that might require flexibility in certain parameters (e.g. electrode names, polarity of recordings, individualised frequency bands etc). The only requirement is that the data should be structured in `BIDS <https://bids-specification.readthedocs.io/en/stable/>`_.



Check out the :doc:`usage` section for further information, including
how to :ref:`installation` the project.

.. note::

   This project is under active development.

Contents
--------

.. toctree::

   usage
   api


