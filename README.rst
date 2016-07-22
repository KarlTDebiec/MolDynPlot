.. image:: https://travis-ci.org/KarlTDebiec/MolDynPlot.svg?branch=master
    :target: https://travis-ci.org/KarlTDebiec/MolDynPlot

.. image:: https://coveralls.io/repos/github/KarlTDebiec/MolDynPlot/badge.svg?branch=master
    :target: https://coveralls.io/github/KarlTDebiec/MolDynPlot?branch=master

Introduction
============

MolDynPlot is a Python package used to design matplotlib-based figures of
Molecular Dynamics (MD) simulation data using the simple text format YAML.

Example Outputs
===============

.. github_start

TimeSeriesFigureManager.py
--------------------------

Root-mean standard deviation:

.. image:: doc/_static/p53/rmsd.png
    :scale: 10
    :align: center

.. raw:: html

   <big>text</big>

Radius of gyration:

.. image:: doc/_static/p53/radgyr.png
    :scale: 10
    :align: center

TimeSeries2DFigureManager.py
----------------------------

Per-residue Root-mean standard deviation:

.. image:: doc/_static/p53/perresrmsd.png
    :scale: 10
    :align: center

Dictionary of secondary structure prediction:

.. image:: doc/_static/p53/dssp.png
    :scale: 10
    :align: center

.. github_end

.. only:: html

    TimeSeriesFigureManager.py
    --------------------------

    Root-mean standard deviation:

    .. image:: _static/p53/rmsd.png
        :scale: 35
        :align: center

    Radius of gyration:

    .. image:: _static/p53/radgyr.png
        :scale: 35
        :align: center

    TimeSeries2DFigureManager.py
    ----------------------------

    Per-residue Root-mean standard deviation:

    .. image:: _static/p53/perresrmsd.png
        :scale: 35
        :align: center

    Dictionary of secondary structure prediction:

    .. image:: _static/p53/dssp.png
        :scale: 35
        :align: center

Dependencies
------------

MolDynPlot supports Python 2.7 and 3.4, and requires the following
packages:

- h5py
- ipython
- matplotlib
- pandas
- pyyaml
- scikit-learn
- scipy

Installation
------------

Change to the ``MolDynPlot`` directory and run::

    python setup.py install

Authorship
----------

MolDynPlot is developed by Karl T. Debiec, a graduate student at the
University of Pittsburgh advised by Professors Lillian T. Chong and Angela M.
Gronenborn.

License
-------

Released under a 3-clause BSD license.
