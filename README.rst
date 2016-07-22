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

TimeSeriesFigureManager.py
--------------------------

Root-mean standard deviation:

.. only:: not html

    .. image:: doc/_static/p53/rmsd.png

.. only:: html

    test

    .. image:: _static/p53/rmsd.png

Radius of gyration:

.. image:: doc/_static/p53/radgyr.png

TimeSeries2DFigureManager.py
----------------------------

Per-residue Root-mean standard deviation:

.. image:: doc/_static/p53/perresrmsd.png

Dictionary of secondary structure prediction:

.. image:: doc/_static/p53/dssp.png

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
