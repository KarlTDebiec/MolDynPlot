.. image:: https://travis-ci.org/KarlTDebiec/MolDynPlot.svg?branch=master
    :target: https://travis-ci.org/KarlTDebiec/MolDynPlot

.. image:: https://coveralls.io/repos/github/KarlTDebiec/MolDynPlot/badge.svg?branch=master
    :target: https://coveralls.io/github/KarlTDebiec/MolDynPlot?branch=master

Introduction
============

MolDynPlot is a Python package used to design matplotlib-based figures of
Molecular Dynamics (MD) simulation data using the simple text format YAML.

Example Output
==============

.. github_start

Time Series:

.. image:: doc/_static/p53/rmsd.png

2D Time Series:

.. image:: doc/_static/p53/dssp.png

Sequence:

.. image:: doc/_static/gb3/s2.png

NMR HSQC:

.. image:: doc/_static/mocvnh3/hsqc.png

.. github_end

.. only:: html

    Time Series:

    .. image:: _static/p53/rmsd.png
        :scale: 35
        :align: center

    2D Time Series:

    .. image:: _static/p53/dssp.png
        :scale: 35
        :align: center

    Sequence:

    .. image:: _static/gb3/s2.png
        :scale: 35
        :align: center
    NMR HSQC:

    .. image:: _static/mocvnh3/hsqc.png
        :scale: 35
        :align: center

Publications Prepared Using MolDynPlot
--------------------------------------

- `Further along the Road Less Traveled: AMBER ff15ipq, an Original Protein
  Force Field Built on a Self-Consistent Physical Model
  <http://pubs.acs.org/doi/abs/10.1021/acs.jctc.6b00567>`_

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
