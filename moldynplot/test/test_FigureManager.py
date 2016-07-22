#!/usr/bin/python
# -*- coding: utf-8 -*-
#   test_FigureManager.py
#
#   Copyright (C) 2015-2016 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
################################### MODULES ###################################
from matplotlib.testing.compare import compare_images
from moldynplot.SequenceFigureManager import SequenceFigureManager
from moldynplot.TimeSeriesFigureManager import TimeSeriesFigureManager
from moldynplot.TimeSeries2DFigureManager import TimeSeries2DFigureManager
################################## FUNCTIONS ##################################
def test_relax():
    sfm = SequenceFigureManager()
    sfm.draw_report(yaml_spec="yaml/gb3/relax.yml")
    compare_images("relax.png", "figure/gb3/relax.png", tol=0)

def test_rmsd():
    tsfm = TimeSeriesFigureManager()
    tsfm.draw_report(yaml_spec="yaml/p53/rmsd.yml")
    compare_images("rmsd.png", "figure/p53/rmsd.png", tol=0)

def test_radgyr():
    tsfm = TimeSeriesFigureManager()
    tsfm.draw_report(yaml_spec="yaml/p53/radgyr.yml")
    compare_images("radgyr.png", "figure/p53/radgyr.png", tol=0)

def test_perresrmsd():
    tsfm = TimeSeries2DFigureManager()
    tsfm.draw_report(yaml_spec="yaml/p53/perresrmsd.yml")
    compare_images("perresrmsd.png", "figure/p53/perresrmsd.png", tol=0)

def test_dssp():
    tsfm = TimeSeries2DFigureManager()
    tsfm.draw_report(yaml_spec="yaml/p53/dssp.yml")
    compare_images("dssp.png", "figure/p53/dssp.png", tol=0)

if __name__ == "__main__":
    test_relax()
    test_rmsd()
    test_radgyr()
    test_perresrmsd()
    test_dssp()
