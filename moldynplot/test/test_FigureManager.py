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
from moldynplot.TimeSeriesFigureManager import TimeSeriesFigureManager
from moldynplot.TimeSeries2DFigureManager import TimeSeries2DFigureManager
################################## FUNCTIONS ##################################
def test_rmsd():
    tsfm = TimeSeriesFigureManager()
    tsfm.draw_report(yaml_spec="yaml/p53/rmsd.yml")
    compare_images("rmsd.png", "figure/p53/rmsd.png", tol=0)

def test_perresrmsd():
    tsfm = TimeSeries2DFigureManager()
    tsfm.draw_report(yaml_spec="yaml/p53/perresrmsd.yml")
    compare_images("perresrmsd.png", "figure/p53/perresrmsd.png", tol=0)

if __name__ == "__main__":
    test_rmsd()
    test_perresrmsd()
