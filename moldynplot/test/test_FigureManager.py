#!/usr/bin/python
# -*- coding: utf-8 -*-
#   test_FigureManager.py
#
#   Copyright (C) 2015-2016 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.

def test_TimeSeriesFigureManager():
    from pandas.util.testing import assert_frame_equal
    from moldynplot.Dataset import TimeSeriesDataset

    # Read cpptraj
    cpptraj_df = TimeSeriesDataset(
      infile="data/p53/rmsd.cpptraj",
      dt=0.1,
      toffset=-0.1).timeseries_df

    # Read text
    text_df = TimeSeriesDataset(
      infile="data/p53/rmsd.dat").timeseries_df

    # Read hdf5
    hdf5_df = TimeSeriesDataset(
      infile="data/p53/rmsd.h5").timeseries_df

    # Compare
    assert_frame_equal(cpptraj_df, hdf5_df)
    assert_frame_equal(cpptraj_df, text_df)
    assert_frame_equal(text_df, hdf5_df)

    # Write text

    # Write hdf5
