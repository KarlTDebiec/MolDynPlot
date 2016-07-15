#!/usr/bin/python
# -*- coding: utf-8 -*-
#   test_timeseries.py
#
#   Copyright (C) 2015-2016 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
def test_TimeSeriesDataset():
    from pandas.util.testing import assert_frame_equal
    from moldynplot.Dataset import TimeSeriesDataset
    text_df = TimeSeriesDataset(infile="test/data/p53/rmsd.dat").timeseries_df
    hdf5_df = TimeSeriesDataset(infile="test/data/p53/rmsd.h5").timeseries_df
    assert_frame_equal(text_df, hdf5_df)
