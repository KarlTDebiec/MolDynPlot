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
    from moldynplot.Dataset import TimeSeriesDataset

    # Read text
    text_df = TimeSeriesDataset(infile="test/data/p53/rmsd.dat").timeseries_df
    assert hash(text_df.to_string()) == 8096571869640597227

    # Read hdf5
    hdf5_df = TimeSeriesDataset(infile="test/data/p53/rmsd.h5").timeseries_df
    assert hash(hdf5_df.to_string()) == 8096571869640597227
