#!/usr/bin/python
# -*- coding: utf-8 -*-
#   test_Dataset.py
#
#   Copyright (C) 2015-2016 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
################################### MODULES ###################################
from pandas.util.testing import assert_frame_equal
from moldynplot.Dataset import SequenceDataset, TimeSeriesDataset
################################## FUNCTIONS ##################################
def test_sequence():
    # Read text
    text_df = SequenceDataset(
      infile="data/gb3/relax_s2.dat").sequence_df

    # Read hdf5
    hdf5_df = SequenceDataset(
      infile="data/gb3/relax_s2.h5").sequence_df

    # Merge text
    merged_text_df = SequenceDataset(
      infiles=["data/gb3/relax.dat", "data/gb3/s2.dat"]).sequence_df

    # Merge hdf5
    merged_hdf5_df = SequenceDataset(
      infiles=["data/gb3/relax.h5", "data/gb3/s2.h5"]).sequence_df

    # Compare
    assert_frame_equal(text_df, hdf5_df)
    assert_frame_equal(text_df, merged_text_df)
    assert_frame_equal(hdf5_df, merged_hdf5_df)

    # Write text

    # Write hdf5

def test_rmsd():
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

    # Write text

    # Write hdf5

def test_radgyr():
    # Read cpptraj
    cpptraj_df = TimeSeriesDataset(
      infile="data/p53/radgyr.cpptraj",
      dt=0.1,
      toffset=-0.1).timeseries_df

    # Read text
    text_df = TimeSeriesDataset(
      infile="data/p53/radgyr.dat").timeseries_df

    # Read hdf5
    hdf5_df = TimeSeriesDataset(
      infile="data/p53/radgyr.h5").timeseries_df

    # Compare
    assert_frame_equal(cpptraj_df, hdf5_df)
    assert_frame_equal(cpptraj_df, text_df)

    # Write text

    # Write hdf5

def test_perresrmsd():
    # Read cpptraj
    cpptraj_df = TimeSeriesDataset(
      infile="data/p53/perresrmsd.cpptraj",
      dt=0.1,
      toffset=-0.1).timeseries_df

    # Read text
    text_df = TimeSeriesDataset(
      infile="data/p53/perresrmsd.dat").timeseries_df

    # Read hdf5
    hdf5_df = TimeSeriesDataset(
      infile="data/p53/perresrmsd.h5").timeseries_df

    # Read legacy hdf5
    lgcy_df = TimeSeriesDataset(
      infile="data/p53/perresrmsd_legacy.h5",
      dt=0.1).timeseries_df

    # Compare
    assert_frame_equal(cpptraj_df, hdf5_df)
    assert_frame_equal(cpptraj_df, text_df)
    assert_frame_equal(cpptraj_df, lgcy_df)

    # Write text

    # Write hdf5

def test_dssp():
    # Read cpptraj
    cpptraj_df = TimeSeriesDataset(
      infile="data/p53/dssp.cpptraj",
      dt=0.1,
      toffset=-0.1).timeseries_df

    # Read text
    text_df = TimeSeriesDataset(
      infile="data/p53/dssp.dat").timeseries_df

    # Read hdf5
    hdf5_df = TimeSeriesDataset(
      infile="data/p53/dssp.h5").timeseries_df

    # Read legacy hdf5
    lgcy_df = TimeSeriesDataset(
      infile="data/p53/dssp_legacy.h5",
      dt=0.1).timeseries_df

    # Compare
    assert_frame_equal(cpptraj_df, hdf5_df)
    assert_frame_equal(cpptraj_df, text_df)
    assert_frame_equal(cpptraj_df, lgcy_df)

    # Write text

    # Write hdf5

if __name__ == "__main__":
    test_sequence()
    test_rmsd()
    test_perresrmsd()
    test_dssp()
