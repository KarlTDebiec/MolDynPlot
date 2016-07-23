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
from filecmp import cmp
import numpy as np
from pandas.util.testing import assert_frame_equal
from moldynplot.Dataset import (HSQCDataset, SequenceDataset,
                                TimeSeriesDataset)
################################## FUNCTIONS ##################################
def h5_cmp(file_1, file_2):
    from os import devnull
    from subprocess import Popen, PIPE
    with open(devnull, "w") as fnull:
        output = Popen("h5diff {0} {1}".format(file_1, file_2),
            stdout=PIPE, stderr=fnull, shell=True).stdout.read().strip()
    if output == "":
        return True
    else:
        return False

#################################### TESTS ####################################
def test_hsqc():
    # Read NMRPipe
    pipe = HSQCDataset(infile="data/mocvnh3/hsqc.ft")

    # Read text
    text = HSQCDataset(infile="data/mocvnh3/hsqc.dat",)
    assert_frame_equal(pipe.hsqc_df, text.hsqc_df)

    # Read hdf5
    hdf5 = HSQCDataset(infile="data/mocvnh3/hsqc.h5")
    assert_frame_equal(pipe.hsqc_df, hdf5.hsqc_df)

    # Write text
    pipe.write(dataframe=pipe.hsqc_df, outfile="hsqc.dat")
    assert(cmp("hsqc.dat", "data/mocvnh3/hsqc.dat") == True)

    # Write hdf5
    pipe.write(dataframe=pipe.hsqc_df, outfile="hsqc.h5")
    assert(h5_cmp("hsqc.h5", "data/mocvnh3/hsqc.h5") == True)

def test_sequence():
    # Read text
    text = SequenceDataset(infile="data/gb3/relax_s2.dat")

    # Read hdf5
    hdf5 = SequenceDataset(infile="data/gb3/relax_s2.h5")
    assert_frame_equal(text.sequence_df, hdf5.sequence_df)

    # Merge text
    merged_text = SequenceDataset(
      infiles=["data/gb3/relax.dat", "data/gb3/s2.dat"])
    assert_frame_equal(text.sequence_df, merged_text.sequence_df)

    # Merge hdf5
    merged_hdf5 = SequenceDataset(
      infiles=["data/gb3/relax.h5", "data/gb3/s2.h5"])
    assert_frame_equal(text.sequence_df, merged_hdf5.sequence_df)

    # Write text
    text.write(dataframe=text.sequence_df, outfile="relax_s2.dat")
    assert(cmp("relax_s2.dat", "data/gb3/relax_s2.dat") == True)

    # Write hdf5
    text.write(dataframe=text.sequence_df, outfile="relax_s2.h5")
    assert(h5_cmp("relax_s2.h5", "data/gb3/relax_s2.h5") == True)

def test_rmsd():
    # Read cpptraj
    cpptraj = TimeSeriesDataset(infile="data/p53/rmsd.cpptraj",
      dt=0.1, toffset=-0.1)

    # Read text
    text = TimeSeriesDataset(infile="data/p53/rmsd.dat")
    assert_frame_equal(cpptraj.timeseries_df, text.timeseries_df)

    # Read hdf5
    hdf5 = TimeSeriesDataset(infile="data/p53/rmsd.h5")
    assert_frame_equal(cpptraj.timeseries_df, hdf5.timeseries_df)

    # Write text
    cpptraj.write(dataframe=text.timeseries_df, outfile="rmsd.dat")
    assert(cmp("rmsd.dat", "data/p53/rmsd.dat") == True)

    # Write hdf5
    cpptraj.write(dataframe=text.timeseries_df, outfile="rmsd.h5")
    assert(h5_cmp("rmsd.h5", "data/p53/rmsd.h5") == True)

def test_radgyr():
    # Read cpptraj
    cpptraj = TimeSeriesDataset(infile="data/p53/radgyr.cpptraj",
      dt=0.1, toffset=-0.1)

    # Read text
    text = TimeSeriesDataset(infile="data/p53/radgyr.dat")
    assert_frame_equal(cpptraj.timeseries_df, text.timeseries_df)

    # Read hdf5
    hdf5 = TimeSeriesDataset(infile="data/p53/radgyr.h5")
    assert_frame_equal(cpptraj.timeseries_df, hdf5.timeseries_df)

    # Write text
    cpptraj.write(dataframe=text.timeseries_df, outfile="radgyr.dat")
    assert(cmp("radgyr.dat", "data/p53/radgyr.dat") == True)

    # Write hdf5
    cpptraj.write(dataframe=text.timeseries_df, outfile="radgyr.h5")
    assert(h5_cmp("radgyr.h5", "data/p53/radgyr.h5") == True)

def test_perresrmsd():
    # Read cpptraj
    cpptraj = TimeSeriesDataset(infile="data/p53/perresrmsd.cpptraj",
      dt=0.1, toffset=-0.1, dtype=np.float32)

    # Read text
    text = TimeSeriesDataset(infile="data/p53/perresrmsd.dat",
      dtype=np.float32)
    assert_frame_equal(cpptraj.timeseries_df, text.timeseries_df)

    # Read hdf5
    hdf5 = TimeSeriesDataset(infile="data/p53/perresrmsd.h5")
    assert_frame_equal(cpptraj.timeseries_df, hdf5.timeseries_df)

    # Read legacy hdf5
    legacy = TimeSeriesDataset(infile="data/p53/perresrmsd_legacy.h5", dt=0.1)
    assert_frame_equal(cpptraj.timeseries_df, legacy.timeseries_df)

    # Write text
    cpptraj.write(dataframe=text.timeseries_df, outfile="perresrmsd.dat")
    assert(cmp("perresrmsd.dat", "data/p53/perresrmsd.dat") == True)

    # Write hdf5
    cpptraj.write(dataframe=text.timeseries_df, outfile="perresrmsd.h5")
    assert(h5_cmp("perresrmsd.h5", "data/p53/perresrmsd.h5") == True)

def test_dssp():
    # Read cpptraj
    cpptraj = TimeSeriesDataset(infile="data/p53/dssp.cpptraj", dt=0.1,
      toffset=-0.1, dtype=np.uint8)

    # Read text
    text = TimeSeriesDataset(infile="data/p53/dssp.dat", dtype=np.uint8)
    assert_frame_equal(cpptraj.timeseries_df, text.timeseries_df)

    # Read hdf5
    hdf5 = TimeSeriesDataset(infile="data/p53/dssp.h5")
    assert_frame_equal(cpptraj.timeseries_df, hdf5.timeseries_df)

    # Read legacy hdf5
    legacy = TimeSeriesDataset(infile="data/p53/dssp_legacy.h5", dt=0.1)
    assert_frame_equal(cpptraj.timeseries_df, legacy.timeseries_df)

    # Write text
    cpptraj.write(dataframe=text.timeseries_df, outfile="dssp.dat")
    assert(cmp("dssp.dat", "data/p53/dssp.dat") == True)

    # Write hdf5
    cpptraj.write(dataframe=text.timeseries_df, outfile="dssp.h5")
    assert(h5_cmp("dssp.h5", "data/p53/dssp.h5") == True)

if __name__ == "__main__":
    test_sequence()
    test_rmsd()
    test_radgyr()
    test_perresrmsd()
    test_dssp()
    test_hsqc()
