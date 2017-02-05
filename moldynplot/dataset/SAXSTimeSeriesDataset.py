#!/usr/bin/python
# -*- coding: utf-8 -*-
#   moldynplot.dataset.SAXSTimeSeriesDataset.py
#
#   Copyright (C) 2015-2017 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Manages Small Angle X-ray Scattering (SAXS) time series datas
"""
################################### MODULES ###################################
from __future__ import (absolute_import, division, print_function,
    unicode_literals)

if __name__ == "__main__":
    __package__ = str("moldynplot.dataset")
    import moldynplot.dataset
from IPython import embed
import h5py
import numpy as np
import pandas as pd
import six
from .SAXSDataset import SAXSDataset
from .TimeSeriesDataset import TimeSeriesDataset
from ..myplotspec import wiprint


################################### CLASSES ###################################
class SAXSTimeSeriesDataset(TimeSeriesDataset, SAXSDataset):
    """
    Manages Small Angle X-ray Scattering (SAXS) time series datas
    """

    def __init__(self, infile, address="saxs", dt=None, toffset=None,
      downsample=None, calc_mean=False, calc_error=True, error_method="std",
      scale=False, outfile=None, interactive=False, **kwargs):
        """
        Arguments:
          infile (str): Path to input file, may contain environment
            variables
          usecols (list): Columns to select from DataFrame, once
            dataframe has already been loaded
          dt (float): Time interval between points; units unspecified
          toffset (float): Time offset to be added to all points (i.e.
            time of first point)
          downsample (int): Interval by which to downsample points
          verbose (int): Level of verbose output
          kwargs (dict): Additional keyword arguments

        .. todo:
          - Calculate error
          - Shift downsampling to superclass
        """
        from os.path import expandvars

        # Arguments
        verbose = kwargs.get("verbose", 1)
        self.dataset_cache = kwargs.get("dataset_cache", None)

        # Read data
        with h5py.File(expandvars(infile)) as h5_file:
            q = ["{0:5.3f}".format(a) for a in
                np.array(h5_file[address + "/q"])]
        self.timeseries_df = self.read(
          infile=infile + ":/" + address + "/intensity",
          dataframe_kw=dict(columns=q), **kwargs)

        # Process data
        if dt:
            self.timeseries_df.set_index(
              self.timeseries_df.index.values * float(dt), inplace=True)
            self.timeseries_df.index.name = "time"
        if toffset:
            index_name = self.timeseries_df.index.name
            self.timeseries_df.set_index(
              self.timeseries_df.index.values + float(toffset), inplace=True)
            self.timeseries_df.index.name = index_name
        if downsample:
            self.downsample(downsample, downsample_mode="mean", **kwargs)

        # Output data
        if verbose >= 2:
            print("Processed timeseries DataFrame:")
            print(self.timeseries_df)
        if outfile is not None:
            self.write(df=self.timeseries_df, outfile=outfile, **kwargs)

        # Calculate mean and standard error
        if calc_mean:
            block_kw = dict(min_n_blocks=2, max_cut=0.1, all_factors=False,
              fit_exp=True, fit_sig=False)
            #            block_kw = dict(min_n_blocks=2, max_cut=0.1,
            # all_factors=True,
            #                            fit_exp=True, fit_sig=False)
            block_kw.update(kwargs.get("block_kw", {}))
            self.mean_df, self.block_averager = self.calc_mean(
              df=self.timeseries_df, mode="se", verbose=verbose, **block_kw)
            self.mean_df.index.name = "q"
            self.mean_df.columns = ["intensity", "intensity se"]
            if verbose >= 2:
                print("Processed mean DataFrame:")
                print(self.mean_df)
            if isinstance(calc_mean, six.string_types):
                self.write(df=self.mean_df, outfile=calc_mean, **kwargs)
            self.df = self.mean_df

        # Scale
        if scale:
            self.scale(scale, **kwargs)


#################################### MAIN #####################################
if __name__ == "__main__":
    SAXSTimeSeriesDataset.main()