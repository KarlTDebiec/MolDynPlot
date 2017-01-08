#!/usr/bin/python
# -*- coding: utf-8 -*-
#   moldynplot.dataset.NatConTimeSeriesDataset.py
#
#   Copyright (C) 2015-2017 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Represents native contacts as a function of time
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
from .TimeSeriesDataset import TimeSeriesDataset
from ..myplotspec import wiprint


################################### CLASSES ###################################
class NatConTimeSeriesDataset(TimeSeriesDataset):
    """
    Represents native contacts as a function of time
    """

    def __init__(self, downsample=None, calc_pdist=True, **kwargs):
        """
        Arguments:
          infile (str): Path to input file, may contain environment
            variables
          usecols (list): Columns to select from DataFrame, once
            dataframe has already been loaded
          dt (float): Time interval between points; units unspecified
          toffset (float): Time offset to be added to all points (i.e.
            time of first point)
          cutoff (float): Minimum distance within which a contact is
            considered to be formed
          downsample (int): Interval by which to downsample points using
            mode
          pdist (bool): Calculate probability distribution
          verbose (int): Level of verbose output
          kwargs (dict): Additional keyword arguments
        """
        verbose = kwargs.get("verbose", 1)

        # Load
        super(NatConTimeSeriesDataset, self).__init__(**kwargs)
        dataframe = self.dataframe
        n_contacts = self.dataframe.shape[1]

        # Convert minimum distances to percent native contacts
        cutoff = kwargs.get("cutoff", 5.5)
        percent = pd.DataFrame(
          data=(dataframe.values <= cutoff).sum(axis=1) / dataframe.shape[1],
          index=dataframe.index, columns=["percent_native_contacts"])
        dataframe = self.dataframe = percent

        # Downsample; flag included in function definition to prevent
        #   superclass from downsampling before applying cutoff
        if downsample is not None:
            from scipy.stats.mstats import mode

            if verbose >= 1:
                print("downsampling by factor of {0} using mode".format(
                  downsample))

            reduced = dataframe.values[
            :dataframe.shape[0] - (dataframe.shape[0] % downsample), :]
            new_shape = (int(reduced.shape[0] / downsample), downsample,
            reduced.shape[1])
            index = np.reshape(dataframe.index.values[
            :dataframe.shape[0] - (dataframe.shape[0] % downsample)],
              new_shape[:-1]).mean(axis=1)
            reduced = np.reshape(reduced, new_shape)
            reduced = np.squeeze(mode(reduced, axis=1)[0])
            reduced = pd.DataFrame(data=reduced, index=index,
              columns=dataframe.columns.values)
            reduced.index.name = "time"
            dataframe = self.dataframe = reduced

        # Calculate probability distribution
        if calc_pdist:
            if verbose >= 1:
                print("calculating probability distribution using histogram")
            bins = np.linspace(0 - ((1 / n_contacts) / 2),
              1 + ((1 / n_contacts) / 2), n_contacts + 1)
            pdist, _ = np.histogram(self.dataframe.values, bins)
            pdist = np.array(pdist, np.float) / pdist.sum()
            pdist_x = np.zeros(bins.size * 2)
            pdist_y = np.zeros(bins.size * 2)
            pdist_x[::2] = pdist_x[1::2] = bins
            pdist_y[1:-1:2] = pdist_y[2:-1:2] = pdist
            self.pdist_x = pdist_x
            self.pdist_y = pdist_y

        self.timeseries = dataframe


#################################### MAIN #####################################
if __name__ == "__main__":
    NatConTimeSeriesDataset.main()