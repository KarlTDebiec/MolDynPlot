# -*- coding: utf-8 -*-
#   moldynplot.CpptrajDataset.py
#
#   Copyright (C) 2015 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Manages cpptraj datasets.
"""
################################### MODULES ###################################
from __future__ import absolute_import,division,print_function,unicode_literals
from .myplotspec.Dataset import Dataset
################################### CLASSES ###################################
class CpptrajDataset(Dataset):

    def __init__(self, verbose=1, debug=0, **kwargs):
        import numpy as np
        import pandas as pd

        super(CpptrajDataset, self).__init__(verbose=verbose, debug=debug,
          **kwargs)
        dataframe = self.dataframe
        dataframe.index.name = "time"

        # Convert from frame index to time
        if "dt" in kwargs:
            dt = kwargs.pop("dt")
            dataframe.index *= dt

        # Offset
        if "toffset" in kwargs:
            toffset = kwargs.pop("toffset")
            dataframe.index += toffset

        # Downsample
        if kwargs.get("downsample") is not None:
            downsample = kwargs.pop("downsample")
            downsample_mode = kwargs.get("downsample_mode", "mean")

            reduced = dataframe.values[
              :dataframe.shape[0]-(dataframe.shape[0] % downsample),:]
            new_shape=(int(reduced.shape[0]/ downsample),
                downsample, reduced.shape[1])
            index = np.reshape(dataframe.index.values[
              :dataframe.shape[0]-(dataframe.shape[0] % downsample)],
              new_shape[:-1]).mean(axis=1)
            reduced = np.reshape(reduced, new_shape)
            reduced = np.squeeze(reduced.mean(axis=1))
            reduced = pd.DataFrame(data=reduced, index=index,
              columns=dataframe.columns.values)
            reduced.index.name = "time"
            dataframe = self.dataframe = reduced

        if kwargs.get("pdist", False):
            from sklearn.neighbors import KernelDensity

            kde_kw = kwargs.get("kde_kw",  {"bandwidth": 0.1})
            grid = kwargs.get("grid", np.linspace(
                dataframe.values.min(), dataframe.values.max(), 100))
            kde = KernelDensity(**kde_kw)
            kde.fit(dataframe["rmsd"][:, np.newaxis])
            pdf = np.exp(kde.score_samples(grid[:, np.newaxis]))
            pdf /= pdf.sum()
            self.pdist_x = grid
            self.pdist_y = pdf

class NatConDataset(CpptrajDataset):

    def __init__(self, downsample=None, pdist=True,
        verbose=1, debug=0, **kwargs):
        import numpy as np
        import pandas as pd

        # Load
        super(NatConDataset, self).__init__(verbose=verbose, debug=debug,
          **kwargs)
        dataframe = self.dataframe
        n_contacts = self.dataframe.shape[1]

        # Convert minimum distances to percent native contacts
        cutoff = kwargs.get("cutoff", 5.5)
        percent = pd.DataFrame(data=(dataframe.values <= cutoff).sum(axis=1)
          / dataframe.shape[1], index=dataframe.index,
          columns=["percent_native_contacts"])
        dataframe = self.dataframe = percent

        # Downsample; flag included in function definition to prevent
        #   superclass from downsampling before applying cutoff
        if downsample is not None:
            from scipy.stats.mstats import mode

            reduced = dataframe.values[
              :dataframe.shape[0]-(dataframe.shape[0] % downsample),:]
            new_shape=(int(reduced.shape[0]/ downsample),
                downsample, reduced.shape[1])
            index = np.reshape(dataframe.index.values[
              :dataframe.shape[0]-(dataframe.shape[0] % downsample)],
              new_shape[:-1]).mean(axis=1)
            reduced = np.reshape(reduced, new_shape)
            reduced = np.squeeze(mode(reduced, axis=1)[0])
            reduced = pd.DataFrame(data=reduced, index=index,
              columns=dataframe.columns.values)
            reduced.index.name = "time"
            dataframe = self.dataframe = reduced

        # Calculate probability distribution
        if pdist:
            bins = np.linspace(0-((1/n_contacts)/2), 1+((1/n_contacts)/2),
              n_contacts+1)
            pdist, _ = np.histogram(self.dataframe.values, bins)
            pdist    =  np.array(pdist, np.float) / pdist.sum()
            pdist_x = np.zeros(bins.size*2)
            pdist_y = np.zeros(bins.size*2)
            pdist_x[::2]    = pdist_x[1::2]   = bins
            pdist_y[1:-1:2] = pdist_y[2:-1:2] = pdist
            self.pdist_x = pdist_x
            self.pdist_y = pdist_y
