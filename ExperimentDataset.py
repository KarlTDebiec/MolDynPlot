# -*- coding: utf-8 -*-
#   moldynplot.ExperimentDataset.py
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
class SAXSDataset(Dataset):

    def __init__(self, infile,
        scale=False,
        verbose=1, debug=0, **kwargs):
        from os.path import expandvars
        import h5py
        import numpy as np
        import pandas as pd

        # Load
        super(SAXSDataset, self).__init__(infile=infile,
          verbose=verbose, debug=debug, **kwargs)
        dataframe = self.dataframe

        if scale:
            from scipy.interpolate import interp1d
            from scipy.optimize import curve_fit

            # Prepare target
            scale_target = expandvars(kwargs.pop("scale_target"))
            target = self.load_dataset(infile=scale_target, loose=True,
              **kwargs)
            target_x = target.dataframe.index.values
            target_y = target.dataframe["intensity"]

            # Prepare own values
            self_x = dataframe.index.values
            self_y = dataframe["intensity"].values
            indexes = np.logical_and(self_x > target_x.min(),
                                     self_x < target_x.max())
            self_x = self_x[indexes]
            self_y = self_y[indexes]

            # Must increase precision to support 
            self_x      = np.array(self_x, np.float64)
            self_y      = np.array(self_y, np.float64)
            target_y    = np.array(target_y, np.float64)

            # Update target
            interp_target_y    = interp1d(target_x, target_y, kind="cubic")
            target_y           = interp_target_y(self_x)

            def scale_y(_, a):
                return a * self_y
            scaling_factor = curve_fit(scale_y, self_x, target_y,
              p0=(1))[0][0]
            print(scaling_factor)
            self.dataframe["intensity"]    *= scaling_factor
            if "intensity_se" in self.dataframe.columns.values:
                self.dataframe["intensity_se"] *= scaling_factor
