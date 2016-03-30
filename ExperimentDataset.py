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

    def __init__(self, infile, log=False, yoffset=None,
        verbose=1, debug=0, **kwargs):
        from os.path import expandvars
        import h5py
        import numpy as np
        import pandas as pd

        # Load
        super(SAXSDataset, self).__init__(infile=infile,
          verbose=verbose, debug=debug, **kwargs)
        dataframe = self.dataframe

#        # Store y
#        self.y = np.array(dataframe.columns, np.float)

        # Log scale
        if log:
            self.dataframe["intensity_se"] = (
             dataframe["intensity_se"] / (self.dataframe["intensity"] *
             np.log(10)))
            self.dataframe["intensity"] = np.log10(dataframe["intensity"])

        # Offset
        if yoffset is not None:
            self.dataframe["intensity"] += yoffset
        print(self.dataframe)
