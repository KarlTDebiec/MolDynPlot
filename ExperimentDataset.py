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
        yoffset=None,
        verbose=1, debug=0, **kwargs):
        from os.path import expandvars
        import h5py
        import numpy as np
        import pandas as pd

        # Load
        super(SAXSDataset, self).__init__(infile=infile,
          verbose=verbose, debug=debug, **kwargs)
        dataframe = self.dataframe

#        # Offset
#        if yoffset is not None:
#            self.dataframe["intensity"] += yoffset
