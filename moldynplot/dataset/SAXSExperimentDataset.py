#!/usr/bin/python
# -*- coding: utf-8 -*-
#   moldynplot.dataset.SAXSExperimentDataset.py
#
#   Copyright (C) 2015-2017 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Represents experimental Small Angle X-ray Scattering (SAXS) data
"""
################################### MODULES ###################################
from __future__ import (absolute_import, division, print_function,
    unicode_literals)

if __name__ == "__main__":
    __package__ = str("moldynplot.dataset")
    import moldynplot.dataset
import numpy as np
import pandas as pd
from .SAXSDataset import SAXSDataset
from ..myplotspec import wiprint


################################### CLASSES ###################################
class SAXSExperimentDataset(SAXSDataset):
    """
    Represents experimental Small Angle X-ray Scattering (SAXS) data
    """

    def __init__(self, scale=False, **kwargs):
        """
        Arguments:
          infile (str): Path to input file, may contain environment
            variables
          verbose (int): Level of verbose output
          kwargs (dict): Additional keyword arguments
        """

        # Load
        super(SAXSExperimentDataset, self).__init__(**kwargs)
        self.df = self.dataframe
        del self.dataframe

        # Scale
        if scale:
            self.scale(scale, **kwargs)


#################################### MAIN #####################################
if __name__ == "__main__":
    SAXSExperimentDataset.main()