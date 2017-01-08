#!/usr/bin/python
# -*- coding: utf-8 -*-
#   moldynplot.dataset.SAXSDiffDataset.py
#
#   Copyright (C) 2015-2017 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Represents Small-Angle X-ray Scattering (SAXS) data
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
class SAXSDiffDataset(SAXSDataset):
    """
    Represents difference Small Angle X-ray Scattering (SAXS) data.
    """

    def __init__(self, dataset_cache=None, **kwargs):
        """
        """
        from ..myplotspec import multi_get_copy

        self.dataset_cache = dataset_cache

        minuend_kw = multi_get_copy(["minuend", "minuend_kw"], kwargs, {})
        subtrahend_kw = multi_get_copy(["subtrahend", "subtrahend_kw"], kwargs,
          {})
        minuend = self.load_dataset(loose=True, **minuend_kw)
        subtrahend = self.load_dataset(loose=True, **subtrahend_kw)
        m_I = minuend.df["intensity"]
        m_I_se = minuend.df["intensity se"]
        s_I = subtrahend.df["intensity"]
        s_I_se = subtrahend.df["intensity se"]
        diff_I = (m_I - s_I)
        diff_I_se = np.sqrt(m_I_se ** 2 + s_I_se ** 2)
        diff_I.name = "intensity"
        diff_I_se.name = "intensity se"
        self.df = pd.concat([diff_I, diff_I_se], axis=1)


#################################### MAIN #####################################
if __name__ == "__main__":
    SAXSDiffDataset.main()