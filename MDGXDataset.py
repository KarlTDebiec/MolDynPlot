# -*- coding: utf-8 -*-
#   moldynplot.MDGXDataset.py
#
#   Copyright (C) 2015-2016 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Manages MDGX datasets
"""
################################### MODULES ###################################
from __future__ import absolute_import,division,print_function,unicode_literals
from .myplotspec.Dataset import Dataset
################################### CLASSES ###################################
class MDGXDataset(Dataset):
    """
    Manages MDGX datasets
    """

    @classmethod
    def get_cache_key(cls, infile, selections=None, *args, **kwargs):
        """
        Generates tuple of arguments to be used as key for dataset
        cache.

        .. todo:
          - Verify that keyword arguments passed to pandas may be safely
            converted to hashable tuple, and if they cannot throw a
            warning and load dataset without memoization
        """
        from os.path import expandvars

        read_csv_kw = []
        for key, value in kwargs.get("read_csv_kw", {}).items():
            if isinstance(value, list):
                value = tuple(value)
            read_csv_kw.append((key, value))
        return (cls, expandvars(infile), tuple(selections), tuple(read_csv_kw))

    def __init__(self, infile,
        selections=None,
        verbose=1, debug=0, **kwargs):
        from os.path import expandvars
        import numpy as np
        import pandas as pd
        pd.set_option('display.width', 1000)

        # Load
        super(MDGXDataset, self).__init__(infile=infile,
          verbose=verbose, debug=debug, **kwargs)
        dataframe = self.dataframe
        dataframe.index.name = "conformation"
        dataframe["error"] = np.abs(dataframe["qm_energy"]
          - dataframe["mm_energy"])

        selection_dataframes = []
        if selections is not None:
            for selection in selections:
                selection_dataframes.append(
                  dataframe[dataframe["topology"].str.endswith(
                  "/" + selection)])
        self.selections = selection_dataframes
