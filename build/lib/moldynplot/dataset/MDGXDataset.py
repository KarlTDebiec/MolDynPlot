#!/usr/bin/python
# -*- coding: utf-8 -*-
#   moldynplot.dataset.MDGXDataset.py
#
#   Copyright (C) 2015-2016 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
.. todo:
  - Bring up to date with other Dataset classes
  - Fix separation and ordering of argument groups: input, action, output
"""
################################### MODULES ###################################
from __future__ import absolute_import,division,print_function,unicode_literals
if __name__ == "__main__":
    __package__ = str("moldynplot")
    import moldynplot
from IPython import embed
import h5py
import numpy as np
import pandas as pd
from ..myplotspec.Dataset import Dataset
from ..myplotspec import sformat, wiprint
################################### CLASSES ###################################
class MDGXDataset(Dataset):
    """
    Represents MDGX force field parameterization data.
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

    def __init__(self, infile, selections=None, **kwargs):
        """
        """
        from os.path import expandvars

        # Load
        super(MDGXDataset, self).__init__(infile=infile, **kwargs)
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

#################################### MAIN #####################################
if __name__ == "__main__":
    import argparse

    # Prepare argument parser
    parser = argparse.ArgumentParser(
      description = """Processes datasets""")
    subparsers = parser.add_subparsers(
      dest        = "mode",
      description = "")

    MDGXDataset.construct_argparser(subparsers)

    kwargs  = vars(parser.parse_args())
    kwargs.pop("cls")(**kwargs)
