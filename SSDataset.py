# -*- coding: utf-8 -*-
#   moldynplot.SSDataset.py
#
#   Copyright (C) 2015 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Manages secondary structure datasets.
"""
################################### MODULES ###################################
from __future__ import absolute_import,division,print_function,unicode_literals
from .myplotspec.Dataset import Dataset
################################### CLASSES ###################################
class SSDataset(Dataset):
    """
    Manages secondary structure datasets.

    Input data should be providied in a whitespace-delimited text file
    output from `cpptraj`'s `secstruct` command::

      #Residue   Para   Anti   3-10  Alpha     Pi   Turn   Bend
             1 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
           ...    ...    ...    ...    ...    ...    ...    ...
    """

    def __init__(self, verbose=1, debug=0, **kwargs):
        """
        Arguments:
          infile (str): Path to text input file, may contain environment
            variables
          verbose (int): Level of verbose output
          debug (int): Level of debug output
          kwargs (dict): Additional keyword arguments
        """
        import numpy as np
        from .myplotspec import multi_get_copy

        # Manage arguments
        read_csv_kw = dict(delim_whitespace=True, index_col=0)
        read_csv_kw.update(kwargs.pop("read_csv_kw", {}))

        # Load data
        dataframe = self.load_dataset(verbose=verbose, debug=debug,
          read_csv_kw=read_csv_kw, **kwargs).data

        self.data = dataframe
