#!/usr/bin/python
# -*- coding: utf-8 -*-
#   moldynplot.dataset.CorrDataset.py
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
    __package__ = str("moldynplot.dataset")
    import moldynplot.dataset
from IPython import embed
import h5py
import numpy as np
import pandas as pd
from ..myplotspec.Dataset import Dataset
from ..myplotspec import sformat, wiprint
################################### CLASSES ###################################
class CorrDataset(Dataset):
    """
    Represents correlations between different datasets.
    """

    @classmethod
    def get_cache_key(cls, *args, **kwargs):
        """
        Generates tuple of arguments to be used as key for dataset
        cache.
        """
        import six
        from ..myplotspec import multi_get_copy

        x_kw = multi_get_copy(["x", "x_kw"], kwargs, {})
        x_cls = x_kw.get("cls", Dataset)
        if isinstance(x_cls, six.string_types):
            mod_name = ".".join(x_cls.split(".")[:-1])
            x_cls_name   = x_cls.split(".")[-1]
            mod = __import__(mod_name, fromlist=[x_cls_name])
            x_cls = getattr(mod, x_cls_name)

        y_kw = multi_get_copy(["y", "y_kw"], kwargs, {})
        y_cls = y_kw.get("cls", Dataset)
        if isinstance(y_cls, six.string_types):
            mod_name = ".".join(y_cls.split(".")[:-1])
            y_cls_name   = y_cls.split(".")[-1]
            mod = __import__(mod_name, fromlist=[y_cls_name])
            y_cls = getattr(mod, y_cls_name)

        return (cls, x_cls.get_cache_key(**x_kw), y_cls.get_cache_key(**y_kw))

    @staticmethod
    def get_cache_message(cache_key):
        """
        Generates message to be used when reloading previously-loaded
        dataset.

        Arguments:
          cache_key (tuple): key with which dataset object is stored in
            dataset cache

        Returns:
          str: message to be used when reloading previously-loaded
          dataset
        """
        return "Dataset previously loaded from '{0}' and '{1}'".format(
          cache_key[1][1], cache_key[2][1])

    def __init__(self, verbose=1, debug=0, **kwargs):
        """
        """
        from ..myplotspec import multi_get_copy

        # Load
        x_kw = multi_get_copy(["x", "x_kw"], kwargs, {})
        x_dataset = self.load_dataset(verbose=verbose, debug=debug, **x_kw)
        x_df = x_dataset.dataframe

        y_kw = multi_get_copy(["y", "y_kw"], kwargs, {})
        y_dataset = self.load_dataset(verbose=verbose, debug=debug, **y_kw)
        y_df = y_dataset.dataframe

        overlap_idx = x_df.index.intersection(y_df.index)
        corr_cols = [c for c in y_df.columns.values
                      if not c.endswith("_se")
                      and c in x_df.columns.values]

        overlap_cols = []
        for c in corr_cols:
            overlap_cols += [(c, "x"), (c, "y")]
            if c + "_se" in x_df.columns.values:
                overlap_cols += [(c + "_se", "x")]
            if c + "_se" in y_df.columns.values:
                overlap_cols += [(c + "_se", "y")]

        corr = pd.DataFrame(0, index=overlap_idx,
          columns=pd.MultiIndex.from_tuples(overlap_cols))
        corr.iloc[:, corr.columns.get_level_values(1)=="x"] = x_df[
          [a[0] for a in overlap_cols if a[1] == "x"]].loc[overlap_idx].values
        corr.iloc[:, corr.columns.get_level_values(1)=="y"] = y_df[
          [a[0] for a in overlap_cols if a[1] == "y"]].loc[overlap_idx].values

        self.dataframe = corr

#################################### MAIN #####################################
if __name__ == "__main__":
    import argparse

    # Prepare argument parser
    parser = argparse.ArgumentParser(
      description = """Processes datasets""")
    subparsers = parser.add_subparsers(
      dest        = "mode",
      description = "")

    CorrDataset.construct_argparser(subparsers)

    kwargs  = vars(parser.parse_args())
    kwargs.pop("cls")(**kwargs)
