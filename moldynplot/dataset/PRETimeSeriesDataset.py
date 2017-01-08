#!/usr/bin/python
# -*- coding: utf-8 -*-
#   moldynplot.dataset.PRETimeSeriesDataset.py
#
#   Copyright (C) 2015-2017 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Represents paramagnetic relaxation enhancement timeseries data
"""
################################### MODULES ###################################
from __future__ import (absolute_import, division, print_function,
    unicode_literals)

if __name__ == "__main__":
    __package__ = str("moldynplot.dataset")
    import moldynplot.dataset
from .RelaxDataset import RelaxDataset
from .TimeSeriesDataset import TimeSeriesDataset


################################### CLASSES ###################################
class PRETimeSeriesDataset(TimeSeriesDataset, RelaxDataset):
    """
    Represents paramagnetic relaxation enhancement timeseries data
    """

    @staticmethod
    def construct_argparser(parser_or_subparsers=None, **kwargs):
        """
        Adds arguments to an existing argument parser, constructs a
        subparser, or constructs a new parser

        Arguments:
          parser_or_subparsers (ArgumentParser, _SubParsersAction,
            optional): If ArgumentParser, existing parser to which
            arguments will be added; if _SubParsersAction, collection of
            subparsers to which a new argument parser will be added; if
            None, a new argument parser will be generated
          kwargs (dict): Additional keyword arguments

        Returns:
          ArgumentParser: Argument parser or subparser
        """
        import argparse

        # Process arguments
        help_message = """Process NMR paramagnetic relaxation enhancement
                       data"""
        if isinstance(parser_or_subparsers, argparse.ArgumentParser):
            parser = parser_or_subparsers
        elif isinstance(parser_or_subparsers, argparse._SubParsersAction):
            parser = parser_or_subparsers.add_parser(name="pre",
              description=help_message, help=help_message)
        elif parser is None:
            parser = argparse.ArgumentParser(description=help_message)

        # Defaults
        if parser.get_default("class_") is None:
            parser.set_defaults(class_=PRETimeSeriesDataset)

        # Arguments unique to this class

        # Arguments inherited from superclass
        TimeSeriesDataset.construct_argparser(parser)

        return parser

    def __init__(self, dt=None, **kwargs):
        """

        Args:
            verbose (int): Level of verbose output
            kwargs (dict): Additional Keyword Arguments
        """

        # Process arguments
        verbose = kwargs.get("verbose", 1)
        self.dataset_cache = kwargs.get("dataset_cache", None)

        # Read data
        if not hasattr(self, "timeseries_df"):
            self.timeseries_df = self.df = self.read(**kwargs)
        if dt:
            self.timeseries_df.set_index(
              self.timeseries_df.index.values * float(dt), inplace=True)
            self.timeseries_df.index.name = "time"

        import pandas as pd
        import numpy as np
        pd.set_option('display.max_rows', 500)
        # print(self.timeseries_df)
        mean_df = pd.DataFrame(data=self.timeseries_df.mean(axis=0),
          dtype=np.double)
        print(mean_df)

        k = 0.0123  # Å^6 ns-2
        tc = 4  # ns
        v = 8e-7  # ns-1
        mean_df = 1 / (mean_df ** 6)  # Å-6
        print(mean_df)
        mean_df = k * mean_df  # ns-2
        print(mean_df)
        mean_df = mean_df * 1e9 * 1e9  # s-2
        print(mean_df)
        const = (4 * 4e-9 + ((3 * 4e-9) / (1 + (800 * (4e-9 ** 2)))))  # s
        print(const)
        mean_df = mean_df * const  # s-1
        print(mean_df)
        mean_df = 1 / mean_df  # s


        # embed()
        # Allow superclass to handle other actions
        # super(PRETimeSeriesDataset, self).__init__(**kwargs)


#################################### MAIN #####################################
if __name__ == "__main__":
    PRETimeSeriesDataset.main()