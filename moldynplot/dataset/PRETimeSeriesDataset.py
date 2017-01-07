#!/usr/bin/python
# -*- coding: utf-8 -*-
#   moldynplot.dataset.PRETimeSeriesDataset.py
#
#   Copyright (C) 2015-2016 Karl T Debiec
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
from IPython import embed
import h5py
import numpy as np
import pandas as pd
from ..myplotspec.Dataset import Dataset
from ..myplotspec import sformat, wiprint
from .SequenceDataset import RelaxDataset
from .TimeSeriesDataset import TimeSeriesDataset


################################### CLASSES ###################################
class PRETimeSeriesDataset(TimeSeriesDataset, RelaxDataset):
    """
    Represents paramagnetic relaxation enhancement data as a function of time
    and residue number.
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
        if parser.get_default("cls") is None:
            parser.set_defaults(cls=PRETimeSeriesDataset)

        # Arguments unique to this class

        # Arguments inherited from superclass
        TimeSeriesDataset.construct_argparser(parser)

        return parser

    def __init__(self, **kwargs):
        """

        Args:
            **kwargs (dict):
        """
        super(PRETimeSeriesDataset, self).__init__(**kwargs)


#################################### MAIN #####################################
def main():
    import argparse

    # Prepare argument parser
    parser = argparse.ArgumentParser(description=__doc__)

    PRETimeSeriesDataset.construct_argparser(parser)

    kwargs = vars(parser.parse_args())
    kwargs.pop("cls")(**kwargs)


if __name__ == "__main__":
    main()
