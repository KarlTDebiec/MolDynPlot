#!/usr/bin/python
# -*- coding: utf-8 -*-
#   moldynplot.dataset.IREDTimeSeriesDataset.py
#
#   Copyright (C) 2015-2017 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Represents iRED NMR relaxation data as a function of time and residue number
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
import six
from .IREDDataset import IREDDataset
from .TimeSeriesDataset import TimeSeriesDataset
from ..myplotspec.Dataset import Dataset
from ..myplotspec import wiprint


################################### CLASSES ###################################
class IREDTimeSeriesDataset(TimeSeriesDataset, IREDDataset):
    """
    Represents iRED NMR relaxation data as a function of time and
    residue number
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
        help_message = """Process NMR relaxation data calculated from MD
          simulation using the iRED method as implemented in cpptraj"""
        if isinstance(parser_or_subparsers, argparse.ArgumentParser):
            parser = parser_or_subparsers
        elif isinstance(parser_or_subparsers, argparse._SubParsersAction):
            parser = parser_or_subparsers.add_parser(name="ired",
              description=help_message, help=help_message)
        elif parser is None:
            parser = argparse.ArgumentParser(description=help_message)

        # Defaults
        if parser.get_default("class_") is None:
            parser.set_defaults(class_=IREDTimeSeriesDataset)

        # Arguments unique to this class

        # Arguments inherited from superclass
        IREDDataset.construct_argparser(parser)
        TimeSeriesDataset.construct_argparser(parser)

        return parser

    @staticmethod
    def concatenate_timeseries(timeseries_dfs=None, relax_dfs=None,
      order_dfs=None, **kwargs):
        """
        Concatenates a series of iRED datasets.

        Arguments:
          timeseries_dfs (list): DataFrames containing data from
            timeseries infiles
          relax_dfs (list): DataFrames containing data from relax infiles
          order_dfs (list): DataFrames containing data from order infiles
          verbose (int): Level of verbose output
          kwargs (dict): Additional keyword arguments

        Returns:
          DataFrame: Concatenated DataFrame including relax and order
        """

        # Process arguments
        verbose = kwargs.get("verbose", 1)

        # Process timeseries
        if len(timeseries_dfs) >= 1:
            if verbose >= 1:
                wiprint("""Concatenating timeseries from {0} timeseries
                        infiles""".format(len(timeseries_dfs)))
            timeseries_df = pd.concat(timeseries_dfs)
        else:
            timeseries_df = None

        # Process relaxation
        if len(relax_dfs) >= 1:
            if verbose >= 1:
                wiprint("""Concatenating timeseries from {0} relaxation
                        infiles""".format(len(relax_dfs)))
            relax_df = pd.concat([rdf.stack() for rdf in relax_dfs],
              axis=1).transpose()
        else:
            relax_df = None

        # Process order parameters
        if len(order_dfs) >= 1:
            if verbose >= 1:
                wiprint("""Concatenating timeseries from {0} order parameter
                        infiles""".format(len(order_dfs)))
            order_df = pd.concat([odf.stack() for odf in order_dfs],
              axis=1).transpose()
        else:
            order_df = None

        # Merge and sort relaxation and order parameters
        if relax_df is not None and order_df is not None:
            df = pd.merge(relax_df, order_df, how="outer", left_index=True,
              right_index=True)
            df = df[sorted(list(set(df.columns.get_level_values(0))),
              key=lambda x: int(x.split(":")[1]))]
        elif relax_df is None and order_df is not None:
            df = order_df
        elif order_df is None and relax_df is not None:
            df = relax_df
        else:
            df = None

        # Append to existing timeseries_df
        if timeseries_df is not None and df is not None:
            df = pd.concat([timeseries_df, df])
        elif df is None:
            df = timeseries_df

        return df

    def read(self, **kwargs):
        """
        Reads iRED time series data from one or more *infiles* into a
        DataFrame.
        """
        import re
        from ..myplotspec import multi_pop_merged

        # Process arguments
        infile_args = multi_pop_merged(["infile", "infiles"], kwargs)
        infiles = self.infiles = self.process_infiles(infiles=infile_args)
        if len(infiles) == 0:
            raise Exception(sformat("""No infiles found matching
            '{0}'""".format(infile_args)))
        re_h5 = re.compile(
          r"^(?P<path>(.+)\.(h5|hdf5))((:)?(/)?(?P<address>.+))?$",
          flags=re.UNICODE)

        # Load data
        timeseries_dfs = []
        relax_dfs = []
        order_dfs = []
        for infile in infiles:
            if re_h5.match(infile):
                df = self._read_hdf5(infile, **kwargs)
            else:
                df = self._read_text(infile, **kwargs)
            if df.columns.nlevels == 2:
                timeseries_dfs.append(df)
            else:
                columns = df.columns.values

                if "r1" in columns and "r2" in columns and "noe" in columns:
                    relax_dfs.append(df)
                if "s2" in columns:
                    order_dfs.append(df)
                if not ((
                              "r1" in columns and "r2" in columns and "noe"
                  in columns) or (
                      "s2" in columns)):
                    raise Exception(sformat("""DataFrame loaded from '{0}' does
                      not appear to contain either relaxation ('r1', 'r2',
                      'noe') or order parameter ('s2')
                      columns""".format(infile)))

        # Concatenate into timeseries
        df = self.concatenate_timeseries(timeseries_dfs, relax_dfs, order_dfs)
        return df


#################################### MAIN #####################################
if __name__ == "__main__":
    IREDTimeSeriesDataset.main()