#!/usr/bin/python
# -*- coding: utf-8 -*-
#   moldynplot.Dataset.TimeSeriesDataset.py
#
#   Copyright (C) 2015-2016 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
.. todo:
  - Fix separation and ordering of argument groups: input, action, output
  - Move in faster text loader from cpptraj2hdf5.py
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
from .myplotspec.Dataset import Dataset
from .myplotspec import sformat, wiprint
################################### CLASSES ###################################
class TimeSeriesDataset(Dataset):
    """
    Represents data as a function of time.

    Attributes:
      timeseries_df (DataFrame): DataFrame whose index corresponds to
        time as represented by frame number or chemical time and whose
        columns are a series of quantities as a function of time.
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
        help_message = """Process data that is a function of time"""
        if isinstance(parser_or_subparsers, argparse.ArgumentParser):
            parser = parser_or_subparsers
        elif isinstance(parser_or_subparsers, argparse._SubParsersAction):
            parser = parser_or_subparsers.add_parser(
              name        = "timeseries",
              description = help_message,
              help        = help_message)
        elif parser is None:
            parser = argparse.ArgumentParser(
              description = help_message)

        # Defaults
        if parser.get_default("cls") is None:
            parser.set_defaults(cls=TimeSeriesDataset)

        # Arguments unique to this class
        arg_groups = {ag.title: ag for ag in parser._action_groups}
        action_group = arg_groups.get("action",
          parser.add_argument_group("action"))

        # Action arguments
        try:
            action_group.add_argument(
              "-dt",
              type     = float,
              help     = """time between frames""")
        except argparse.ArgumentError:
            pass
        try:
            action_group.add_argument(
              "-toffset",
              type     = float,
              help     = """offset to add to index (time or frame number)""")
        except argparse.ArgumentError:
            pass
        try:
            action_group.add_argument(
              "-downsample",
              type     = int,
              help     = """factor by which to downsample data""")
        except argparse.ArgumentError:
            pass
        try:
            action_group.add_argument(
              "--pdist",
              action   = "store_true",
              dest     = "calc_pdist",
              help     = """calculate probability distribution over timeseries
                         """)
        except argparse.ArgumentError:
            pass

        # Arguments inherited from superclass
        Dataset.construct_argparser(parser)

        return parser

    def __init__(self, dt=None, toffset=None, downsample=None,
        calc_pdist=False, outfile=None, interactive=False, **kwargs):
        """
        Arguments:
          infile{s} (list): Path(s) to input file(s); may contain
            environment variables and wildcards
          dt (float): Time interval between points; units unspecified
          toffset (float): Time offset to be added to all points (i.e.
            time of first point)
          downsample (int): Interval by which to downsample points
          downsample_mode (str): Method of downsampling; may be 'mean'
            or 'mode'
          calc_pdist (bool): Calculate probability distribution
          pdist_key (str): Column of which to calculate probability
            distribution
          kde_kw (dict): Keyword arguments passed to
            sklearn.neighbors.KernelDensity; key argument is 'bandwidth'
          grid (ndarray): Grid on which to calculate probability
            distribution
          interactive (bool): Provide iPython prompt and reading and
            processing data
          verbose (int): Level of verbose output
          kwargs (dict): Additional keyword arguments

          .. todo:
            - Calculate pdist using histogram
            - Verbose pdist
        """

        # Process arguments
        verbose = kwargs.get("verbose", 1)

        # Load
        self.timeseries_df = self.read(**kwargs)

        # Convert from frame index to time
        if dt:
            self.timeseries_df.set_index(self.timeseries_df.index.values *
              float(dt), inplace=True)
            self.timeseries_df.index.name = "time"

        # Offset time
        if toffset:
            index_name = self.timeseries_df.index.name
            self.timeseries_df.set_index(self.timeseries_df.index.values +
              float(toffset), inplace=True)
            self.timeseries_df.index.name = index_name

        # Downsample
        if downsample:
            self.timeseries_df = self.downsample(downsample, **kwargs)

        # Calculate probability distibution
        if calc_pdist:
            self.pdist_df = self.calc_pdist(**kwargs)

        # Output to screen
        if verbose >= 2:
            print("Processed timeseries DataFrame:")
            print(self.timeseries_df)
            if calc_pdist:
                print("Processed pdist DataFrame:")
                print(self.pdist_df)

        # Write data
        if outfile is not None:
            self.write(df=self.timeseries_df, outfile=outfile, **kwargs)

        # Interactive prompt
        if interactive:
            embed()

    def downsample(self, downsample, downsample_mode="mean", **kwargs):
        """
        Downsamples time series.

        Arguments:
          downsample (int): Interval by which to downsample points
          downsample_mode (str): Method of downsampling; may be 'mean'
            or 'mode'
          verbose (int): Level of verbose output
          kwargs (dict): Additional keyword arguments
        """
        from scipy.stats.mstats import mode

        # Process rguments
        verbose = kwargs.get("verbose", 1)
        df      = kwargs.get("df")
        if df is None:
            if hasattr(self, "timeseries_df"):
                df = self.timeseries_df
            else:
                raise()

        # Truncate dataset
        reduced = df.values[:df.shape[0] - (df.shape[0] % downsample),:]
        new_shape = (int(reduced.shape[0]/downsample), downsample,
          reduced.shape[1])
        index = np.reshape(df.index.values[
          :df.shape[0]-(df.shape[0] % downsample)],
          new_shape[:-1]).mean(axis=1)
        reduced = np.reshape(reduced, new_shape)

        # Downsample
        if downsample_mode == "mean":
            if verbose >= 1:
                wiprint("downsampling by factor of {0} using mean".format(
                  downsample))
            reduced = np.squeeze(reduced.mean(axis=1))
        elif downsample_mode == "mode":
            if verbose >= 1:
                wiprint("downsampling by factor of {0} using mode".format(
                  downsample))
            reduced = np.squeeze(mode(reduced, axis=1)[0])

        # Store downsampled time series
        reduced = pd.DataFrame(data=reduced, index=index,
          columns=df.columns.values)
        reduced.index.name = "time"
        df = reduced

        return df

#    def calc_error(self, error_method="std", **kwargs):
#        """
#        Calculates standard error using time series data.
#
#        .. todo:
#          - Support breaking into blocks (essentially downsampling,
#            then calculating standard error)
#        """
#
#        # Arguments
#        verbose = kwargs.get("verbose", 1)
#        timeseries=self.timeseries
#
#        # Calculate standard error
#        if error_method == "std":
#            if verbose >= 1:
#                print("calculating standard error using standard deviation")
#            se = timeseries.std()
#        elif error_method == "block":
#            from .fpblockaverager.FPBlockAverager import FPBlockAverager
#            if verbose >= 1:
#                print("calculating standard error using block averaging")
#            ba = FPBlockAverager(timeseries, **kwargs)
#            se = ba.parameters.loc["exp", "a (se)"]
#        else:
#            if verbose >= 1:
#                print("error_method '{0}' not understood, ".format(scale) +
#                      "must be one of 'std', 'block'; not calculating error.")
#            return
#
#        return se

    def calc_pdist(self, **kwargs):
        """
        Calcualtes probability distribution of time series.

        Arguments:
          pdist_kw (dict): Keyword arguments used to configure
            probability distribution calculation
          verbose (int): Level of verbose output
          kwargs (dict): Additional keyword arguments
        """
        from collections import OrderedDict
        from sklearn.neighbors import KernelDensity

        # Process arguments
        verbose = kwargs.get("verbose", 1)
        if verbose >= 1:
            wiprint("""Calculating probability distribution over timeseries""")
        df      = kwargs.get("df")
        if df is None:
            if hasattr(self, "timeseries_df"):
                df = self.timeseries_df
            else:
                raise()

        pdist_kw = kwargs.get("pdist_kw", {"bandwidth": 0.1})
        mode = "kde"
        if mode == "kde":

            # Prepare bandwidths
            bandwidth = pdist_kw.pop("bandwidth", None)
            if bandwidth is None:
                all_bandwidth = None
                bandwidth = {}
            elif isinstance(bandwidth, float):
                all_bandwidth = bandwidth
                bandwidth = {}
            elif isinstance(bandwidth, dict):
                all_bandwidth = None
                pass
            for column, series in df.iteritems():
                if column in bandwidth:
                    bandwidth[column] = float(bandwidth[column])
                elif all_bandwidth is not None:
                    bandwidth[column] = all_bandwidth
                else:
                    bandwidth[column] = series.std() / 10.0

            # Prepare grids
            grid = pdist_kw.pop("grid", None)
            if grid is None:
                all_grid = None
                grid = {}
            elif isinstance(grid, list) or isinstance(grid, np.ndarray):
                all_grid = np.array(grid)
                grid = {}
            elif isinstance(grid, dict):
                all_grid = None
                pass
            for column, series in df.iteritems():
                if column in grid:
                    grid[column] = np.array(grid[column])
                elif all_grid is not None:
                    grid[column] = all_grid
                else:
                    grid[column] = np.linspace(series.min() - series.std(),
                                      series.max() + series.std(), 100)

            # Calculate probability distributions
            kde_kw = pdist_kw.get("kde_kw", {})
            pdist = OrderedDict()
            for column, series in df.iteritems():
                if verbose >= 1:
                   wiprint("calculating probability distribution of "
                    "{0} using a kernel density estimate".format(column))
                kde = KernelDensity(bandwidth=bandwidth[column], **kde_kw)
                kde.fit(series[:, np.newaxis])
                pdf = np.exp(kde.score_samples(grid[column][:, np.newaxis]))
                pdf /= pdf.sum()
                series_pdist = pd.DataFrame(pdf, index=grid[column],
                  columns=["probability"])
                series_pdist.index.name = column
                pdist[column] = series_pdist

            return pdist

    def timeseries_to_sequence(self, **kwargs):
        """
        Calculates the mean and standard error over a timeseries.

        Arguments:
          {timeseries_}d{ata}f{rame} (DataFrame): Timeseries over which
            to calculate mean and standard error; if omitted looks for
            :attr:`timeseries_df`
          block_kw (dict): Keyword arguments passed to
            :class:`fpblockaverager.FPBlockAverager`
          block_kw[all_factors] (bool): Use all factors by which the
          dataset is divisible rather than only factors of two
          block_kw[min_n_blocks] (int): Minimum number of blocks after
            transformation
          block_kw[max_cut] (float): Maximum proportion of dataset of
            omit in transformation
          block_kw[fit_exp] (bool): Fit exponential curve
          block_kw[fit_sig] (bool): Fit sigmoid curve
          verbose (int): Level of verbose output
          kwargs (dict): Additional keyword arguments

        Returns:
          DataFrame: Sequence dataframe including mean and standard
          error for each column in *timeseries_df*

        .. todo:
          - Make general, not specific to sequence data
        """
        from .myplotspec import multi_get

        # Process arguments
        verbose = kwargs.get("verbose", 1)
        timeseries_df = multi_get(["timeseries_dataframe", "timeseries_df",
              "dataframe", "df"], kwargs)
        if timeseries_df is None:
            if hasattr(self, "timeseries_df"):
                timeseries_df = self.timeseries_df
            else:
                raise Exception("Cannot find timeseries DataFrame")
        if verbose >= 1:
            wiprint("""Calculating mean and standard error over timeseries""")

        sequence_df = pd.DataFrame(data=timeseries_df.mean(axis=0))

        from .fpblockaverager.FPBlockAverager import FPBlockAverager
        block_kw = dict(min_n_blocks=2, max_cut=0.1, all_factors=True,
                        fit_exp=True, fit_sig=False)
        block_kw.update(kwargs.get("block_kw", {}))
        block_averager = FPBlockAverager(timeseries_df, **block_kw)
        errors       = block_averager.parameters.loc[("exp", "a (se)")]
        errors.index = pd.MultiIndex.from_tuples(map(eval, errors.index.values))
        errors.index = errors.index.set_levels([c+" se" for c in
                        errors.index.levels[1].values], level=1)

        a = sequence_df.squeeze().unstack()
        b = errors.unstack()
        c = a.join(b)
        c = c[["r1", "r1 se", "r2", "r2 se", "noe", "noe se", "s2", "s2 se"]]
        self.sequence_df = c
        c = c.loc[sorted(c.index.values, key=lambda x: int(x.split(":")[1]))]
        sequence_df = c

        self.block_averager = block_averager
        return sequence_df

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
        help_message = """Process simulated paramagnetic relaxation enhancement
          data from MD simulation """
        if isinstance(parser_or_subparsers, argparse.ArgumentParser):
            parser = parser_or_subparsers
        elif isinstance(parser_or_subparsers, argparse._SubParsersAction):
            parser = parser_or_subparsers.add_parser(
              name        = "pre_timeseries",
              description = help_message,
              help        = help_message)
        elif parser is None:
            parser = argparse.ArgumentParser(
              description = help_message)

        # Defaults
        if parser.get_default("cls") is None:
            parser.set_defaults(cls=PRETimeSeriesDataset)

        # Arguments unique to this class

        # Arguments inherited from superclass
        TimeSeriesDataset.construct_argparser(parser)

        return parser

    def __init__(self, dt=None, toffset=None, downsample=None,
        calc_pdist=False, outfile=None, interactive=False, **kwargs):
        """
        Arguments:
          infile{s} (list): Path(s) to input file(s); may contain
            environment variables and wildcards
          dt (float): Time interval between points; units unspecified
          toffset (float): Time offset to be added to all points (i.e.
            time of first point)
          downsample (int): Interval by which to downsample points
          downsample_mode (str): Method of downsampling; may be 'mean'
            or 'mode'
          calc_pdist (bool): Calculate probability distribution
          pdist_key (str): Column of which to calculate probability
            distribution
          kde_kw (dict): Keyword arguments passed to
            sklearn.neighbors.KernelDensity; key argument is 'bandwidth'
          grid (ndarray): Grid on which to calculate probability
            distribution
          interactive (bool): Provide iPython prompt and reading and
            processing data
          verbose (int): Level of verbose output
          kwargs (dict): Additional keyword arguments
        """

        # Process arguments
        verbose = kwargs.get("verbose", 1)

        # Load
        self.timeseries_df = self.read(**kwargs)

        # Convert from frame index to time
        if dt:
            self.timeseries_df.set_index(self.timeseries_df.index.values *
              float(dt), inplace=True)
            self.timeseries_df.index.name = "time"

        # Offset time
        if toffset:
            index_name = self.timeseries_df.index.name
            self.timeseries_df.set_index(self.timeseries_df.index.values +
              float(toffset), inplace=True)
            self.timeseries_df.index.name = index_name

        # Downsample
        if downsample:
            self.timeseries_df = self.downsample(downsample, **kwargs)

        # Calculate probability distibution
        if calc_pdist:
            self.pdist_df = self.calc_pdist(**kwargs)

        # Output to screen
        if verbose >= 2:
            print("Processed timeseries DataFrame:")
            print(self.timeseries_df)
            if calc_pdist:
                print("Processed pdist DataFrame:")
                print(self.pdist_df)

        # Write data
        if outfile is not None:
            self.write(df=self.timeseries_df, outfile=outfile, **kwargs)

        # Interactive prompt
        if interactive:
            embed()
    
class IREDTimeSeriesDataset(TimeSeriesDataset, IREDDataset):
    """
    Represents iRED NMR relaxation data as a function of time and
    residue number.
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
        help_message = """Process relaxation data calculated from MD
          simulation using the iRED method as implemented in cpptraj;
          treat infiles as consecutive (potentially overlapping)
          excerpts of a longer simulation; processed results are a
          timeseries and the average across the timeseries, including
          standard errors calculated using block averaging"""
        if isinstance(parser_or_subparsers, argparse.ArgumentParser):
            parser = parser_or_subparsers
        elif isinstance(parser_or_subparsers, argparse._SubParsersAction):
            parser = parser_or_subparsers.add_parser(
              name        = "ired_timeseries",
              description = help_message,
              help        = help_message)
        elif parser is None:
            parser = argparse.ArgumentParser(
              description = help_message)

        # Defaults
        if parser.get_default("cls") is None:
            parser.set_defaults(cls=IREDTimeSeriesDataset)

        # Arguments unique to this class
        arg_groups = {ag.title: ag for ag in parser._action_groups}

        # Action arguments
        action_group = arg_groups.get("action",
          parser.add_argument_group("action"))
        try:
            action_group.add_argument(
              "--mean",
              action   = "store_true",
              dest     = "calc_mean",
              help     = """Calculate mean and standard error over timeseries
                         """)
        except argparse.ArgumentError:
            pass

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
          relax_dfs (list): DataFrames containing data from relax infiles
          order_dfs (list): DataFrames containing data from order infiles
          verbose (int): Level of verbose output
          kwargs (dict): Additional keyword arguments

        Returns:
          df (DataFrame): Averaged DataFrame including relax and order
        """

        # Process arguments
        verbose = kwargs.get("verbose", 1)

        # Process timeseris
        if len(timeseries_dfs) >= 1:
            if verbose >= 1:
                wiprint("""Concatenating timeseries from {0} timeseries
                        infiles""".format(len(timeseries_dfs)))
            df = pd.concat(timeseries_dfs)
        else:
            df = pd.DataFrame()

        # Process relaxation
        if len(relax_dfs) >= 1:
            if verbose >= 1:
                wiprint("""Concatenating timeseries from {0} relaxation
                        infiles""".format(len(relax_dfs)))
            relax_df = pd.concat([df.stack() for df in relax_dfs],
                         axis=1).transpose()
        else:
            relax_df = None

        # Process order parameters
        if len(order_dfs) >= 1:
            if verbose >= 1:
                wiprint("""Concatenating timeseries from {0} order parameter
                        infiles""".format(len(order_dfs)))
            order_df = pd.concat([df.stack() for df in order_dfs],
                         axis=1).transpose()
        else:
            order_df = None

        # Merge and sort
        if relax_df is not None and order_df is not None:
            df = pd.merge(relax_df, order_df, how="outer", left_index=True,
                           right_index=True)
            df = df[sorted(list(set(df.columns.get_level_values(0))),
                   key=lambda x: int(x.split(":")[1]))]
        elif relax_df is None and order_df is not None:
            df = order_df
        elif order_df is None and relax_df is not None:
            df = relax_df

        return df

    def __init__(self, dt=None, toffset=None, downsample=None,
        calc_pdist=False, calc_mean=False, outfile=None, interactive=False,
        **kwargs):
        """
        Arguments:
          infile{s} (list): Path(s) to input file(s); may contain
            environment variables and wildcards
          dt (float): Time interval between points; units unspecified
          toffset (float): Time offset to be added to all points (i.e.
            time of first point)
          downsample (int): Interval by which to downsample points
          downsample_mode (str): Method of downsampling; may be 'mean'
            or 'mode'
          calc_pdist (bool): Calculate probability distribution
          interactive (bool): Provide iPython prompt and reading and
            processing data
          verbose (int): Level of verbose output
          kwargs (dict): Additional keyword arguments
        """
        # Process arguments
        verbose = kwargs.get("verbose", 1)
        self.dataset_cache = kwargs.get("dataset_cache", None)

        # Load
        self.timeseries_df = self.read(**kwargs)

        # Convert from frame index to time
        if dt:
            self.timeseries_df.set_index(self.timeseries_df.index.values *
              float(dt), inplace=True)
            self.timeseries_df.index.name = "time"

        # Offset time
        if toffset:
            index_name = self.timeseries_df.index.name
            self.timeseries_df.set_index(self.timeseries_df.index.values +
              float(toffset), inplace=True)
            self.timeseries_df.index.name = index_name

        # Downsample
        if downsample:
            self.timeseries_df = self.downsample(downsample, **kwargs)

        # Calculate probability distibution
        if calc_pdist:
            self.pdist_df = self.calc_pdist(**kwargs)

        # Output to screen
        if verbose >= 2:
            print("Processed timeseries DataFrame:")
            print(self.timeseries_df)

        # Prepare sequence dataframe using averages and block standard errors
        if calc_mean:
            self.sequence_df = self.timeseries_to_sequence(
                                 df=self.timeseries_df, **kwargs)
            # Output to screen
            if verbose >= 2:
                print("Processed sequence DataFrame:")
                print(self.sequence_df)

        # Write data
        if outfile is not None:
            self.write(df=self.timeseries_df, outfile=outfile, **kwargs)

        # Interactive prompt
        if interactive:
            embed()

    def read(self, **kwargs):
        """
        Reads iRED time series data from one or more *infiles* into a
        DataFrame.
        """
        import re
        from .myplotspec import multi_pop_merged

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
                if not (("r1" in columns and "r2" in columns
                and "noe" in columns) or ("s2" in columns)):
                    raise Exception(sformat("""DataFrame loaded from '{0}' does
                      not appear to contain either relaxation ('r1', 'r2',
                      'noe') or order parameter ('s2')
                      columns""".format(infile)))

        # Concatenate into timeseries
        df = self.concatenate_timeseries(timeseries_dfs, relax_dfs, order_dfs)
        df.index.name = "frame"
        return df

class NatConTimeSeriesDataset(TimeSeriesDataset):
    """
    Manages native contact datasets.
    """

    def __init__(self, downsample=None, calc_pdist=True, **kwargs):
        """
        Arguments:
          infile (str): Path to input file, may contain environment
            variables
          usecols (list): Columns to select from DataFrame, once
            dataframe has already been loaded
          dt (float): Time interval between points; units unspecified
          toffset (float): Time offset to be added to all points (i.e.
            time of first point)
          cutoff (float): Minimum distance within which a contact is
            considered to be formed
          downsample (int): Interval by which to downsample points using
            mode
          pdist (bool): Calculate probability distribution
          verbose (int): Level of verbose output
          kwargs (dict): Additional keyword arguments
        """
        verbose = kwargs.get("verbose", 1)

        # Load
        super(NatConTimeSeriesDataset, self).__init__(**kwargs)
        dataframe = self.dataframe
        n_contacts = self.dataframe.shape[1]

        # Convert minimum distances to percent native contacts
        cutoff = kwargs.get("cutoff", 5.5)
        percent = pd.DataFrame(data=(dataframe.values <= cutoff).sum(axis=1)
          / dataframe.shape[1], index=dataframe.index,
          columns=["percent_native_contacts"])
        dataframe = self.dataframe = percent

        # Downsample; flag included in function definition to prevent
        #   superclass from downsampling before applying cutoff
        if downsample is not None:
            from scipy.stats.mstats import mode

            if verbose >= 1:
                print("downsampling by factor of {0} using mode".format(
                  downsample))

            reduced = dataframe.values[
              :dataframe.shape[0]-(dataframe.shape[0] % downsample),:]
            new_shape=(int(reduced.shape[0]/ downsample),
                downsample, reduced.shape[1])
            index = np.reshape(dataframe.index.values[
              :dataframe.shape[0]-(dataframe.shape[0] % downsample)],
              new_shape[:-1]).mean(axis=1)
            reduced = np.reshape(reduced, new_shape)
            reduced = np.squeeze(mode(reduced, axis=1)[0])
            reduced = pd.DataFrame(data=reduced, index=index,
              columns=dataframe.columns.values)
            reduced.index.name = "time"
            dataframe = self.dataframe = reduced

        # Calculate probability distribution
        if calc_pdist:
            if verbose >= 1:
                print("calculating probability distribution using histogram")
            bins = np.linspace(0-((1/n_contacts)/2), 1+((1/n_contacts)/2),
              n_contacts+1)
            pdist, _ = np.histogram(self.dataframe.values, bins)
            pdist    =  np.array(pdist, np.float) / pdist.sum()
            pdist_x = np.zeros(bins.size*2)
            pdist_y = np.zeros(bins.size*2)
            pdist_x[::2]    = pdist_x[1::2]   = bins
            pdist_y[1:-1:2] = pdist_y[2:-1:2] = pdist
            self.pdist_x = pdist_x
            self.pdist_y = pdist_y

        self.timeseries = dataframe

class SAXSTimeSeriesDataset(TimeSeriesDataset, SAXSDataset):
    """
    Manages Small Angle X-ray Scattering time series datasets.
    """

    def __init__(self, infile, address="saxs", downsample=None,
        calc_mean=False, calc_error=True, error_method="std", scale=False,
        **kwargs):
        """
        Arguments:
          infile (str): Path to input file, may contain environment
            variables
          usecols (list): Columns to select from DataFrame, once
            dataframe has already been loaded
          dt (float): Time interval between points; units unspecified
          toffset (float): Time offset to be added to all points (i.e.
            time of first point)
          downsample (int): Interval by which to downsample points
          verbose (int): Level of verbose output
          kwargs (dict): Additional keyword arguments

        .. todo:
          - Calculate error
          - Shift downsampling to superclass
        """
        from os.path import expandvars

#        # Arguments
#        verbose = kwargs.get("verbose", 1)
#
#        # Load
#        with h5py.File(expandvars(infile)) as h5_file:
#            q = ["{0:5.3f}".format(a) for a in np.array(h5_file[address+"/q"])]
#        super(SAXSTimeSeriesDataset, self).__init__(infile=infile,
#          address=address+"/intensity", dataframe_kw=dict(columns=q), **kwargs)
#        timeseries = self.timeseries = self.dataframe
#
#        # Downsample
#        if downsample:
#            self.downsample(downsample, downsample_mode="mean", **kwargs)
#
#        # Average over time series
#        if calc_mean:
#            self.dataframe = dataframe = pd.DataFrame(
#              data=timeseries.mean(axis=0), columns=["intensity"])
#            dataframe.index.name = "q"
#            dataframe.index = np.array(timeseries.columns.values, np.float)
#
#            # Scale
#            if scale:
##                curve_fit_kw = dict(p0=(2e-9), bounds=(0.0,0.35))
#                curve_fit_kw = dict(p0=(2e-9))  # Not clear why bounds broke
#                curve_fit_kw.update(kwargs.get("curve_fit_kw", {}))
#                scale = self.scale(scale, curve_fit_kw=curve_fit_kw, **kwargs)
#                self.timeseries *= scale
#        elif scale:
#            self.timeseries *= scale
#        if calc_error:
#            se = self.calc_error(error_method="block", **kwargs)
#            se.name = "intensity_se"
#            dataframe = self.dataframe = pd.concat([dataframe, se], axis=1)

#################################### MAIN #####################################
if __name__ == "__main__":
    import argparse

    # Prepare argument parser
    parser = argparse.ArgumentParser(
      description = """Processes datasets""")
    subparsers = parser.add_subparsers(
      dest        = "mode",
      description = "")

    TimeSeriesDataset.construct_argparser(subparsers)
    IREDTimeSeriesDataset.construct_argparser(subparsers)
    PRETimeSeriesDataset.construct_argparser(subparsers)

    kwargs  = vars(parser.parse_args())
    kwargs.pop("cls")(**kwargs)
