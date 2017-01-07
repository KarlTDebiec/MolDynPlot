#!/usr/bin/python
# -*- coding: utf-8 -*-
#   moldynplot.dataset.TimeSeriesDataset.py
#
#   Copyright (C) 2015-2016 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Processes data that is a function of time

.. todo:
  - Fix separation and ordering of argument groups: input, action, output
  - Move in faster text loader from cpptraj2hdf5.py
  - Re-implement NatConTimeSeriesDataset
  - Re-implement SAXSTimeSeriesDataset
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
from .SequenceDataset import RelaxDataset, IREDDataset
from .SAXSDataset import SAXSDataset
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
        help_message = """Process standard data"""
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

        # Action arguments
        action_group = arg_groups.get("action",
          parser.add_argument_group("action"))
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
        Dataset.construct_argparser(parser)

        return parser

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
          calc_pdist (bool): Calculate probability distribution using
            :method:`calc_pdist`; store in instance variable `pdist_df`
          pdist_kw (dict): Keyword arguments used to configure
            probability distribution calculation
          calc_mean (bool): Calculate mean and standard error using
            :method:`calc_mean`; store in instance variable `sequence_df`
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
        self.dataset_cache = kwargs.get("dataset_cache", None)

        # Read data
        self.timeseries_df = self.read(**kwargs)

        # Process data
        if dt:
            self.timeseries_df.set_index(self.timeseries_df.index.values *
              float(dt), inplace=True)
            self.timeseries_df.index.name = "time"
        if toffset:
            index_name = self.timeseries_df.index.name
            self.timeseries_df.set_index(self.timeseries_df.index.values +
              float(toffset), inplace=True)
            self.timeseries_df.index.name = index_name
        if downsample:
            self.timeseries_df = self.downsample(df=self.timeseries_df, 
              downsample=downsample, **kwargs)

        # Output data
        if verbose >= 2:
            print("Processed timeseries DataFrame:")
            print(self.timeseries_df)
        if outfile is not None:
            self.write(df=self.timeseries_df, outfile=outfile, **kwargs)

        # Calculate probability distibution
        if calc_pdist:
            pdist_kw = kwargs.get("pdist_kw", {})
            self.pdist_df = self.calc_pdist(df=self.timeseries_df,
              verbose=verbose, **pdist_kw)
            if verbose >= 2:
                print("Processed pdist DataFrame:")
                print(self.pdist_df)
            # WRITE IF STRING

        # Calculate mean and standard error
        if calc_mean:
            block_kw = dict(min_n_blocks=2, max_cut=0.1, all_factors=False,
                            fit_exp=True, fit_sig=False)
            block_kw.update(kwargs.get("block_kw", {}))
            self.mean_df, self.block_averager = self.calc_mean(
              df=self.timeseries_df, verbose=verbose, **block_kw)
            if verbose >= 2:
                print("Processed mean DataFrame:")
                print(self.mean_df)
            # WRITE IF STRING

        # Interactive prompt
        if interactive:
            embed()

    @staticmethod
    def downsample(df, downsample, downsample_mode="mean", **kwargs):
        """
        Downsamples time series.

        Arguments:
          df (DataFrame): Timeseries DataFrame to downsample
          downsample (int): Interval by which to downsample points
          downsample_mode (str): Method of downsampling; may be 'mean'
            or 'mode'
          verbose (int): Level of verbose output
          kwargs (dict): Additional keyword arguments
        """
        from scipy.stats.mstats import mode

        # Process arguments
        verbose = kwargs.get("verbose", 1)

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

    @staticmethod
    def calc_mean(df, **kwargs):
        """
        Calculates the mean and standard error over a timeseries.

        Arguments:
          df (DataFrame): Timeseries DataFrame over which to calculate mean and
            standard error of each column over rows
          all_factors (bool): Use all factors by which the
            dataset is divisible rather than only factors of two
          min_n_blocks (int): Minimum number of blocks after
            transformation
          max_cut (float): Maximum proportion of dataset of
            omit in transformation
          fit_exp (bool): Fit exponential curve
          fit_sig (bool): Fit sigmoid curve
          verbose (int): Level of verbose output
          kwargs (dict): Additional keyword arguments

        Returns:
          DataFrame: DataFrame including mean and standard error for each
          column in *timeseries_df*
        """
        from ..fpblockaverager.FPBlockAverager import FPBlockAverager

        # Process arguments
        verbose = kwargs.get("verbose", 1)
        fit_exp = kwargs.get("fit_exp", True)
        fit_sig = kwargs.get("fit_sig", True)
        if verbose >= 1:
            wiprint("""Calculating mean and standard error over timeseries""")

        # Single-level columns
        if df.columns.nlevels == 1:
            mean_df = pd.DataFrame(data=df.mean(axis=0))
            block_averager = FPBlockAverager(df, **kwargs)
            if fit_exp and not fit_sig:
                errors = block_averager.parameters.loc[("exp", "a (se)")]
            elif fit_sig and not fit_exp:
                errors = block_averager.parameters.loc[("sig", "b (se)")]
            elif fit_exp and fit_sig:
                errors = block_averager.parameters.loc[("exp", "a (se)")]
            else:
                raise Exception()
            mean_df = mean_df.join(errors)
            mean_df.index.name = "q"
            mean_df.columns = ["intensity", "intensity se"]
        # Double-level columns
        elif df.columns.nlevels == 2:
            mean_df = pd.DataFrame(data=df.mean(axis=0))
            block_averager = FPBlockAverager(df, **kwargs)
            if fit_exp and not fit_sig:
                errors = block_averager.parameters.loc[("exp", "a (se)")]
            elif fit_sig and not fit_exp:
                errors = block_averager.parameters.loc[("sig", "b (se)")]
            elif fit_exp and fit_sig:
                errors = block_averager.parameters.loc[("exp", "a (se)")]
            else:
                raise Exception()
            errors.index = pd.MultiIndex.from_tuples(map(eval,
                             errors.index.values))
            errors.index = errors.index.set_levels([c+" se" for c in
                             errors.index.levels[1].values], level=1)
            mean_df = mean_df.squeeze().unstack().join(errors.unstack())
            mean_df = mean_df[["r1", "r1 se", "r2", "r2 se",
                               "noe", "noe se", "s2", "s2 se"]]
            mean_df = mean_df.loc[sorted(mean_df.index.values,
              key=lambda x: int(x.split(":")[1]))]
        else:
            raise Exception("Additional MultiIndex Levels not tested")

        return mean_df, block_averager


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
        help_message = """Process NMR relaxation data calculated from MD
          simulation using the iRED method as implemented in cpptraj"""
        if isinstance(parser_or_subparsers, argparse.ArgumentParser):
            parser = parser_or_subparsers
        elif isinstance(parser_or_subparsers, argparse._SubParsersAction):
            parser = parser_or_subparsers.add_parser(
              name        = "ired",
              description = help_message,
              help        = help_message)
        elif parser is None:
            parser = argparse.ArgumentParser(
              description = help_message)

        # Defaults
        if parser.get_default("cls") is None:
            parser.set_defaults(cls=IREDTimeSeriesDataset)

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
                if not (("r1" in columns and "r2" in columns
                and "noe" in columns) or ("s2" in columns)):
                    raise Exception(sformat("""DataFrame loaded from '{0}' does
                      not appear to contain either relaxation ('r1', 'r2',
                      'noe') or order parameter ('s2')
                      columns""".format(infile)))

        # Concatenate into timeseries
        df = self.concatenate_timeseries(timeseries_dfs, relax_dfs, order_dfs)
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

    def __init__(self, infile, address="saxs", dt=None, toffset=None,
        downsample=None, calc_mean=False, calc_error=True, error_method="std",
        scale=False, outfile=None, interactive=False, **kwargs):
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

        # Arguments
        verbose = kwargs.get("verbose", 1)
        self.dataset_cache = kwargs.get("dataset_cache", None)

        # Read data
        with h5py.File(expandvars(infile)) as h5_file:
            q = ["{0:5.3f}".format(a) for a in np.array(h5_file[address+"/q"])]
        self.timeseries_df = self.read(
          infile=infile+":/"+address+"/intensity",
          dataframe_kw=dict(columns=q), **kwargs)

        # Process data
        if dt:
            self.timeseries_df.set_index(self.timeseries_df.index.values *
              float(dt), inplace=True)
            self.timeseries_df.index.name = "time"
        if toffset:
            index_name = self.timeseries_df.index.name
            self.timeseries_df.set_index(self.timeseries_df.index.values +
              float(toffset), inplace=True)
            self.timeseries_df.index.name = index_name
        if downsample:
            self.downsample(downsample, downsample_mode="mean", **kwargs)

        # Output data
        if verbose >= 2:
            print("Processed timeseries DataFrame:")
            print(self.timeseries_df)
        if outfile is not None:
            self.write(df=self.timeseries_df, outfile=outfile, **kwargs)

        # Calculate mean and standard error
        if calc_mean:
            block_kw = dict(min_n_blocks=2, max_cut=0.1, all_factors=False,
                            fit_exp=True, fit_sig=False)
#            block_kw = dict(min_n_blocks=2, max_cut=0.1, all_factors=True,
#                            fit_exp=True, fit_sig=False)
            block_kw.update(kwargs.get("block_kw", {}))
            self.mean_df, self.block_averager = self.calc_mean(
              df=self.timeseries_df, verbose=verbose, **block_kw)
            if verbose >= 2:
                print("Processed mean DataFrame:")
                print(self.mean_df)
            self.df = self.mean_df
            # WRITE IF STRING

        # Scale
        if scale:
            self.scale(scale, **kwargs)

#################################### MAIN #####################################
def main():
    import argparse

    # Prepare argument parser
    parser = argparse.ArgumentParser(
      description = __doc__)
    subparsers = parser.add_subparsers(
      dest        = "mode",
      description = "")

    TimeSeriesDataset.construct_argparser(subparsers)
    IREDTimeSeriesDataset.construct_argparser(subparsers)
    PRETimeSeriesDataset.construct_argparser(subparsers)
#    NatConTimeSeriesDataset.construct_argparser(subparsers)
#    SAXSTimeSeriesDataset.construct_argparser(subparsers)

    kwargs  = vars(parser.parse_args())
    kwargs.pop("cls")(**kwargs)

if __name__ == "__main__":
    main()
