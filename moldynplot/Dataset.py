#!/usr/bin/python
# -*- coding: utf-8 -*-
#   moldynplot.Dataset.py
#
#   Copyright (C) 2015-2016 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Moldynplot includes several dataset classes that build on
:class:`Dataset<myplotspec.Dataset.Dataset>` with additions
specific for molecular dynamics simulation data.

.. todo:
  - FIX SEPARATION AND ORDERING OF ARGUMENT GROUPS: input, action, output
  - Move relaxation error here
  - Move relaxation heteronuclear noe here
  - Move relaxation pre ratio here
  - Move SAXS here
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
class SequenceDataset(Dataset):
    """
    Represents data that is a function of amino acid sequence.

    Attributes:
      sequence_df (DataFrame): DataFrame whose index corresponds to
        amino acid residue number in the form ``XAA:#`` and whose
        columns are a series of quantities specific to each residue.
        Standard errors of these quantities may be represented by
        adjacent columns with ' se' appended to their names. ::

                         r1     r1 se        r2     r2 se  ...
          residue
          GLN:2    2.451434  0.003734  5.041334  0.024776  ...
          TYR:3    2.443613  0.004040  5.138383  0.025376  ...
          LYS:4    2.511626  0.004341  5.589428  0.026236  ...
          ...      ...       ...       ...       ...       ...
    """

    default_hdf5_address = "/"
    default_hdf5_kw = dict(
      chunks      = True,
      compression = "gzip",
      dtype       = np.float32,
      scaleoffset = 5)

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
        help_message = """Process data that is a function of amino acid
          sequence"""
        if isinstance(parser_or_subparsers, argparse.ArgumentParser):
            parser = parser_or_subparsers
        elif isinstance(parser_or_subparsers, argparse._SubParsersAction):
            parser = parser_or_subparsers.add_parser(
              name        = "sequence",
              description = help_message,
              help        = help_message)
        elif parser is None:
            parser = argparse.ArgumentParser(
              description = help_message)

        # Defaults
        if parser.get_default("cls") is None:
            parser.set_defaults(cls=SequenceDataset)

        # Arguments unique to this class
        arg_groups = {ag.title: ag for ag in parser._action_groups}

        # Input arguments
        input_group  = arg_groups.get("input",
          parser.add_argument_group("input"))
        try:
            input_group.add_argument(
              "-indexfile",
              required = False,
              type     = str,
              help     = """text file from which to load residue names; should
                         list amino acids in the form 'XAA:#' separated by
                         whitespace; if omitted will be taken from rows of
                         first infile; may contain environment variables""")
        except argparse.ArgumentError:
            pass

        # Arguments inherited from superclass
        Dataset.construct_argparser(parser)

        return parser

    @classmethod
    def get_cache_key(cls, **kwargs):
        """
        Generates key for dataset cache.

        See :class:`SequenceDataset<moldynplot.Dataset.SequenceDataset>`
        for argument details.

        Returns:
          tuple: Cache key; contains arguments sufficient to reconstruct
          dataset
        """
        from .myplotspec import multi_pop_merged

        # Process arguments
        infiles = multi_pop_merged(["infile", "infiles"], kwargs)
        infiles = SequenceDataset.process_infiles(infiles=infiles)
        if infiles is None:
            return None
        read_csv_kw = []
        if "use_indexes" in kwargs:
            use_indexes = tuple(kwargs.get("use_indexes"))
        else:
            use_indexes = None
        for key, value in kwargs.get("read_csv_kw", {}).items():
            if isinstance(value, list):
                value = tuple(value)
            read_csv_kw.append((key, value))
        return (cls, tuple(infiles), use_indexes, tuple(read_csv_kw))

    def __init__(self, calc_pdist=False, outfile=None,
        interactive=False, **kwargs):
        """
        Arguments:
          infile{s} (list): Path(s) to input file(s); may contain
            environment variables and wildcards
          use_indexes (list): Residue indexes to select from DataFrame,
            once DataFrame has already been loaded
          calc_pdist (bool): Calculate probability distribution
            using :meth:`calc_pdist`
          dataset_cache (dict): Cache of previously-loaded Datasets
          interactive (bool): Provide iPython prompt and reading and
            processing data
          verbose (int): Level of verbose output
          kwargs (dict): Additional keyword arguments
        """

        # Process arguments
        verbose = kwargs.get("verbose", 1)
        self.dataset_cache = kwargs.get("dataset_cache", None)

        # Read data
        self.sequence_df = self.read(**kwargs)

        # Cut data
        if "use_indexes" in kwargs:
            use_indexes = np.array(kwargs.pop("use_indexes"), np.int)
            res_index = np.array([int(i.split(":")[1])
              for i in self.sequence_df.index.values])
            self.sequence_df = self.sequence_df[np.in1d(res_index,use_indexes)]

        # Calculate probability distribution
        if calc_pdist:
            self.pdist_df = self.calc_pdist(df=self.sequence_df, **kwargs)

        # Output data
        if verbose >= 2:
            wiprint("Processed sequence DataFrame:")
            print(self.sequence_df)
            if calc_pdist:
                print("Processed pdist DataFrame:")
                print(self.pdist_df)

        # Write data
        if outfile is not None:
            self.write(df=self.sequence_df, outfile=outfile, **kwargs)

        # Interactive prompt
        if interactive:
            embed()

    def _read_index(self, df, indexfile=None, **kwargs):
        """
        Reads index for sequence DataFrame.

        Arguments:
          df (DataFrame): Nascent sequence DataFrame
          indexfile (str): Path to index file; may contain environment
            variables
          verbose (int): Level of verbose output
          kwargs (dict): Additional keyword arguments

        Returns:
          DataFrame: sequence DataFrame with updated index
        """
        from os.path import expandvars
        import re

        # Process arguments
        verbose = kwargs.get("verbose", 1)
        re_res = re.compile("[a-zA-Z]+:?[0-9]+")

        if indexfile is not None:
            indexfile = expandvars(indexfile)
            if verbose >= 1:
                wiprint("""Loading residue indexes from '{0}'
                        """.format(indexfile))
            res_index = np.loadtxt(indexfile, dtype=np.str).flatten()
            df.set_index(res_index, inplace=True)
            df.index.name = "residue"
        elif (sum([1 for a in [re_res.match(str(b)) for b in df.index.values]
              if a is not None]) == len(df.index.values)):
            df.index.name = "residue"
        else:
            df.index.name = "index"

        return df

    def read(self, **kwargs):
        """
        Reads sequence from one or more *infiles* into a DataFrame.

        Extends :class:`Dataset<myplotspec.Dataset.Dataset>` with
        option to read in residue indexes.
        """
        # Load Data
        df = super(SequenceDataset, self).read(**kwargs)

        # Load index, if applicable
        df = self._read_index(df, **kwargs)

        # Sort
        if df.index.name == "residue":
            df = df.loc[sorted(df.index.values,
                   key=lambda x: int(x.split(":")[1]))]
        else:
            df = df.loc[sorted(df.index.values)]

        return df

    def calc_pdist(self, **kwargs):
        """
        Calculates probability distribution across sequence.

        Arguments:
          df (DataFrame): DataFrame; probability distribution will be
            calculated for each column using rows as data points
          pdist_kw (dict): Keyword arguments used to configure
            probability distribution calculation
          pdist_kw[columns] (list): Columns for which to calculate
            probability distribution
          pdist_kw[mode] (ndarray, dict): Method of calculating
            probability distribution; eventually will support 'hist' for
            histogram and 'kde' for kernel density estimate, though
            presently only kde is implremented
          pdist_kw[grid] (ndarray, dict, optional): Grid on which to
            calculate kernel density estimate; may be a single ndarray
            that will be applied to all columns or a dictionary whose
            keys are column names and values are ndarrays corresponding
            to the grid for each column; for any column for which *grid*
            is not specified, a grid of 1000 points between the minimum
            value minus three times the standard deviation and the
            maximum value plots three times the standard deviation will
            be used
          pdist_kw[bandwidth] (float, dict, str, optional): Bandwidth to
            use for kernel density estimates; may be a single float that
            will be applied to all columns or a dictionary whose keys
            are column names and values are floats corresponding to the
            bandwidth for each column; for any column for which
            *bandwidth* is not specified, the standard deviation will be
            used; alternatively may be 'se', in which case the standard
            error of each value will be used
          verbose (int): Level of verbose output
          kwargs (dict): Additional keyword arguments

        Returns:
          dict: Dictionary whose keys are columns in *df* and values are
          DataFrames whose indexes are the *grid* for that column and
          contain a single column 'probability' containing the
          normalized probability at each grid point

        .. todo:
            - Implement flag to return single dataframe with single grid
        """
        from collections import OrderedDict
        from scipy.stats import norm
        import six

        # Process arguments
        verbose = kwargs.get("verbose", 1)
        df      = kwargs.get("df")
        if df is None:
            if hasattr(self, "sequence_df"):
                df = self.sequence_df
            else:
                raise()
        pdist_kw = kwargs.get("pdist_kw", {})
        columns = pdist_kw.get("columns",
          [a for a in df.columns.values
           if not a.endswith(" se")
           and str(df[a].dtype).startswith("float")])
        if isinstance(columns, six.string_types):
            columns = [columns]
        mode = pdist_kw.get("mode", "kde")

        if mode == "kde":

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
            for column, series in df[columns].iteritems():
                if column in grid:
                    grid[column] = np.array(grid[column])
                elif all_grid is not None:
                    grid[column] = all_grid
                else:
                    grid[column] = np.linspace(
                      series.min() - 3 * series.std(),
                      series.max() + 3 * series.std(), 1000)

            # Prepare bandwidths:
            bandwidth = pdist_kw.pop("bandwidth", None)
            if bandwidth is None:
                all_bandwidth = None
                bandwidth = {}
            elif (isinstance(bandwidth, six.string_types)
            and   bandwidth.lower() == "se"):
                all_bandwidth = None
                bandwidth = "se"
            elif isinstance(bandwidth, float):
                all_bandwidth = float(bandwidth)
                bandwidth = {}
            elif isinstance(bandwidth, dict):
                all_bandwidth = None
                pass
            for column, series in df[columns].iteritems():
                if column in bandwidth:
                    bandwidth[column] = float(bandwidth[column])
                elif all_bandwidth is not None:
                    bandwidth[column] = all_bandwidth
                else:
                    bandwidth[column] = series.std()

            # Calculate probability distributions
            pdist = OrderedDict()
            for column in columns:
                if verbose >= 1:
                    print("calculating probability distribution of "
                    "{0} using a kernel density estimate".format(column))
                g = grid[column]
                b = bandwidth[column]
                pdf = np.zeros_like(g)
                for residue, row in df.iterrows():
                    if np.isnan(row[column]):
                        continue
#                    c = norm(loc=b[column], scale=b[column+" se"])
                    pdf += norm(loc=row[column], scale=b).pdf(g)
                pdf /= pdf.sum()
                series_pdist = pd.DataFrame(pdf, index=g,
                  columns=["probability"])
                series_pdist.index.name = column
                pdist[column] = series_pdist
        else:
            raise Exception("only kde is currently supported")


        return pdist

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

class SAXSDataset(Dataset):
    """
    Represents Small-Angle X-ray Scattering Data.
    """

    def scale(self, scale, **kwargs):
        """
        Scales SAXS intensity, either by a constant or to match the
        intensity of a target dataset.

        Arguments:
          scale (float, str): If float, proportion by which to scale
            intensity; if str, path to input file to which intensity
            will be scaled, may contain environment variables
          curve_fit_kw (dict): Keyword arguments passed to
            scipy.optimize.curve_fit (scale to match target dataset
            only)
          verbose (int): Level of verbose output
          kwargs (dict): Additional keyword arguments
        """
        from os.path import expandvars, isfile
        from scipy.interpolate import interp1d
        from scipy.optimize import curve_fit
        import six

        verbose = kwargs.get("verbose", 1)

        # Scale by constant
        if (isinstance(scale, float)
        or (isinstance(scale, int) and not isinstance(scale, bool))):
            scale = float(scale)
        # Scale to match target
        elif isinstance(scale, six.string_types):
            if not isfile(expandvars(scale)):
                if verbose >= 1:
                   wiprint("scale target '{0}' ".format(scale) +
                          "not found, not scaling.")
                return

            # Prepare target
            target = self.load_dataset(infile=expandvars(scale),
                       loose=True).dataframe
            if "intensity_se" in target.columns.values:
                scale_se = True
            else:
                scale_se = False
            target_q = np.array(target.index.values, np.float64)
            target_I = np.array(target["intensity"], np.float64)
            if scale_se:
                target_Ise = np.array(target["intensity_se"], np.float64)

            # Prepare own values over x range of target
            template = self.dataframe
            template_q  = np.array(template.index.values, np.float64)
            template_I  = np.array(template["intensity"].values, np.float64)
            indexes = np.logical_and(template_q > target_q.min(),
                                     template_q < target_q.max())
            template_q = template_q[indexes]
            template_I = template_I[indexes]

            # Update target
            target_I = interp1d(target_q, target_I, kind="cubic")(template_q)
            if scale_se:
                target_Ise = interp1d(target_q, target_Ise,
                               kind="cubic")(template_q)

            def scale_I(_, a):
                return a * template_I
#            curve_fit_kw = dict(p0=(1), bounds=(0.0,0.35))
            curve_fit_kw = dict(p0=(1))     # Not clear why bounds broke
            curve_fit_kw.update(kwargs.get("curve_fit_kw", {}))
            if scale_se:
                curve_fit_kw["sigma"] = target_Ise
            scale = curve_fit(scale_I, template_q, target_I,
                      **curve_fit_kw)[0][0]
        # 'scale' argument not understood
        else:
            if verbose >= 1:
                wiprint("scale target '{0}' ".format(scale) +
                      "not understood, not scaling.")
            return

        if verbose >= 1:
            wiprint("scaling by factor of {0}".format(scale))
        self.dataframe["intensity"] *= scale
        if "intensity_se" in self.dataframe.columns.values:
            self.dataframe["intensity_se"] *= scale

        return scale

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
        from .myplotspec import multi_get_copy

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
        from .myplotspec import multi_get_copy

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

class HSQCDataset(Dataset):
    """
    Represents NMR HSQC data.

    Attributes:
      hsqc_df (DataFrame): DataFrame whose two-dimensional index
        corresponds to hydrogen and nitrogen chemical shift in ppm and
        whose columns correspond to intensity
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
        help_message = """Process HSQC data"""
        if isinstance(parser_or_subparsers, argparse.ArgumentParser):
            parser = parser_or_subparsers
        elif isinstance(parser_or_subparsers, argparse._SubParsersAction):
            parser = parser_or_subparsers.add_parser(
              name        = "hsqc",
              description = help_message,
              help        = help_message)
        elif parser is None:
            parser = argparse.ArgumentParser(
              description = help_message)

        # Defaults
        if parser.get_default("cls") is None:
            parser.set_defaults(cls=HSQCDataset)

        # Arguments unique to this class
        arg_groups = {ag.title: ag for ag in parser._action_groups}

        # Input arguments
        action_group = arg_groups.get("action",
          parser.add_argument_group("action"))
        try:
            action_group.add_argument(
              "-hoffset",
              required = False,
              type     = float,
              default  = 0,
              help     = """Offset added to 1H dimension """)
        except argparse.ArgumentError:
            pass
        try:
            action_group.add_argument(
              "-noffset",
              required = False,
              type     = float,
              default  = 0,
              help     = """Offset added to 15N dimension """)
        except argparse.ArgumentError:
            pass

        # Arguments inherited from superclass
        Dataset.construct_argparser(parser)

        return parser

    def __init__(self, hoffset=0, noffset=0, outfile=None,
        interactive=False, **kwargs):
        """
        Arguments:
          infile{s} (list): Path(s) to input file(s); may contain
            environment variables and wildcards
          hoffset (float, optional): Offset added to 1H dimension
          noffset (float, optional): Offset added to 15N dimension
          outfile (str, optional): Path to output file; may contain
            environment variables
          interactive (bool): Provide iPython prompt and reading and
            processing data
          verbose (int): Level of verbose output
          kwargs (dict): Additional keyword arguments
        """

        # Process arguments
        verbose = kwargs.get("verbose", 1)

        # Load
        self.hsqc_df = self.read(**kwargs)

        # Offset 1H and 15N
        hydrogen = np.array(self.hsqc_df.index.levels[0].values)
        nitrogen = np.array(self.hsqc_df.index.levels[1].values)
        hydrogen += float(hoffset)
        nitrogen += float(noffset)
        self.hsqc_df.index = pd.MultiIndex.from_product(
          [hydrogen, nitrogen], names=["1H", "15N"])

        # Output to screen
        if verbose >= 2:
            print("Processed HSQC DataFrame:")
            print(self.hsqc_df)

        # Write data
        if outfile is not None:
            self.write(df=self.hsqc_df, outfile=outfile, **kwargs)

        # Interactive prompt
        if interactive:
            embed()

    def _read_nmr(self, infile, **kwargs):
        """
        Reads DataFrame from NMR formats supported by nmrglue.

        Arguments:
          infile (str): Path to input file; may contain environment
            variables
          read_csv_kw (dict): Keyword arguments passed to
            :func:`read_csv<pandas.read_csv>`
          verbose (int): Level of verbose output
          kwargs (dict): Additional keyword arguments

        Returns:
          DataFrame: DataFrame
        """
        import nmrglue
        from os.path import expandvars

        # Process arguments
        verbose = kwargs.get("verbose", 1)
        infile = expandvars(infile)

        # Read DataFrame
        if verbose >= 1:
            wiprint("""Reading DataFrame from '{0}' """.format(infile))
        parameters, intensity = nmrglue.pipe.read(infile)
        hydrogen  = np.array(nmrglue.pipe.make_uc(parameters, intensity,
                       dim=1).ppm_scale(), np.float32)
        nitrogen  = np.array(nmrglue.pipe.make_uc(parameters, intensity,
                        dim=0).ppm_scale(), np.float32)

        index = pd.MultiIndex.from_product([nitrogen, hydrogen],
          names=["15N", "1H"])
        df = pd.DataFrame(data=intensity.flatten(), index=index,
          columns=["intensity"])
        df = df.swaplevel(0, 1)
        df = df.sortlevel()

        return df

    def read(self, **kwargs):
        """
        Reads HSQC data from one or more *infiles* into a DataFrame.

        Arguments:
          infile{s} (str): Path(s) to input file(s); may contain
            environment variables and wildcards
          dataframe_kw (dict): Keyword arguments passed to
            :class:`DataFrame<pandas.DataFrame>` (hdf5 only)
          read_csv_kw (dict): Keyword arguments passed to
            :func:`read_csv<pandas.read_csv>` (text only)
          verbose (int): Level of verbose output
          kwargs (dict): Additional keyword arguments

        Returns:
          DataFrame: DataFrame
        """
        import re
        from .myplotspec import multi_pop_merged

        # Process arguments
        infile_args = multi_pop_merged(["infile", "infiles"], kwargs)
        infiles = self.infiles = self.process_infiles(infiles=infile_args)
        if len(infiles) == 0:
            raise Exception(sformat("""No infiles found matching
            '{0}'""".format(infile_args)))
        elif len(infiles) > 1:
            raise Exception("HSQCDataset only supports a single infile")
        re_h5 = re.compile(
          r"^(?P<path>(.+)\.(h5|hdf5))((:)?(/)?(?P<address>.+))?$",
          flags=re.UNICODE)

        # Load Data
        infile = infiles[0]
        if infile.endswith(".ft"):
            df = self._read_nmr(infile, **kwargs)
        elif re_h5.match(infile):
            df = self._read_hdf5(infile, **kwargs)
        else:
            df = self._read_text(infile, **kwargs)
        df = df.astype(np.float32)

        return df

class PeakListDataset(SequenceDataset):
    """
    Represents an NMR peak list data.
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
        help_message = """Process NMR peak list data"""
        if isinstance(parser_or_subparsers, argparse.ArgumentParser):
            parser = parser_or_subparsers
        elif isinstance(parser_or_subparsers, argparse._SubParsersAction):
            parser = parser_or_subparsers.add_parser(
              name        = "peaklist",
              description = help_message,
              help        = help_message)
        elif parser is None:
            parser = argparse.ArgumentParser(
              description = help_message)

        # Defaults
        if parser.get_default("cls") is None:
            parser.set_defaults(cls=PeakListDataset)

        # Arguments unique to this class
        arg_groups = {ag.title: ag for ag in parser._action_groups}

        # Input arguments
        input_group  = arg_groups.get("input",
          parser.add_argument_group("input"))
        try:
            input_group.add_argument(
              "-delays",
              dest     = "delays",
              metavar  = "DELAY",
              nargs    = "+",
              type     = float,
              help     = """delays for each infile, if infiles represent a
                         series; number of delays must match number of
                         infiles""")
        except argparse.ArgumentError:
            pass

        # Action arguments
        action_group = arg_groups.get("action",
          parser.add_argument_group("action"))
        try:
            action_group.add_argument(
              "-relax",
              dest     = "calc_relax",
              type     = str,
              nargs    = "?",
              default  = None,
              const    = "r1",
              help     = """Calculate relaxation rates and standard errors; may
                         additionally specify kind of relaxation being
                         measured (e.g. r1, r2)""")
        except argparse.ArgumentError:
            pass


        # Arguments inherited from superclass
        SequenceDataset.construct_argparser(parser)

        return parser

    def __init__(self, delays=None, calc_relax=False, calc_pdist=False,
        outfile=None, interactive=False, **kwargs):
        """
        Arguments:
          infile{s} (list): Path(s) to input file(s); may contain
            environment variables and wildcards
          delays (list): Delays corresponding to series of infiles; used to
            name columns of merged sequence DataFrame
          use_indexes (list): Residue indexes to select from DataFrame,
            once DataFrame has already been loaded
          calc_pdist (bool): Calculate probability distribution
          pdist_kw (dict): Keyword arguments used to configure
            probability distribution calculation
          dataset_cache (dict): Cache of previously-loaded Datasets
          verbose (int): Level of verbose output
          kwargs (dict): Additional keyword arguments
        """

        # Process arguments
        verbose = kwargs.get("verbose", 1)
        self.dataset_cache = kwargs.get("dataset_cache", None)

        # Read data
        if delays is not None:
            kwargs["infile_column_prefixes"] = ["{0:3.1f} ms".format(delay)
              for delay in delays]
        self.sequence_df = self.read(**kwargs)

        # Cut data
        if "use_indexes" in kwargs:
            use_indexes = np.array(kwargs.pop("use_indexes"))
            res_index = np.array([int(i.split(":")[1])
              for i in self.sequence_df.index.values])
            self.sequence_df = self.sequence_df[np.in1d(res_index,use_indexes)]

        # Calculate relaxation
        if calc_relax:
            relax_kw = kwargs.pop("relax_kw", {})
            relax_kw["kind"] = calc_relax
            self.sequence_df = self.calc_relax(df=self.sequence_df,
              relax_kw=relax_kw, **kwargs)

        # Calculate probability distribution
        if calc_pdist:
            self.pdist_df = self.calc_pdist(df=self.sequence_df, **kwargs)

        # Output data
        if verbose >= 2:
            print("Processed sequence DataFrame:")
            print(self.sequence_df)
            if calc_pdist:
                print("Processed pdist DataFrame:")
                print(self.pdist_df)

        # Write data
        if outfile is not None:
            self.write(df=self.sequence_df, outfile=outfile, **kwargs)

        # Interactive prompt
        if interactive:
            embed()

    def read(self, **kwargs):
        """
        Reads sequence from one or more *infiles* into a DataFrame.

        Extends :class:`Dataset<myplotspec.Dataset.Dataset>` with
        option to read in residue indexes.
        """
        from os import devnull
        import re
        from subprocess import Popen, PIPE
        from .myplotspec import multi_pop_merged

        # Functions
        def convert_name(name):
            return "{0}:{1}".format(name[-4:-1].upper(),name[2:-4])

        # Process arguments
        infile_args = multi_pop_merged(["infile", "infiles"], kwargs)
        infiles = self.infiles = self.process_infiles(infiles=infile_args)
        if len(infiles) == 0:
            raise Exception(sformat("""No infiles found matching
            '{0}'""".format(infile_args)))
        re_h5 = re.compile(
          r"^(?P<path>(.+)\.(h5|hdf5))((:)?(/)?(?P<address>.+))?$",
          flags=re.UNICODE)
        infile_column_prefixes = kwargs.get("infile_column_prefixes",
          range(len(infiles)))

        # Load Data
        dfs = []
        for infile in infiles:
            if re_h5.match(infile):
                df = self._read_hdf5(infile, **kwargs)
            else:
                with open(devnull, "w") as fnull:
                    header = " ".join(Popen("head -n 1 {0}".format(infile),
                      stdout=PIPE, stderr=fnull, shell=True
                      ).stdout.read().strip().split("\t"))
                ccpnmr_header = sformat("""Number # Position F1 Position F2
                    Assign F1 Assign F2 Height Volume Line Width F1 (Hz) Line
                    Width F2 (Hz) Merit Details Fit Method Vol. Method""")
                if (header == ccpnmr_header):
                    read_csv_kw = dict(
                      index_col = None,
                      delimiter = "\t",
                      dtype = {
                        "Position F1": np.float32,
                        "Position F2": np.float32,
                        "Assign F1":   np.str,
                        "Height":      np.float32,
                        "Volume":      np.float32},
                      usecols = [
                        "Position F1",
                        "Position F2",
                        "Assign F1",
                        "Height",
                        "Volume"],
                      converters = {
                        "Assign F1": convert_name})
                    read_csv_kw.update(kwargs.get("read_csv_kw", {}))
                    kwargs["read_csv_kw"] = read_csv_kw
                    df = self._read_text(infile, **kwargs)
                    df.columns = ["1H", "15N", "residue", "height", "volume"]
                    df.set_index("residue", inplace=True)
                else:
                    df = self._read_text(infile, **kwargs)
            dfs.append(df)
        if len(dfs) == 1:
            df = dfs[0]
        else:
            df = dfs[0][["1H", "15N"]]
            if len(dfs) != len(infile_column_prefixes):
                raise Exception(sformat("""Numb of infile column prefixes
                  must match number of provided infiles"""))
            for df_i, prefix_i in zip(dfs, infile_column_prefixes):
                df["{0} height".format(prefix_i)] = df_i["height"]
                df["{0} volume".format(prefix_i)] = df_i["volume"]
        self.dfs = dfs

        # Sort
        if df.index.name == "residue":
            df = df.loc[sorted(df.index.values,
                   key=lambda x: int(x.split(":")[1]))]
        else:
            df = df.loc[sorted(df.index.values)]

        return df

    def calc_relax(self, **kwargs):
        """
        Calculates relaxation rates.

        Arguments:
          df (DataFrame): DataFrame; probability distribution will be
            calculated for each column using rows as data points
          relax_kw (dict): Keyword arguments used to configure
            relaxation rate calculation
          relax_kw[kind] (str): Kind of relaxation rate being
            calculated; will be used to name column
          relax_kw[intensity_method] (str): Metric to use for peak
            instensity; may be 'height' (default) or 'volume'
          relax_kw[error_method] (str): Metric to use for error
            calculation; may be 'rmse' for root-mean-square error
            (default) or 'mae' for mean absolute error
          relax_kw[n_synth_datasets] (int): Number of synthetic datasets
            to use for error calculation

        Returns:
          DataFrame: Sequence DataFrame with additional columns for
          relaxation rate and standard error
        """
        import re
        from scipy.optimize import curve_fit

        # Process arguments
        verbose = kwargs.get("verbose", 1)
        df      = kwargs.get("df")
        if df is None:
            if hasattr(self, "sequence_df"):
                df = self.sequence_df
            else:
                raise()
        relax_kw = kwargs.get("relax_kw", {})
        kind             = relax_kw.get("kind",             "r1")
        intensity_method = relax_kw.get("intensity_method", "height")
        error_method     = relax_kw.get("error_method",     "mae")
        n_synth_datasets = relax_kw.get("n_synth_datasets", 1000)

        # Calculate relaxation rates
        re_column = re.compile("^(?P<delay>\d+\.?\d*?) ms {0}".format(
          intensity_method))
        columns = [c for c in df.columns.values if re_column.match(c)]
        delays = np.array([re.match(re_column, c).groupdict()["delay"]
                   for c in columns], np.float) / 1000

        def calc_relax_rate(residue, **kwargs):
            """
            """
            from . import multiprocess_map

            if verbose >= 1:
                wiprint("""Calculating {0} relaxation rate for {1}""".format(
                kind, residue.name))

            def model_function(delay, intensity, relaxation):
                return intensity * np.exp(-1 * delay * relaxation)

            I = np.array(residue.filter(columns, np.float64))
            I0, R = curve_fit(model_function, delays, I, p0=(I[0], 1.0))[0]

            # Calculate error
            if error_method == "rmse":
                error = np.sqrt(np.mean((I-model_function(delays, I0, R))**2))
            elif error_method == "mae":
                error = np.mean(np.sqrt((I-model_function(delays, I0, R))**2))

            # Construct synthetic relaxation profiles
            synth_datasets = np.zeros((n_synth_datasets, I.size))
            for i, I_mean in enumerate(model_function(delays, I0, R)):
                synth_datasets[:, i] = np.random.normal(I_mean, error,
                                         n_synth_datasets)

            def synth_fit_decay(synth_intensity):
                try:
                    synth_I0, synth_R = curve_fit(model_function, delays,
                      synth_intensity, p0=(I0, R))[0]
                    return synth_R
                except RuntimeError:
                    if verbose >= 1:
                        wiprint("""Unable to calculate standard error for {0}
                                """.format(residue.name))
                    return np.nan

            # Calculate standard error
            synth_Rs = multiprocess_map(synth_fit_decay, synth_datasets, 16)
            R_se = np.std(synth_Rs)

            return pd.Series([I0, R, R_se])

        # Calculate relaxation rates and standard errors
        fit = df.apply(calc_relax_rate, axis=1)

        # Format and return
        fit.columns = ["I0", kind, kind + " se"]
        df = df.join(fit)
        return df

class RelaxSequenceDataset(SequenceDataset):
    """
    Represents NMR relaxation data as a function of residue number.
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
        help_message = """Process relaxation data"""
        if isinstance(parser_or_subparsers, argparse.ArgumentParser):
            parser = parser_or_subparsers
        elif isinstance(parser_or_subparsers, argparse._SubParsersAction):
            parser = parser_or_subparsers.add_parser(
              name        = "relax",
              description = help_message,
              help        = help_message)
        elif parser is None:
            parser = argparse.ArgumentParser(
              description = help_message)

        # Defaults
        if parser.get_default("cls") is None:
            parser.set_defaults(cls=RelaxSequenceDataset)

        # Arguments unique to this class

        # Arguments inherited from superclass
        SequenceDataset.construct_argparser(parser)

        return parser

    def __init__(self, calc_pdist=False, outfile=None, interactive=False,
        **kwargs):
        """
        Arguments:
          infile{s} (list): Path(s) to input file(s); may contain
            environment variables and wildcards
          use_indexes (list): Residue indexes to select from DataFrame,
            once DataFrame has already been loaded
          calc_pdist (bool): Calculate probability distribution
          pdist_kw (dict): Keyword arguments used to configure
            probability distribution calculation
          dataset_cache (dict): Cache of previously-loaded Datasets
          verbose (int): Level of verbose output
          kwargs (dict): Additional keyword arguments
        """

        # Process arguments
        verbose = kwargs.get("verbose", 1)
        self.dataset_cache = kwargs.get("dataset_cache", None)

        # Read data
        self.sequence_df = self.read(**kwargs)

        # Cut data
        if "use_indexes" in kwargs:
            use_indexes = np.array(kwargs.pop("use_indexes"), np.int)
            res_index = np.array([int(i.split(":")[1])
              for i in self.sequence_df.index.values])
            self.sequence_df = self.sequence_df[np.in1d(res_index,use_indexes)]

        # Calculate r2/r1 ratio
        if "r1" in self.sequence_df and "r2" in self.sequence_df:
            self.sequence_df["r2/r1"] = self.sequence_df["r2"] / \
                                        self.sequence_df["r1"]
            self.sequence_df["r2/r1 se"] = np.sqrt((self.sequence_df["r2 se"] /
              self.sequence_df["r2"]) ** 2 + (self.sequence_df["r1 se"] /
              self.sequence_df["r1"]) ** 2) * self.sequence_df["r2/r1"]

        # Calculate probability distribution
        if calc_pdist:
            self.pdist_df = self.calc_pdist(df=self.sequence_df, **kwargs)

        # Output data
        if verbose >= 2:
            if verbose >= 1:
                print("Processed sequence DataFrame:")
                print(self.sequence_df)

        # Write data
        if outfile is not None:
            self.write(df=self.sequence_df, outfile=outfile, **kwargs)

        # Interactive prompt
        if interactive:
            embed()

    def write_for_relax(self, outfile, **kwargs):
        """
        Writes sequence DataFrame in format readable by relax.
        """
        from os.path import expandvars
        from . import three_one

        # Process arguments
        df = kwargs.get("df")
        if df is None:
            if hasattr(self, "sequence_df"):
                df = self.sequence_df
            else:
                raise()
        outfile = expandvars(outfile)
        res_index = np.array([int(i.split(":")[1])
          for i in df.index.values])
        res_code = np.array([three_one(i.split(":")[0])
          for i in df.index.values])

        df["index"] = res_index
        df["code"] = res_code

        with open(outfile + "r1", "w") as r1_file:
            with open(outfile + "r2", "w") as r2_file:
                with open(outfile + "noe", "w") as noe_file:
                    for _, row in df.iterrows():
                        r1_file.write("{0}\t{1}\t{2}\t{3}\n".format(
                          row["code"],row["index"],row["r1"],row["r1 se"]))
                        r2_file.write("{0}\t{1}\t{2}\t{3}\n".format(
                          row["code"],row["index"],row["r2"],row["r2 se"]))
                        noe_file.write("{0}\t{1}\t{2}\t{3}\n".format(
                          row["code"],row["index"],row["noe"],row["noe se"]))

class IREDRelaxDataset(RelaxSequenceDataset):
    """
    Represents iRED NMR relaxation data as a function of residue number.
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
        help_message = """Process relaxation data calculated from MD simulation
            using the iRED method as implemented in cpptraj; treat infiles as
            independent simulations; processed results are the average across
            the simulations, including standard errors calculated using
            standard deviation"""
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
            parser.set_defaults(cls=IREDRelaxDataset)

        # Arguments unique to this class
        arg_groups = {ag.title: ag for ag in parser._action_groups}

        # Input arguments
        input_group = arg_groups.get("input",
          parser.add_argument_group("input"))
        try:
            input_group.add_argument(
              "-infiles",
              required = True,
              dest     = "infiles",
              metavar  = "INFILE",
              nargs    = "+",
              type     = str,
              help     = """File(s) from which to load data; may be text or
                         hdf5; if text, may be pandas-formatted DataFrames, or
                         may be cpptraj-formatted iRED output; may contain
                         environment variables and wildcards""")
        except argparse.ArgumentError:
            pass

        # Arguments inherited from superclass
        RelaxSequenceDataset.construct_argparser(parser)

        return parser

    @staticmethod
    def average_independent(relax_dfs=None, order_dfs=None, **kwargs):
        """
        Calculates the average and standard error of a set of independent
        iRED datasets.

        Arguments:
          relax_dfs (list): DataFrames containing data from relax infiles
          order_dfs (list): DataFrames containing data from order infiles
          kwargs (dict): Additional keyword arguments

        Returns:
          df (DataFrame): Averaged DataFrame including relax and order
        """

        # Process arguments
        verbose = kwargs.get("verbose", 1)
        df = pd.DataFrame()

        # Process relaxation
        if len(relax_dfs) == 1:
            if verbose >= 1:
                wiprint("""Single relaxation infile provided; skipping error
                        calculation""")
            df["r1"]         = relax_dfs[0]["r1"]
            if "r1 se" in relax_dfs[0]:
                df["r1 se"]  = relax_dfs[0]["r1 se"]
            df["r2"]         = relax_dfs[0]["r2"]
            if "r2 se" in relax_dfs[0]:
                df["r2 se"]  = relax_dfs[0]["r2 se"]
            df["noe"]        = relax_dfs[0]["noe"]
            if "noe se" in relax_dfs[0]:
                df["noe se"] = relax_dfs[0]["noe se"]
        elif len(relax_dfs) >= 2:
            if verbose >= 1:
                wiprint("""Calculating mean and standard error of {0}
                        relaxation infiles""".format(len(relax_dfs)))
            n_relax_dfs  = len(relax_dfs)
            relax_dfs    = pd.concat(relax_dfs)
            relax_mean   = relax_dfs.groupby(level=0).mean()
            relax_se     = relax_dfs.groupby(level=0).std() / \
                             np.sqrt(n_relax_dfs)
            df["r1"]     = relax_mean["r1"]
            df["r1 se"]  = relax_se["r1"]
            df["r2"]     = relax_mean["r2"]
            df["r2 se"]  = relax_se["r2"]
            df["noe"]    = relax_mean["noe"]
            df["noe se"] = relax_se["noe"]

        # Process order parameters
        if len(order_dfs) == 1:
            if verbose >= 1:
                wiprint("""Single order parameter infile provided; skipping
                        error calculation""")
            df["s2"] = order_dfs[0]["s2"]
            if "s2 se" in order_dfs[0]:
                df["s2 se"] = order_dfs[0]["s2 se"]
        elif len(order_dfs) >= 2:
            if verbose >= 1:
                wiprint("""Calculating mean and standard error of {0} order
                        parameter infiles""".format(len(order_dfs)))
            n_order_dfs = len(order_dfs)
            order_dfs   = pd.concat(order_dfs)
            order_mean  = order_dfs.groupby(level=0).mean()
            order_se    = order_dfs.groupby(level=0).std() / \
                          np.sqrt(n_order_dfs)
            df["s2"]    = order_mean["s2"]
            df["s2 se"] = order_se["s2"]

        # Sort by index
        if df.index.name == "residue":
            df = df.loc[sorted(df.index.values,
                   key=lambda x: int(x.split(":")[1]))]
        else:
            df = df.loc[sorted(df.index.values)]

        return df

    @staticmethod
    def _identify_infile(infile, **kwargs):
        """
        Determines if an infile contains iRED relaxation data, iRED
        order parameters, or neither

        Arguments:
          infile (str): Path to input file

        Returns:
          str: Kind of data in *infile*; may be 'ired_relax',
          'ired_order', or 'other'

        .. todo:
          - Identify pandas files and hdf5 files
        """
        import six
        from os import devnull
        import re
        from subprocess import Popen, PIPE

        re_t1t2noe = re.compile(
          r"^#Vec\s+[\w_]+\[T1\]\s+[\w_]+\[\T2\]\s+[\w_]+\[NOE\]$",
          flags=re.UNICODE)
        re_order = re.compile(
          r"^#Vec\s+[\w_]+\[S2\]$", flags=re.UNICODE)
        re_h5 = re.compile(
          r"^(?P<path>(.+)\.(h5|hdf5))((:)?(/)?(?P<address>.+))?$",
          flags=re.UNICODE)

        with open(devnull, "w") as fnull:
            header = Popen("head -n 1 {0}".format(infile),
              stdout=PIPE, stderr=fnull, shell=True).stdout.read().strip()
            if six.PY3:
                header = str(header, "utf-8")

        if re.match(re_t1t2noe, header):
            return "ired_relax"
        elif re.match(re_order, header):
            return "ired_order"
        elif re.match(re_h5, infile):
            return "hdf5"
        else:
            return "other"

    def _read_text(self, infile, **kwargs):
        """
        Reads iRED sequence DataFrame from text.

        Arguments:
          infile (str): Path to input file; may contain environment
            variables
          read_csv_kw (dict): Keyword arguments passed to
            :func:`read_csv<pandas.read_csv>`
          verbose (int): Level of verbose output
          kwargs (dict): Additional keyword arguments

        Returns:
          DataFrame: Sequence DataFrame
        """
        from os.path import expandvars

        # Process arguments
        verbose = kwargs.get("verbose", 1)
        infile = expandvars(infile)
        kind = IREDRelaxDataset._identify_infile(infile)

        if kind == "ired_relax":                    # Parse relaxation
            if verbose >= 1:
                wiprint("""Loading iRED relaxation data from '{0}'
                        """.format(infile))
            df = pd.read_csv(infile, delim_whitespace=True, header=0,
              index_col=0, names=["r1","r2","noe"])
            df["r1"] = 1 / df["r1"]
            df["r2"] = 1 / df["r2"]
        elif kind == "ired_order":                  # Parse order parameters
            if verbose >= 1:
                wiprint("""Loading iRED order parameter data from '{0}'
                        """.format(infile))
            df = pd.read_csv(infile, delim_whitespace=True, header=0,
              index_col=0, names=["s2"])
        else:                                       # Parse other
            df = super(IREDRelaxDataset, self)._read_text(infile, **kwargs)
        df = self._read_index(df, **kwargs)

        return df

    def read(self, **kwargs):
        """
        Reads iRED sequence data from one or more *infiles* into a
        DataFrame.

        *infiles* may contain relaxation data, order parameters, or
        both. If more than one *infile* is provided, the resulting
        DataFrame will contain their average, and the standard error
        will be calculated assuming the *infiles* represent independent
        samples.

        After generating the DataFrame from *infiles*, the index may be
        set by loading a list of residue names and numbers in the form
        ``XAA:#`` from *indexfile*. This is useful when loading data
        from files that do not specify residue names.

        Arguments:
          infile{s} (list): Path(s) to input file(s); may contain
            environment variables and wildcards
          dataframe_kw (dict): Keyword arguments passed to
            :class:`DataFrame<pandas.DataFrame>` (hdf5 only)
          read_csv_kw (dict): Keyword arguments passed to
            :func:`read_csv<pandas.read_csv>` (text only)
          indexfile (str): Path to index file; may contain environment
            variables
          verbose (int): Level of verbose output
          kwargs (dict): Additional keyword arguments

        Returns:
          df (DataFrame): iRED sequence DataFrame
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
        relax_dfs = []
        order_dfs = []
        for infile in infiles:
            if re_h5.match(infile):
                df = self._read_hdf5(infile, **kwargs)
            else:
                df = self._read_text(infile, **kwargs)
            columns = df.columns.values
            if "r1" in columns and "r2" in columns and "noe" in columns:
                relax_dfs.append(df)
            if "s2" in columns:
                order_dfs.append(df)
            if not (("r1" in columns and "r2" in columns and "noe" in columns)
            or      ("s2" in columns)):
                raise Exception(sformat("""DataFrame loaded from '{0}' does not
                  appear to contain either relaxation ('r1', 'r2', 'noe') or
                  order parameter ('s2') columns""".format(infile)))

        # Average, if applicable
        df = self.average_independent(relax_dfs, order_dfs)

        return df

class IREDTimeSeriesDataset(TimeSeriesDataset, IREDRelaxDataset):
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
        IREDRelaxDataset.construct_argparser(parser)
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

class SAXSExperimentDataset(SAXSDataset):
    """
    Manages Small Angle X-ray Scattering experimental datasets.
    """

    def __init__(self, scale=False, **kwargs):
        """
        Arguments:
          infile (str): Path to input file, may contain environment
            variables
          verbose (int): Level of verbose output
          kwargs (dict): Additional keyword arguments
        """
        from os.path import expandvars

        # Load
        super(SAXSExperimentDataset, self).__init__(**kwargs)
        dataframe = self.dataframe

        # Scale
        if scale:
            self.scale(scale, **kwargs)

class SAXSDiffDataset(SAXSDataset):
    """
    Manages Small Angle X-ray Scattering difference datasets.
    """

#    @classmethod
#    def get_cache_key(cls, dataset_classes=None, *args, **kwargs):
#        """
#        Generates tuple of arguments to be used as key for dataset
#        cache.
#
#        Arguments documented under :func:`__init__`.
#        """
#        from .myplotspec import multi_get_copy
#
#        minuend_kw = multi_get_copy(["minuend", "minuend_kw"], kwargs, {})
#        minuend_class = dataset_classes[minuend_kw["kind"].lower()]
#        key = [cls, mask_cutoff, minuend_class.get_cache_key(**minuend_kw)]
#
#        subtrahend_kw = multi_get_copy(["subtrahend", "subtrahend_kw"],
#          kwargs, {})
#        if isinstance(subtrahend_kw, dict):
#            subtrahend_kw = [subtrahend_kw]
#        for sh_kw in subtrahend_kw:
#            sh_class = dataset_classes[sh_kw.pop("kind").lower()]
#            key.append(sh_class.get_cache_key(**sh_kw))
#
#        return tuple(key)

    def __init__(self, dataset_cache=None, **kwargs):
        """
        """
        from sys import exit
        from .myplotspec import multi_get_copy

        self.dataset_cache = dataset_cache

        minuend_kw    = multi_get_copy(["minuend", "minuend_kw"],
                          kwargs, {})
        subtrahend_kw = multi_get_copy(["subtrahend", "subtrahend_kw"],
                          kwargs, {})
        minuend    = self.load_dataset(loose=True, **minuend_kw)
        subtrahend = self.load_dataset(loose=True, **subtrahend_kw)
        m_I    = minuend.dataframe["intensity"]
        m_I_se = minuend.dataframe["intensity_se"]
        s_I    = subtrahend.dataframe["intensity"]
        s_I_se = subtrahend.dataframe["intensity_se"]
        diff_I    = (m_I - s_I)
        diff_I_se = np.sqrt(m_I_se**2 +s_I_se**2)
        diff_I.name = "intensity"
        diff_I_se.name = "intensity_se"
        self.dataframe = pd.concat([diff_I, diff_I_se], axis=1)

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

class H5Dataset(object):
    """
    Class for managing hdf5 datasets

    .. warning::
      Will be reimplemented or removed eventually
    """

    def __init__(self, **kwargs):
        """
        Arguments:
          infiles (list): List of infiles
          infile (str): Alternatively, single infile
        """
        self.default_address = kwargs.get("default_address", "")
        self.default_key     = kwargs.get("default_key",     "key")
        self.datasets        = {}
        self.attrs           = {}

        if   "infiles" in kwargs:
            self.load(infiles = kwargs.pop("infiles"))
        elif "infile"  in kwargs:
            self.load(infiles = [kwargs.pop("infile")])

    def load(self, infiles, **kwargs):
        """
        Loads data from h5 files.

        Arguments:
          infiles (list): infiles
        """
        from os.path import expandvars, isfile
        from h5py import File as h5
        import numpy as np
        import six

        for infile in infiles:
            if isinstance(infile, six.string_types):
                path    = expandvars(infile)
                address = self.default_address
                key     = self.default_key
            elif isinstance(infile, dict):
                path    = expandvars(infile.pop("path"))
                address = infile.pop("address", self.default_address)
                key     = infile.pop("key",     self.default_key)
            elif isinstance(infile, list):
                if len(infile) >= 1:
                    path = expandvars(infile[0])
                else:
                    raise OSError("Path to infile not provided")
                if len(infile) >= 2:
                    address = infile[1]
                else:
                    address = self.default_address
                if len(infile) >= 3:
                    key     = infile[2]
                else:
                    key     = self.default_key

            if not isfile(path):
                raise OSError("h5 file '{0}' does not exist".format(path))

            with h5(path) as in_h5:
                if address not in in_h5:
                    raise KeyError("Dataset {0}[{1}] not found".format(path,
                      address))
                dataset            = in_h5[address]
                self.datasets[key] = np.array(dataset)
                self.attrs[key]    = dict(dataset.attrs)
            print("Loaded Dataset {0}[{1}]; Stored at {2}".format(
              path, address, key))

#################################### MAIN #####################################
if __name__ == "__main__":
    import argparse

    # Prepare argument parser
    parser = argparse.ArgumentParser(
      description = """Processes datasets""")
    subparsers = parser.add_subparsers(
      dest        = "mode",
      description = "")

    Dataset.construct_argparser(subparsers)
    SequenceDataset.construct_argparser(subparsers)
    TimeSeriesDataset.construct_argparser(subparsers)
    RelaxSequenceDataset.construct_argparser(subparsers)
    IREDRelaxDataset.construct_argparser(subparsers)
    IREDTimeSeriesDataset.construct_argparser(subparsers)
    HSQCDataset.construct_argparser(subparsers)
    PeakListDataset.construct_argparser(subparsers)

    kwargs  = vars(parser.parse_args())
    kwargs.pop("cls")(**kwargs)
