#!/usr/bin/python
# -*- coding: utf-8 -*-
#   moldynplot.dataset.TimeSeriesYSpecDataset.py
#
#   Copyright (C) 2015-2017 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Processes and represents data that is a function of time
"""
################################### MODULES ###################################
from __future__ import (absolute_import, division, print_function,
    unicode_literals)

if __name__ == "__main__":
    __package__ = str("moldynplot.dataset")
    import moldynplot.dataset
from IPython import embed
import numpy as np
import pandas as pd
from ..myplotspec.YSpecDataset import YSpecDataset


################################### CLASSES ###################################
class TimeSeriesYSpecDataset(YSpecDataset):
    """
    Processes and represents data that is a function of time

    Attributes:
      timeseries_df (DataFrame): DataFrame whose index corresponds to
        time as represented by frame number or chemical time and whose
        columns are a series of quantities as a function of time.
    """
    help_groups = ["spec"]

    @classmethod
    def construct_argparser(cls, **kwargs):
        """
        Adds arguments to a nascent argument parser

        Arguments:
          kwargs (dict): Additional keyword arguments

        Returns:
          ArgumentParser: Argument parser or subparser
        """

        # Process arguments
        parser = cls.get_argparser(grouped_help=True, **kwargs)
        if parser.get_default("cls") is None:
            parser.set_defaults(cls=cls)
        help_groups = kwargs.get("help_groups")
        arg_groups = {ag.title: ag for ag in parser._action_groups}
        input_group = arg_groups.get("input",
            parser.add_argument_group("input"))
        action_group = arg_groups.get("action",
            parser.add_argument_group("action"))
        output_group = arg_groups.get("output",
            parser.add_argument_group("output"))

        # Input arguments
        cls.add_argument(input_group, "-spec", dest="source_spec",
            metavar="SPEC", type=str, help="""file from which to load
            specification; see '--help
                     spec' for more information""")

        # Action arguments
        cls.add_argument(action_group, "-dt", type=float,
            help="time between frames")
        cls.add_argument(action_group, "-toffset", type=float,
            help="offset to add to index (time or frame number)")
        cls.add_argument(action_group, "-downsample", type=int,
            help="factor by which to downsample data")
        cls.add_argument(action_group, "--pdist", action="store_true",
            dest="calc_pdist",
            help="calculate probability distribution over timeseries")
        cls.add_argument(action_group, "--mean", action="store_true",
            dest="calc_mean",
            help="calculate mean and standard error over timeseries")

        # Arguments inherited from superclass
        super(TimeSeriesYSpecDataset, cls).construct_argparser(
            parser=parser, **kwargs)

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
        print(dt)

        # Read data
        self.timeseries_df = self.read(**kwargs)

        # Process data
        if dt:
            self.timeseries_df.set_index(
                self.timeseries_df.index.values * float(dt), inplace=True)
            self.timeseries_df.index.name = "time"
        if toffset:
            index_name = self.timeseries_df.index.name
            self.timeseries_df.set_index(
                self.timeseries_df.index.values + float(toffset), inplace=True)
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
        reduced = df.values[:df.shape[0] - (df.shape[0] % downsample), :]
        new_shape = (
            int(reduced.shape[0] / downsample), downsample, reduced.shape[1])
        index = np.reshape(
            df.index.values[:df.shape[0] - (df.shape[0] % downsample)],
            new_shape[:-1]).mean(axis=1)
        reduced = np.reshape(reduced, new_shape)

        # Downsample
        if downsample_mode == "mean":
            if verbose >= 1:
                wrapprint("downsampling by factor of {0} using mean".format(
                    downsample))
            reduced = np.squeeze(reduced.mean(axis=1))
        elif downsample_mode == "mode":
            if verbose >= 1:
                wrapprint("downsampling by factor of {0} using mode".format(
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
            wrapprint(
                """Calculating mean and standard error over timeseries""")

        # Single-level columns
        if df.columns.nlevels == 1:
            mean_df = pd.DataFrame(data=df.mean(axis=0))[0]
            block_averager = FPBlockAverager(df, **kwargs)
            if fit_exp and not fit_sig:
                errors = block_averager.parameters.loc[("exp", "a (se)")]
            elif fit_sig and not fit_exp:
                errors = block_averager.parameters.loc[("sig", "b (se)")]
            elif fit_exp and fit_sig:
                errors = block_averager.parameters.loc[("exp", "a (se)")]
            else:
                raise Exception()
            errors.index = errors.index.values + " se"
            mean_df = pd.concat([mean_df, errors])
            mean_df = mean_df[np.array(
                [[c, c + " se"] for c in df.columns.values]).flatten()]
            print("####################################")
            print(mean_df)
            print(df.columns.values)
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
            errors.index = pd.MultiIndex.from_tuples(
                map(eval, errors.index.values))
            errors.index = errors.index.set_levels(
                [c + " se" for c in errors.index.levels[1].values], level=1)
            mean_df = mean_df.squeeze().unstack().join(errors.unstack)
            print("####################################")
            a = mean_df.squeeze().unstack()
            print(a)
            b = errors.unstack()
            print(b)
            c = a.join(b)
            print(c)
        # Additional levels not tested
        else:
            raise Exception()

        # c = c[["r1", "r1 se", "r2", "r2 se", "noe", "noe se", "s2", "s2 se"]]
        # c = c.loc[sorted(c.index.values, key=lambda x: int(x.split(":")[1]))]
        # mean_df = c

        return mean_df, block_averager


#################################### MAIN #####################################
if __name__ == "__main__":
    TimeSeriesYSpecDataset.main()
