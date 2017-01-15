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
Represents paramagnetic relaxation enhancement (PRE) timeseries data
"""
################################### MODULES ###################################
from __future__ import (absolute_import, division, print_function,
    unicode_literals)

if __name__ == "__main__":
    __package__ = str("moldynplot.dataset")
    import moldynplot.dataset
from .SequenceDataset import SequenceDataset
from .TimeSeriesDataset import TimeSeriesDataset


################################### CLASSES ###################################
class PRETimeSeriesDataset(TimeSeriesDataset, SequenceDataset):
    """
    Represents paramagnetic relaxation enhancement (PRE) timeseries data
    """

    @staticmethod
    def construct_argparser(parser_or_subparsers=None, **kwargs):
        """
        Adds arguments to an existing argument parser,
        constructs a subparser, or constructs a new parser

        Arguments:
          parser_or_subparsers (ArgumentParser, _SubParsersAction,
            optional): If ArgumentParser, existing parser to which
            arguments will be added; if _SubParsersAction, collection
            of subparsers to which a new argument parser will be
            added; if None, a new argument parser will be generated
          kwargs (dict): Additional keyword arguments

        Returns:
          ArgumentParser: Argument parser or subparser
        """
        import argparse

        # Process arguments
        help_message = """Process NMR paramagnetic relaxation enhancement (PRE)
          timeseries data"""
        if isinstance(parser_or_subparsers, argparse.ArgumentParser):
            parser = parser_or_subparsers
        elif isinstance(parser_or_subparsers, argparse._SubParsersAction):
            parser = parser_or_subparsers.add_parser(name="pre",
              description=help_message, help=help_message)
        elif parser_or_subparsers is None:
            parser = argparse.ArgumentParser(description=help_message)

        # Defaults
        if parser.get_default("cls") is None:
            parser.set_defaults(cls=PRETimeSeriesDataset)

        # Arguments inherited from superclass
        TimeSeriesDataset.construct_argparser(parser)
        SequenceDataset.construct_argparser(parser)

        return parser

    def __init__(self, dt=None, downsample=None, outfile=None, **kwargs):
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
        if downsample:
            self.timeseries_df = self.downsample(df=self.timeseries_df,
              downsample=downsample, **kwargs)

        import pandas as pd
        pd.set_option('display.width', 200)
        import numpy as np

        distance_df = pd.DataFrame(data=self.timeseries_df.values,
          dtype=np.double, columns=self.timeseries_df.columns)

        k = 0.0123 # Å6 ns-2
        w = 800e6 # s-1
        tc = 4e-9 # s
        r2 = 13 # s
        t = 0.0005 # s

        rho2_df = k / (distance_df ** 6) * 1e9 * 1e9 # s-2
        rho2_df *= (4*tc + ((3*tc) / (1 + (w*(tc**2)))))
        I_I0_df = (r2 * np.exp(-1 * rho2_df * t)) / (r2 + rho2_df)

        block_kw = dict(min_n_blocks=2, max_cut=0.1, all_factors=False,
          fit_exp=True, fit_sig=False)
        block_kw.update(kwargs.get("block_kw", {}))

        distance_mean_df, block_averager = self.calc_mean(df=distance_df,
          verbose=verbose, **block_kw)
        distance_mean_df.columns = ["distance", "distance se"]
        rho2_mean_df, block_averager = self.calc_mean(df=rho2_df,
          verbose=verbose, **block_kw)
        rho2_mean_df.columns = ["rho2", "rho2 se"]
        I_I0_mean_df, block_averager = self.calc_mean(df=I_I0_df,
          verbose=verbose, **block_kw)
        I_I0_mean_df.columns = ["I/I0", "I/I0 se"]
        mean_df = distance_mean_df.join(I_I0_mean_df).join(rho2_mean_df)
        mean_df = self._read_index(df=mean_df,
          indexfile="/Volumes/KDebiecSSD/AMBER15IPQ/protein/MOCVNHLYSM"
                    "/ff15ipq10.3_spceb_T_tau_10.0ps/analysis"
                    "/MOCVNHLYSM_index.dat")
        print(mean_df)

        if outfile is not None:
            with open(outfile, "w") as out:
                out.write(
                  "#residue        distance distance se        I/I0     I/I0 "
                  "se        rho2     rho2 se\n")
                for residue in mean_df.index:
                    row = mean_df.loc[residue]
                    out.write(
                      "{0:12s} {1:11.3f} {2:11.3f} {3:11.3f} {4:11.3f} {"
                      "5:11.2f} {6:11.2f}\n".format(residue, row["distance"],
                        row["distance se"], row["I/I0"], row["I/I0 se"],
                        row["rho2"], row["rho2 se"]))


#################################### MAIN #####################################
if __name__ == "__main__":
    PRETimeSeriesDataset.main()
