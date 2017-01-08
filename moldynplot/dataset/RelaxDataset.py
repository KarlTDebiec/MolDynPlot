#!/usr/bin/python
# -*- coding: utf-8 -*-
#   moldynplot.dataset.RelaxDataset.py
#
#   Copyright (C) 2015-2017 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Represents NMR relaxation data as a function of residue number

.. todo:
  - Move relaxation error here
  - Move relaxation heteronuclear noe here
  - Move relaxation pre ratio here
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
from .SequenceDataset import SequenceDataset
from ..myplotspec import sformat, wiprint


################################### CLASSES ###################################
class RelaxDataset(SequenceDataset):
    """
    Represents NMR relaxation data as a function of residue number
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
        help_message = """Process NMR relaxation data"""
        if isinstance(parser_or_subparsers, argparse.ArgumentParser):
            parser = parser_or_subparsers
        elif isinstance(parser_or_subparsers, argparse._SubParsersAction):
            parser = parser_or_subparsers.add_parser(name="relax",
              description=help_message, help=help_message)
        elif parser_or_subparsers is None:
            parser = argparse.ArgumentParser(description=help_message)

        # Defaults
        if parser.get_default("cls") is None:
            parser.set_defaults(cls=RelaxDataset)

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
            res_index = np.array(
              [int(i.split(":")[1]) for i in self.sequence_df.index.values])
            self.sequence_df = self.sequence_df[
                np.in1d(res_index, use_indexes)]

        # Calculate r2/r1 ratio
        if "r1" in self.sequence_df and "r2" in self.sequence_df:
            self.sequence_df["r2/r1"] = self.sequence_df["r2"] / \
                                        self.sequence_df["r1"]
            self.sequence_df["r2/r1 se"] = np.sqrt(
              (self.sequence_df["r2 se"] / self.sequence_df["r2"]) ** 2 + (
                                                                              self.sequence_df[
                                                                                  "r1 se"] /
                                                                              self.sequence_df[
                                                                                  "r1"]) ** 2) * \
                                           self.sequence_df["r2/r1"]

        # Calculate probability distribution
        if calc_pdist:
            pdist_kw = kwargs.get("pdist_kw", {})
            self.pdist_df = self.calc_pdist(df=self.sequence_df,
              verbose=verbose, **pdist_kw)
            if verbose >= 2:
                print("Processed pdist DataFrame:")
                print(self.pdist_df)

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
                raise ()
        outfile = expandvars(outfile)
        res_index = np.array([int(i.split(":")[1]) for i in df.index.values])
        res_code = np.array(
          [three_one(i.split(":")[0]) for i in df.index.values])

        df["index"] = res_index
        df["code"] = res_code

        with open(outfile + "r1", "w") as r1_file:
            with open(outfile + "r2", "w") as r2_file:
                with open(outfile + "noe", "w") as noe_file:
                    for _, row in df.iterrows():
                        r1_file.write(
                          "{0}\t{1}\t{2}\t{3}\n".format(row["code"],
                            row["index"], row["r1"], row["r1 se"]))
                        r2_file.write(
                          "{0}\t{1}\t{2}\t{3}\n".format(row["code"],
                            row["index"], row["r2"], row["r2 se"]))
                        noe_file.write(
                          "{0}\t{1}\t{2}\t{3}\n".format(row["code"],
                            row["index"], row["noe"], row["noe se"]))


#################################### MAIN #####################################
if __name__ == "__main__":
    RelaxDataset.main()
