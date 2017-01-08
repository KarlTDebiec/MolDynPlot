#!/usr/bin/python
# -*- coding: utf-8 -*-
#   moldynplot.dataset.IREDDataset.py
#
#   Copyright (C) 2015-2017 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Represents iRED NMR relaxation data as a function of residue number

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
from .RelaxDataset import RelaxDataset
from ..myplotspec import sformat, wiprint


################################### CLASSES ###################################
class IREDDataset(RelaxDataset):
    """
    Represents iRED NMR relaxation data as a function of residue number
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
                       simulation using the iRED method as implemented in
                       cpptraj"""
        if isinstance(parser_or_subparsers, argparse.ArgumentParser):
            parser = parser_or_subparsers
        elif isinstance(parser_or_subparsers, argparse._SubParsersAction):
            parser = parser_or_subparsers.add_parser(name="ired",
              description=help_message, help=help_message)
        elif parser_or_subparsers is None:
            parser = argparse.ArgumentParser(description=help_message)

        # Defaults
        if parser.get_default("class_") is None:
            parser.set_defaults(class_=IREDDataset)

        # Arguments unique to this class
        arg_groups = {ag.title: ag for ag in parser._action_groups}

        # Input arguments
        input_group = arg_groups.get("input",
          parser.add_argument_group("input"))
        try:
            input_group.add_argument("-infiles", required=True, dest="infiles",
              metavar="INFILE", nargs="+", type=str, help="""File(s) from
                which to load data; may be text or
                         hdf5; if text, may be pandas-formatted DataFrames, or
                         may be cpptraj-formatted iRED output; may contain
                         environment variables and wildcards""")
        except argparse.ArgumentError:
            pass

        # Arguments inherited from superclass
        RelaxDataset.construct_argparser(parser)

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
            df["r1"] = relax_dfs[0]["r1"]
            if "r1 se" in relax_dfs[0]:
                df["r1 se"] = relax_dfs[0]["r1 se"]
            df["r2"] = relax_dfs[0]["r2"]
            if "r2 se" in relax_dfs[0]:
                df["r2 se"] = relax_dfs[0]["r2 se"]
            df["noe"] = relax_dfs[0]["noe"]
            if "noe se" in relax_dfs[0]:
                df["noe se"] = relax_dfs[0]["noe se"]
        elif len(relax_dfs) >= 2:
            if verbose >= 1:
                wiprint("""Calculating mean and standard error of {0}
                        relaxation infiles""".format(len(relax_dfs)))
            n_relax_dfs = len(relax_dfs)
            relax_dfs = pd.concat(relax_dfs)
            relax_mean = relax_dfs.groupby(level=0).mean()
            relax_se = relax_dfs.groupby(level=0).std() / np.sqrt(n_relax_dfs)
            df["r1"] = relax_mean["r1"]
            df["r1 se"] = relax_se["r1"]
            df["r2"] = relax_mean["r2"]
            df["r2 se"] = relax_se["r2"]
            df["noe"] = relax_mean["noe"]
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
            order_dfs = pd.concat(order_dfs)
            order_mean = order_dfs.groupby(level=0).mean()
            order_se = order_dfs.groupby(level=0).std() / np.sqrt(n_order_dfs)
            df["s2"] = order_mean["s2"]
            df["s2 se"] = order_se["s2"]

        # Sort by index
        if df.index.name == "residue":
            df = df.loc[
                sorted(df.index.values, key=lambda x: int(x.split(":")[1]))]
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
        re_order = re.compile(r"^#Vec\s+[\w_]+\[S2\]$", flags=re.UNICODE)
        re_h5 = re.compile(
          r"^(?P<path>(.+)\.(h5|hdf5))((:)?(/)?(?P<address>.+))?$",
          flags=re.UNICODE)

        with open(devnull, "w") as fnull:
            header = Popen("head -n 1 {0}".format(infile), stdout=PIPE,
              stderr=fnull, shell=True).stdout.read().strip()
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
        kind = IREDDataset._identify_infile(infile)

        if kind == "ired_relax":  # Parse relaxation
            if verbose >= 1:
                wiprint("""Loading iRED relaxation data from '{0}'
                        """.format(infile))
            df = pd.read_csv(infile, delim_whitespace=True, header=0,
              index_col=0, names=["r1", "r2", "noe"])
            df["r1"] = 1 / df["r1"]
            df["r2"] = 1 / df["r2"]
        elif kind == "ired_order":  # Parse order parameters
            if verbose >= 1:
                wiprint("""Loading iRED order parameter data from '{0}'
                        """.format(infile))
            df = pd.read_csv(infile, delim_whitespace=True, header=0,
              index_col=0, names=["s2"])
        else:  # Parse other
            df = super(IREDDataset, self)._read_text(infile, **kwargs)
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
            if not (
                ("r1" in columns and "r2" in columns and "noe" in columns) or (
                  "s2" in columns)):
                raise Exception(sformat("""DataFrame loaded from '{0}' does not
                  appear to contain either relaxation ('r1', 'r2', 'noe') or
                  order parameter ('s2') columns""".format(infile)))

        # Average, if applicable
        df = self.average_independent(relax_dfs, order_dfs)

        return df


#################################### MAIN #####################################
if __name__ == "__main__":
    IREDDataset.main()