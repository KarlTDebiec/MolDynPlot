#!/usr/bin/python
# -*- coding: utf-8 -*-
#   moldynplot.dataset.HSQCDataset.py
#
#   Copyright (C) 2015-2017 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
.. todo:
  - Fix separation and ordering of argument groups: input, action, output
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
from ..myplotspec.Dataset import Dataset
from ..myplotspec import sformat, wiprint


################################### CLASSES ###################################
class HSQCDataset(Dataset):
    """
    Represents two-dimensional NMR data.

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
            parser = parser_or_subparsers.add_parser(name="hsqc",
                description=help_message, help=help_message)
        elif parser_or_subparsers is None:
            parser = argparse.ArgumentParser(description=help_message)

        # Defaults
        if parser.get_default("class_") is None:
            parser.set_defaults(class_=HSQCDataset)

        # Arguments unique to this class
        arg_groups = {ag.title: ag for ag in parser._action_groups}

        # Input arguments
        action_group = arg_groups.get("action",
            parser.add_argument_group("action"))
        try:
            action_group.add_argument("-hoffset", required=False, type=float,
                default=0, help="""Offset added to 1H dimension """)
        except argparse.ArgumentError:
            pass
        try:
            action_group.add_argument("-noffset", required=False, type=float,
                default=0, help="""Offset added to 15N dimension """)
        except argparse.ArgumentError:
            pass

        # Arguments inherited from superclass
        Dataset.construct_argparser(parser)

        return parser

    def __init__(self, hoffset=0, noffset=0, outfile=None, interactive=False,
            **kwargs):
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
        self.hsqc_df.index = pd.MultiIndex.from_product([hydrogen, nitrogen],
            names=["1H", "15N"])

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
        hydrogen = np.array(
            nmrglue.pipe.make_uc(parameters, intensity, dim=1).ppm_scale(),
            np.float32)
        nitrogen = np.array(
            nmrglue.pipe.make_uc(parameters, intensity, dim=0).ppm_scale(),
            np.float32)

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
        from ..myplotspec import multi_pop_merged

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


#################################### MAIN #####################################
if __name__ == "__main__":
    HSQCDataset.main()