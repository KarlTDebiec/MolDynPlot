#!/usr/bin/python
# -*- coding: utf-8 -*-
#   moldynplot.dataset.SequenceDataset.py
#
#   Copyright (C) 2015-2017 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Represents data that is a function of amino acid sequence
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
class SequenceDataset(Dataset):
    """
    Represents data that is a function of amino acid sequence

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
    default_hdf5_kw = dict(chunks=True, compression="gzip", dtype=np.float32,
      scaleoffset=5)

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
            parser = parser_or_subparsers.add_parser(name="sequence",
              description=help_message, help=help_message)
        elif parser_or_subparsers is None:
            parser = argparse.ArgumentParser(description=help_message)

        # Defaults
        if parser.get_default("class_") is None:
            parser.set_defaults(class_=SequenceDataset)

        # Arguments unique to this class
        arg_groups = {ag.title: ag for ag in parser._action_groups}

        # Input arguments
        input_group = arg_groups.get("input",
          parser.add_argument_group("input"))
        try:
            input_group.add_argument("-indexfile", required=False, type=str,
              help="""text file from which to load residue names; should
                         list amino acids in the form 'XAA:#' separated by
                         whitespace; if omitted will be taken from rows of
                         first infile; may contain environment variables""")
        except argparse.ArgumentError:
            pass

        # Arguments inherited from superclass
        Dataset.construct_argparser(parser)

        return parser

    @classmethod
    def get_cache_key(class_, **kwargs):
        """
        Generates key for dataset cache.

        See :class:`SequenceDataset<moldynplot.Dataset.SequenceDataset>`
        for argument details.

        Returns:
          tuple: Cache key; contains arguments sufficient to reconstruct
          dataset
        """
        from ..myplotspec import multi_pop_merged

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
        return (class_, tuple(infiles), use_indexes, tuple(read_csv_kw))

    def __init__(self, calc_pdist=False, outfile=None, interactive=False,
      **kwargs):
        """
        Arguments:
          infile{s} (list): Path(s) to input file(s); may contain
            environment variables and wildcards
          use_indexes (list): Residue indexes to select from DataFrame,
            once DataFrame has already been loaded
          calc_pdist (bool): Calculate probability distribution
            using :meth:`calc_pdist`
          pdist_kw (dict): Keyword arguments used to configure
            probability distribution calculation
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

        # Process data
        if "use_indexes" in kwargs:
            use_indexes = np.array(kwargs.pop("use_indexes"), np.int)
            res_index = np.array(
              [int(i.split(":")[1]) for i in self.sequence_df.index.values])
            self.sequence_df = self.sequence_df[
                np.in1d(res_index, use_indexes)]

        # Output data
        if verbose >= 2:
            wiprint("Processed sequence DataFrame:")
            print(self.sequence_df)
        if outfile is not None:
            self.write(df=self.sequence_df, outfile=outfile, **kwargs)

        # Calculate probability distribution
        if calc_pdist:
            pdist_kw = kwargs.get("pdist_kw", {})
            self.pdist_df = self.calc_pdist(df=self.sequence_df,
              verbose=verbose, **pdist_kw)
            if verbose >= 2:
                print("Processed pdist DataFrame:")
                print(self.pdist_df)

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
        elif (sum(
          [1 for a in [re_res.match(str(b)) for b in df.index.values] if
                  a is not None]) == len(df.index.values)):
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
            df = df.loc[
                sorted(df.index.values, key=lambda x: int(x.split(":")[1]))]
        else:
            df = df.loc[sorted(df.index.values)]

        return df


#################################### MAIN #####################################
if __name__ == "__main__":
    SequenceDataset.main()