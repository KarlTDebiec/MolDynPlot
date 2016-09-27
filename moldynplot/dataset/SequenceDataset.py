#!/usr/bin/python
# -*- coding: utf-8 -*-
#   moldynplot.dataset.SequenceDataset.py
#
#   Copyright (C) 2015-2016 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Processes data that is a function of amino acid sequence

.. todo:
  - Fix separation and ordering of argument groups: input, action, output
  - Move relaxation error here
  - Move relaxation heteronuclear noe here
  - Move relaxation pre ratio here
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
        help_message = """Process standard data"""
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
            res_index = np.array([int(i.split(":")[1])
              for i in self.sequence_df.index.values])
            self.sequence_df = self.sequence_df[np.in1d(res_index,use_indexes)]

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


class ChemicalShiftDataset(SequenceDataset):
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
        help_message = """Process NMR chemical shift data"""
        if isinstance(parser_or_subparsers, argparse.ArgumentParser):
            parser = parser_or_subparsers
        elif isinstance(parser_or_subparsers, argparse._SubParsersAction):
            parser = parser_or_subparsers.add_parser(
              name        = "chemical_shift",
              description = help_message,
              help        = help_message)
        elif parser is None:
            parser = argparse.ArgumentParser(
              description = help_message)

        # Defaults
        if parser.get_default("cls") is None:
            parser.set_defaults(cls=ChemicalShiftDataset)

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
        from ..myplotspec import multi_pop_merged

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


class RelaxDataset(SequenceDataset):
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
        help_message = """Process NMR relaxation data"""
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


class IREDDataset(RelaxDataset):
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
        help_message = """Process NMR relaxation data calculated from MD
                       simulation using the iRED method as implemented in
                       cpptraj"""
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
            parser.set_defaults(cls=IREDDataset)

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
        kind = IREDDataset._identify_infile(infile)

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
            if not (("r1" in columns and "r2" in columns and "noe" in columns)
            or      ("s2" in columns)):
                raise Exception(sformat("""DataFrame loaded from '{0}' does not
                  appear to contain either relaxation ('r1', 'r2', 'noe') or
                  order parameter ('s2') columns""".format(infile)))

        # Average, if applicable
        df = self.average_independent(relax_dfs, order_dfs)

        return df


#################################### MAIN #####################################
def main():
    import argparse

    # Prepare argument parser
    parser = argparse.ArgumentParser(
      description = __doc__)
    subparsers = parser.add_subparsers(
      dest        = "mode",
      description = "")

    SequenceDataset.construct_argparser(subparsers)
    ChemicalShiftDataset.construct_argparser(subparsers)
    RelaxDataset.construct_argparser(subparsers)
    IREDDataset.construct_argparser(subparsers)

    kwargs  = vars(parser.parse_args())
    kwargs.pop("cls")(**kwargs)

if __name__ == "__main__":
    main()
