#!/usr/bin/python
# -*- coding: utf-8 -*-
#   moldynplot.dataset.ChemicalShiftDataset.py
#
#   Copyright (C) 2015-2017 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Represents an NMR chemical shift data
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
class ChemicalShiftDataset(SequenceDataset):
    """
    Represents an NMR chemical shift data
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
            parser = parser_or_subparsers.add_parser(name="chemical_shift",
              description=help_message, help=help_message)
        elif parser is None:
            parser = argparse.ArgumentParser(description=help_message)

        # Defaults
        if parser.get_default("class_") is None:
            parser.set_defaults(class_=ChemicalShiftDataset)

        # Arguments unique to this class
        arg_groups = {ag.title: ag for ag in parser._action_groups}

        # Input arguments
        input_group = arg_groups.get("input",
          parser.add_argument_group("input"))
        try:
            input_group.add_argument("-delays", dest="delays", metavar="DELAY",
              nargs="+", type=float, help="""delays for each infile,
                if infiles represent a
                         series; number of delays must match number of
                         infiles""")
        except argparse.ArgumentError:
            pass

        # Action arguments
        action_group = arg_groups.get("action",
          parser.add_argument_group("action"))
        try:
            action_group.add_argument("-relax", dest="calc_relax", type=str,
              nargs="?", default=None, const="r1", help="""Calculate
                relaxation rates and standard errors; may
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
            kwargs["infile_column_prefixes"] = ["{0:3.1f} ms".format(delay) for
                delay in delays]
        self.sequence_df = self.read(**kwargs)

        # Cut data
        if "use_indexes" in kwargs:
            use_indexes = np.array(kwargs.pop("use_indexes"))
            res_index = np.array(
              [int(i.split(":")[1]) for i in self.sequence_df.index.values])
            self.sequence_df = self.sequence_df[
                np.in1d(res_index, use_indexes)]

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
            return "{0}:{1}".format(name[-4:-1].upper(), name[2:-4])

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
                    header = " ".join(
                      Popen("head -n 1 {0}".format(infile), stdout=PIPE,
                        stderr=fnull, shell=True).stdout.read().strip().split(
                        "\t"))
                ccpnmr_header = sformat("""Number # Position F1 Position F2
                    Assign F1 Assign F2 Height Volume Line Width F1 (Hz) Line
                    Width F2 (Hz) Merit Details Fit Method Vol. Method""")
                if (header == ccpnmr_header):
                    read_csv_kw = dict(index_col=None, delimiter="\t",
                      dtype={"Position F1": np.float32,
                          "Position F2": np.float32, "Assign F1": np.str,
                          "Height": np.float32, "Volume": np.float32},
                      usecols=[str("Position F1"), str("Position F2"),
                          str("Assign F1"), str("Height"), str("Volume")],
                      converters={"Assign F1": convert_name})
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
            df = df.loc[
                sorted(df.index.values, key=lambda x: int(x.split(":")[1]))]
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
        df = kwargs.get("df")
        if df is None:
            if hasattr(self, "sequence_df"):
                df = self.sequence_df
            else:
                raise ()
        relax_kw = kwargs.get("relax_kw", {})
        kind = relax_kw.get("kind", "r1")
        intensity_method = relax_kw.get("intensity_method", "height")
        error_method = relax_kw.get("error_method", "mae")
        n_synth_datasets = relax_kw.get("n_synth_datasets", 1000)

        # Calculate relaxation rates
        re_column = re.compile(
          "^(?P<delay>\d+\.?\d*?) ms {0}".format(intensity_method))
        columns = [c for c in df.columns.values if re_column.match(c)]
        delays = np.array(
          [re.match(re_column, c).groupdict()["delay"] for c in columns],
          np.float) / 1000

        def calc_relax_rate(residue, **kwargs):
            """
            """
            from . import multiprocess_map

            if verbose >= 1:
                wiprint(
                  """Calculating {0} relaxation rate for {1}""".format(kind,
                    residue.name))

            def model_function(delay, intensity, relaxation):
                return intensity * np.exp(-1 * delay * relaxation)

            I = np.array(residue.filter(columns, np.float64))
            I0, R = curve_fit(model_function, delays, I, p0=(I[0], 1.0))[0]

            # Calculate error
            if error_method == "rmse":
                error = np.sqrt(
                  np.mean((I - model_function(delays, I0, R)) ** 2))
            elif error_method == "mae":
                error = np.mean(
                  np.sqrt((I - model_function(delays, I0, R)) ** 2))

            # Construct synthetic relaxation profiles
            synth_datasets = np.zeros((n_synth_datasets, I.size))
            for i, I_mean in enumerate(model_function(delays, I0, R)):
                synth_datasets[:, i] = np.random.normal(I_mean, error,
                  n_synth_datasets)

            def synth_fit_decay(synth_intensity):
                try:
                    synth_I0, synth_R = \
                        curve_fit(model_function, delays, synth_intensity,
                          p0=(I0, R))[0]
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


#################################### MAIN #####################################
if __name__ == "__main__":
    ChemicalShiftDataset.main()