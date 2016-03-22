#!/usr/bin/python
# -*- coding: utf-8 -*-
#   moldynplot.relaxation.py
#
#   Copyright (C) 2012-2015 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Processes NMR relaxation-related data output by cpptraj
"""
################################### MODULES ####################################
from __future__ import absolute_import,division,print_function,unicode_literals
##################################### MAIN #####################################
#def process_acf():
#    from os import devnull
#    from subprocess import Popen, PIPE
#    from scipy.interpolate import interp1d
#    import numpy as np
#    headerfile = kwargs["headerfile"]
#    infiles  = kwargs["infile"]
#    interpolation = kwargs["interpolation"]
#    outfile = kwargs["outfile"]
#    freqs_MHz = np.squeeze(np.array(kwargs["frequency"]))
#    freqs_rad_ps = freqs_MHz*1e6*1e-12*2*np.pi
#    if headerfile is None:
#        headerfile = infiles[0]
#    header = None
#    data = None
#
#    for i, infile in enumerate(infiles):
#        print("Loading from '{0}'".format(infile))
#
#        # Load column names from first infile (residue names and numbers)
#        if header is None:
#            with open(devnull, "w") as fnull:
#                command = "head -n 1 {0}".format(headerfile)
#                process = Popen(command, stdout=PIPE, stderr=fnull, shell=True)
#            fields = process.stdout.read().split()
#            if fields[0] in ["#", "#Frame"]:
#                fields.pop(0)
#            header = "#           "
#            for field in fields:
#                header += "{0:>12}".format(field)
#            header += "\n"
#
#        # Read input file(s) containing autocorrelation function from cpptraj
#        raw_data = np.loadtxt(infile, dtype=np.float32)
#        acf = raw_data[:,1:]
#        dt = np.mean(raw_data[1:,0] - raw_data[:-1,0])
#        
#        # Initialize data, now that output is known
#        if data is None:
#            data = np.zeros((len(infiles), freqs_rad_ps.size, acf.shape[1]),
#                     np.float64)
#
#        # Calculate fft of autocorrelation function to obtain spectral density
#        sdf = np.real(np.fft.rfft(acf, axis=0))
#
#        # Calculate appropriate frequencies for desired spectral density
#        sdf_freq = np.fft.fftfreq(sdf.shape[0], d=1/dt)
#
#        # Interpolate calculated spectral density 
#        for k in range(sdf.shape[1]):
#            max_index = np.abs(sdf_freq - freqs_rad_ps.max()).argmin() + 10
#            interpolated_sdf = interp1d(sdf_freq[:max_index],
#              sdf[:max_index,k], kind=interpolation)
#            data[i,:,k] = interpolated_sdf(freqs_rad_ps)
#    mean = np.mean(data, axis=0)
#    std = np.std(data, axis=0)
#
#    with open(outfile, "w") as output:
#        print("Writing results to '{0}'".format(outfile))
#        output.write(header)
#        for i, freq_MHz in enumerate(freqs_MHz):
#            output.write("{0:<12}".format("J({0})".format(freq_MHz)))
#            for j in range(mean.shape[1]):
#                output.write("{0:>12.2f}".format(mean[i,j]))
#            output.write("\n")
#            output.write("{0:<12}".format("J({0})se".format(freq_MHz)))
#            for j in range(mean.shape[1]):
#                output.write("{0:>12.2f}".format(std[i,j]))
#            output.write("\n")

def process_ired(infiles, outfile, indexfile=None, **kwargs):
    from os import devnull
    import re
    from subprocess import Popen, PIPE
    import pandas as pd
    pd.set_option("display.width", None)
    import numpy as np

    r1r2noe_datasets = []
    s2_datasets = []

    # Load data
    for i, infile in enumerate(infiles):
        with open(devnull, "w") as fnull:
            fields = Popen("head -n 1 {0}".format(infile),
              stdout=PIPE, stderr=fnull, shell=True).stdout.read().strip()
        re_t1t2noe = re.compile(
          "^#Vec\s+[\w_]+\[T1\]\s+[\w_]+\[\T2\]\s+[\w_]+\[NOE\]$")
        re_s2 = re.compile("^#Vec\s+[\w_]+\[S2\]$")
        if re.match(re_t1t2noe, fields):
            raw_data = np.loadtxt(infile, dtype=np.float32)
            read_csv_kw = kwargs.get("read_csv_kw", dict(delim_whitespace=True,
              header=0, index_col=0, names=["r1","r2","noe"]))
            raw_data = pd.read_csv(infile, **read_csv_kw)
            raw_data["r1"] = 1 / raw_data["r1"]
            raw_data["r2"] = 1 / raw_data["r2"]
            r1r2noe_datasets.append(raw_data)
        elif re.match(re_s2, fields):
            raw_data = np.loadtxt(infile, dtype=np.float32)
            read_csv_kw = kwargs.get("read_csv_kw", dict(delim_whitespace=True,
              header=0, index_col=0, names=["s2"]))
            raw_data = pd.read_csv(infile, **read_csv_kw)
            s2_datasets.append(raw_data)
        else:
            raise Exception()
    if indexfile is not None:
        residue = np.loadtxt(indexfile, dtype=np.str).flatten()

    # Process data
    items = []
    fmt = []
    if indexfile is not None:
        items.append(("residue", residue))
        fmt.append("%12s")
    else:
        fmt.append("%12d")
    if   len(r1r2noe_datasets) >= 2:
        r1r2noe_mean = pd.concat(r1r2noe_datasets).groupby(level=0).mean()
        r1r2noe_std  = pd.concat(r1r2noe_datasets).groupby(level=0).std()
        items.extend(
            [("r1",  r1r2noe_mean["r1"]),  ("r1_se",  r1r2noe_std["r1"]),
             ("r2",  r1r2noe_mean["r2"]),  ("r2_se",  r1r2noe_std["r2"]),
             ("noe", r1r2noe_mean["noe"]), ("noe_se", r1r2noe_std["noe"])])
        fmt.extend(["%11.5f", "%11.5f", "%11.5f",
                    "%11.5f", "%11.5f", "%11.5f"])
    elif len(r1r2noe_datasets) == 1:
        r1r2noe_mean = r1r2noe_datasets[0]
        items.extend(
            [("r1",  r1r2noe_mean["r1"]),
             ("r2",  r1r2noe_mean["r2"]),
             ("noe", r1r2noe_mean["noe"])])
        fmt.extend(["%11.5f", "%11.5f", "%11.5f"])
    if   len(s2_datasets) >= 2:
        s2_mean      = pd.concat(s2_datasets).groupby(level=0).mean()
        s2_std       = pd.concat(s2_datasets).groupby(level=0).std()
        items.extend(
            [("s2",  s2_mean["s2"]),       ("s2_se",  s2_std["s2"])])
        fmt.extend(["%11.5f", "%11.5f"])
    elif len(s2_datasets) == 1:
        s2_mean      = s2_datasets[0]
        items.extend(
            [("s2",  s2_mean["s2"])])
        fmt.extend(["%11.5f"])

    data = pd.DataFrame.from_items(items)
    if indexfile is not None:
        data.set_index("residue", inplace=True)
    else:
        data.index.name = "vector"
    columns = [data.index.name] + list(data.columns.values)
    header = "{0:<10s}".format(columns.pop(0))
    for column in columns:
        header += "{0:>12s}".format(column)

    np.savetxt(outfile, np.column_stack((data.index.values, data.values)),
      fmt=fmt, header=header)

if __name__ == "__main__":
    import argparse

    # Prepare argument parser
    parser            = argparse.ArgumentParser(
      description     = __doc__,
      formatter_class = argparse.RawTextHelpFormatter)
    subparsers = parser.add_subparsers(
                   dest            = "mode",
                   description     = "")

    # Prepare autocorrelation function subparser
    acf_subparser  = subparsers.add_parser(
      name     = "acf",
      help     = "Process N-H vector autocorrelation function data")
    input_group  = acf_subparser.add_argument_group("input")
    action_group = acf_subparser.add_argument_group("action")
    output_group = acf_subparser.add_argument_group("output")

    input_group.add_argument(
      "-infile",
      required = True,
      nargs    = "+",
      type     = str,
      help     = "cpptraj output file(s) from which to load datasets; " +
                 "may be plain text or compressed")
    input_group.add_argument(
      "-headerfile",
      required = False,
      type     = str,
      help     = "Text file from which to load column names; if omitted " +
                 "will be taken from columns of first infile")
    action_group.add_argument(
      "-frequency",
      required = True,
      nargs    = "+",
      type     = int,
      action   = "append",
      help     = "Frequency(ies) at which to calculate spectral density (MHz)")
    action_group.add_argument(
      "-interpolation",
      required = False,
      type     = str,
      choices  = ["linear", "nearest", "zero",
                  "slinear", "quadratic", "cubic"],
      default  = "cubic",
      help     = "Method of interpolating spectral density function; " +
                 "passed to scipy.interpolate.interp1d")
    output_group.add_argument(
      "-outfile",
      required = True,
      type     = str,
      help     = "Text file to which processed data will be output")

    # Prepare iRED subparser
    ired_subparser  = subparsers.add_parser(
      name     = "ired",
      help     = "Process iRED data")
    input_group  = ired_subparser.add_argument_group("input")
    action_group = ired_subparser.add_argument_group("action")
    output_group = ired_subparser.add_argument_group("output")
    input_group.add_argument(
      "-infile",
      required = True,
      dest     = "infiles",
      nargs    = "+",
      type     = str,
      help     = "cpptraj output file(s) from which to load datasets; " +
                 "may be plain text or compressed")
    input_group.add_argument(
      "-indexfile",
      required = False,
      type     = str,
      help     = "Text file from which to load residue names; if omitted " +
                 "will be taken from columns of first infile")
    output_group.add_argument(
      "-outfile",
      required = True,
      type     = str,
      help     = "Text file to which processed data will be output")

    # Parse arguments and run selected function
    kwargs  = vars(parser.parse_args())
    if kwargs.pop("mode") == "ired":
        process_ired(**kwargs)
