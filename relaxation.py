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
"""
################################### MODULES ####################################
from __future__ import absolute_import,division,print_function,unicode_literals
if __name__ == "__main__":
    __package__ = str("moldynplot")
    import moldynplot
##################################### MAIN #####################################
if __name__ == "__main__":
    import argparse
    from os import devnull
    from subprocess import Popen, PIPE
    from sys import exit
    from scipy.interpolate import interp1d
    import numpy as np

    # Prepare argument parser
    parser            = argparse.ArgumentParser(
      description     = __doc__,
      formatter_class = argparse.RawTextHelpFormatter)
    parser.add_argument(
      "-frequency",
      required = True,
      nargs    = "+",
      type     = int,
      action   = "append",
      help     = "Frequency at which to calculate spectral density (MHz)")
    parser.add_argument(
      "-interpolation",
      required = False,
      type     = str,
      choices  = ["linear", "nearest", "zero",
                  "slinear", "quadratic", "cubic"],
      default  = "cubic",
      help     = "Method of interpolating spectral density function; " +
                 "passed to scipy.interpolate.interp1d")
    parser.add_argument(
      "-headerfile",
      required = False,
      type     = str,
      help     = "Text file from which to load column names; by default " +
                 "will be taken from first infile")
    parser.add_argument(
      "-infile",
      required = True,
      nargs    = "+",
      type     = str,
      help     = "cpptraj output file from which to load datasets; " +
                 "may be plain text or gzip")
    parser.add_argument(
      "-outfile",
      required = True,
      type     = str,
      help     = "Text file to which to dataset will be output")

    # Parse arguments
    kwargs  = vars(parser.parse_args())
    headerfile = kwargs["headerfile"]
    infiles  = kwargs["infile"]
    interpolation = kwargs["interpolation"]
    outfile = kwargs["outfile"]
    freqs_MHz = np.squeeze(np.array(kwargs["frequency"]))
    freqs_rad_ps = freqs_MHz*1e6*1e-12*2*np.pi
    if headerfile is None:
        headerfile = infiles[0]
    header = None
    data = None

    for i, infile in enumerate(infiles):
        print("Loading from '{0}'".format(infile))

        # Load column names from first infile (residue names and numbers)
        if header is None:
            with open(devnull, "w") as fnull:
                command = "head -n 1 {0}".format(headerfile)
                process = Popen(command, stdout=PIPE, stderr=fnull, shell=True)
            fields = process.stdout.read().split()
            if fields[0] in ["#", "#Frame"]:
                fields.pop(0)
            header = "#           "
            for field in fields:
                header += "{0:>12}".format(field)
            header += "\n"

        # Read input file(s) containing autocorrelation function from cpptraj
        raw_data = np.loadtxt(infile, dtype=np.float32)
        acf = raw_data[:,1:]
        dt = np.mean(raw_data[1:,0] - raw_data[:-1,0])
        
        # Initialize data, now that output is known
        if data is None:
            data = np.zeros((len(infiles), freqs_rad_ps.size, acf.shape[1]),
                     np.float64)

        # Calculate fft of autocorrelation function to obtain spectral density
        sdf = np.real(np.fft.rfft(acf, axis=0))

        # Calculate appropriate frequencies for desired spectral density
        sdf_freq = np.fft.fftfreq(sdf.shape[0], d=1/dt)

        # Interpolate calculated spectral density 
        for k in range(sdf.shape[1]):
            max_index = np.abs(sdf_freq - freqs_rad_ps.max()).argmin() + 10
            interpolated_sdf = interp1d(sdf_freq[:max_index],
              sdf[:max_index,k], kind=interpolation)
            data[i,:,k] = interpolated_sdf(freqs_rad_ps)
    mean = np.mean(data, axis=0)
    std = np.std(data, axis=0)

    with open(outfile, "w") as output:
        print("Writing results to '{0}'".format(outfile))
        output.write(header)
        for i, freq_MHz in enumerate(freqs_MHz):
            output.write("{0:<12}".format("J({0})".format(freq_MHz)))
            for j in range(mean.shape[1]):
                output.write("{0:>12.2f}".format(mean[i,j]))
            output.write("\n")
            output.write("{0:<12}".format("J({0})se".format(freq_MHz)))
            for j in range(mean.shape[1]):
                output.write("{0:>12.2f}".format(std[i,j]))
            output.write("\n")
