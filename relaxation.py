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
      "-infile",
      required = True,
      nargs    = "+",
      type     = str,
      help     = "cpptraj output file from which to load dataset; " +
                 "may be plain text or gzip")
    parser.add_argument(
      "-outfile",
      required = True,
      type     = str,
      help     = "Text or hdf5 file to which to dataset will be output")

    # Parse arguments
    kwargs  = vars(parser.parse_args())
    infiles  = kwargs["infile"]
    outfile = kwargs["outfile"]
    freqs_MHz = np.squeeze(np.array(kwargs["frequency"]))
    freqs_rad_ps = freqs_MHz*1e6*1e-12*2*np.pi
    fields = None
    data = None

    for i, infile in enumerate(infiles):
        print(infile)

        # Load column names from first infile (residue names and numbers)
        if fields is None:
            with open(devnull, "w") as fnull:
                command = "head -n 1 {0}".format(infile)
                process = Popen(command, stdout=PIPE, stderr=fnull, shell=True)
                fields = process.stdout.read().split()
                if fields.pop(0) != "#Frame":
                    pass

        # Read input file(s) containing autocorrelation function from cpptraj
        acf = np.loadtxt(infile, dtype=np.float32)[:,1:]
        
        # Initialize data, now that output is known
        if data is None:
            data = np.zeros((len(infiles), freqs_rad_ps.size, acf.shape[1]),
                     np.float64)

        # Calculate fft of autocorrelation function to obtain spectral density
        sdf = np.real(np.fft.rfft(acf, axis=0))

        # Calculate appropriate frequencies for desired spectral density
        sdf_freq = np.fft.fftfreq(sdf.shape[0])

        # Interpolate calculated spectral density 
        for k in range(sdf.shape[1]):
            interpolated_sdf = interp1d(sdf_freq, sdf[:,k])
            data[i,:,k] = interpolated_sdf(freqs_rad_ps)
    mean = np.mean(data, axis=0)
    std = np.std(data, axis=0)

    with open(outfile, "w") as output:
        output.write("#           ")
        for field in fields:
            output.write("{0:>12}".format(field))
        output.write("\n")
        for i, freq_MHz in enumerate(freqs_MHz):
            output.write("{0:<12}".format("J({0})".format(freq_MHz)))
            for j in range(mean.shape[1]):
                output.write("{0:>12.2f}".format(mean[i,j]))
            output.write("\n")
            output.write("{0:<12}".format("J({0}) se".format(freq_MHz)))
            for j in range(mean.shape[1]):
                output.write("{0:>12.2f}".format(std[i,j]))
            output.write("\n")
