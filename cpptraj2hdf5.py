#!/usr/bin/python
# -*- coding: utf-8 -*-
#   moldynplot.cpptraj2hdf5.py
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
    import h5py
    import numpy as np

    # Prepare argument parser
    parser            = argparse.ArgumentParser(
      description     = __doc__,
      formatter_class = argparse.RawTextHelpFormatter)
    parser.add_argument(
      "kind",
      type     = str,
      choices  = [str("dihedral"),   str("hbond"), 
                  str("perresrmsd"), str("secstruct")],
      help     = "kind of dataset")
    parser.add_argument(
      "infile",
      type     = str,
      help     = "cpptraj output file from which to load dataset; " +
                 "may be plain text or gzip")
    parser.add_argument(
      "outfile",
      type     = str,
      help     = "HDF5 file to which to dataset will be output")

    # Parse arguments
    kwargs  = vars(parser.parse_args())
    infile  = kwargs["infile"]
    outfile = kwargs["outfile"]
    if kwargs["kind"] == "dihedral":
        dtype = np.float32
        scaleoffset = 4
        address = "dihedral"
    elif kwargs["kind"] == "hbond":
        dtype = np.uint8
        scaleoffset = 1
        address = "hbond"
    elif kwargs["kind"] == "perresrmsd":
        dtype = np.float32
        scaleoffset = 4
        address = "perresrmsd"
    elif kwargs["kind"] == "secstruct":
        dtype = np.float32
        scaleoffset = 3
        address = "secstruct"

    # Determine dimensions of dataset
    with open(devnull, "w") as fnull:

        # Determine field names and prepare dtype
        command = "head -n 1 {0}".format(infile)
        process = Popen(command, stdout=PIPE, stderr=fnull, shell=True)
        fields = process.stdout.read().split()
        if fields.pop(0) != "#Frame":
            raise TypeError("cpptraj2hdf5 currently only supports input in " +
                            "the form of '#Frame field_1 field_2 ...'")
        n_fields = len(fields)

    # Open hdf5 file
    with h5py.File(outfile) as hdf5_file:
        # Instantiate resizable dataset

        def iter_func():
            with open(infile, 'r') as open_file:
                next(open_file)
                for line in open_file:
                    line = line.split()[1:]
                    for item in line:
                        yield dtype(item)

        data = np.fromiter(iter_func(), dtype=dtype)
        data = data.reshape((-1, n_fields))
        print(data, data.shape)
        hdf5_file.create_dataset(address, data=data, dtype=dtype,
          chunks=True, compression="gzip", scaleoffset=scaleoffset)
