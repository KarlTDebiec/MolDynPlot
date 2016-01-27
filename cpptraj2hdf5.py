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
      choices  = [str("dihedral"), str("hbond"), str("secstruct")],
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
    kwargs      = vars(parser.parse_args())
    infile      = kwargs["infile"]
    outfile     = kwargs["outfile"]
    
    if kwargs["kind"] == "dihedral":
        field_dtype = "f4"
        address     = "dihedral"
    elif kwargs["kind"] == "hbond":
        field_dtype = "i1"
        address     = "hbond"
    elif kwargs["kind"] == "secstruct":
        field_dtype = "i1"
        address     = "secstruct"

    # Load input file
    with open(devnull, "w") as fnull:

        # Determine field names and prepare dtype
        command = "head -n 1 {0}".format(infile)
        process = Popen(command, stdout=PIPE, stderr=fnull, shell=True)
        fields = process.stdout.read().split()
        if fields.pop(0) != "#Frame":
            raise TypeError("cpptraj2hdf5 currently only supports input in " +
                            "the form of '#Frame field_1 field_2 ...'")
        dtype = [(str("Frame"), np.int64)]
        for field in fields:
            dtype += [(str(field), field_dtype)]

    # Open hdf5 file
    with h5py.File(outfile) as hdf5_file:

        # Just use loadtxt for now
        data = np.loadtxt(infile, dtype=dtype)

        hdf5_file.create_dataset(address, data=data,
          chunks=True, compression="gzip", maxshape=(None,))


    # May need to use code below for large files that cannot fit in memory
    # Very slow; there is probably a much better solution
#            # Instantiate resizable dataset
#            hdf5_file.create_dataset(address, data=np.empty(1, dtype),
#              chunks=True, compression="gzip", maxshape=(None,))
#
#            first_line = 2
#            first_index = 0
#            while True:
#                last_line = first_line + chunk_size - 1
#
#                # Load frame indexes for chunk
#                command = "sed -n '{0},{1}p' {2} ".format(first_line, last_line, infile) + \
#                          "| cut -c1-8 | tr -d '\n'"
#                process = Popen(command, stdout=PIPE, stderr=fnull, shell=True)
#                input_bytes = bytearray(process.stdout.read())
#                frame = np.array(np.frombuffer(input_bytes, dtype="S8", count=int((len(input_bytes)) / 8)), np.int64)
#                last_index = first_index + frame.size
#
#                # Load data for chunk
#                command = "sed -n '{0},{1}p' {2} | cut -c9- | tr -d '\n'".format(first_line, last_line, infile)
#                process = Popen(command, stdout=PIPE, stderr=fnull, shell=True)
#                input_bytes = bytearray(process.stdout.read())
#                data = np.array(np.frombuffer(input_bytes, dtype="S13", count=int((len(input_bytes)) / 13)), np.int8)
#
#                # Store data in hdf5 file
#                hdf5_file[address].resize(size=(last_index,))
#                hdf5_file[address][first_index:last_index] = np.rec.array(frame, dtype=[dtype[0]])
#                hdf5_file[address][first_index:last_index] = np.rec.array(data,  dtype=dtype[1:])
#
#                first_line  += chunk_size
#                first_index += chunk_size
#                if frame.size < chunk_size:
#                    break
