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
################################### MODULES ###################################
from __future__ import absolute_import,division,print_function,unicode_literals
if __name__ == "__main__":
    __package__ = str("moldynplot")
    import moldynplot
import numpy as np
################################## FUNCTIONS ##################################
def saxs(address="saxs", dtype=np.float32, scaleoffset=3, **kwargs):

    infiles = kwargs["infile"]
    outfile = kwargs["outfile"]
    q = None
    data = None
    for i, infile in enumerate(infiles):
        datum = np.loadtxt(infile, skiprows=3)
        if q is None:
            q = datum[:,0]
        else:
            if not (datum[:,0] == q).all():
                raise()
        if data is None:
            data = np.zeros((len(infiles), q.size))
        data[i,:] = datum[:,1]

    # Open hdf5 file
    with h5py.File(outfile) as hdf5_file:

        print(data, data.shape)
        hdf5_file.create_dataset(address + "/q", data=q, dtype=dtype,
          chunks=True, compression="gzip", scaleoffset=scaleoffset)
        hdf5_file.create_dataset(address + "/intensity", data=data,
          dtype=dtype, chunks=True, compression="gzip",
          scaleoffset=scaleoffset)

#################################### MAIN #####################################
if __name__ == "__main__":
    import argparse
    from os import devnull
    from subprocess import Popen, PIPE
    import h5py

    # Prepare argument parser
    parser            = argparse.ArgumentParser(
      description     = __doc__,
      formatter_class = argparse.RawTextHelpFormatter)
    parser.add_argument(
      "kind",
      type     = str,
      choices  = [str("dihedral"), str("hbond"),      str("jcoupling"),
                  str("natcon"),   str("perresrmsd"), str("saxs"),
                  str("secstruct")],
      help     = "kind of dataset")
    parser.add_argument(
      "infile",
      nargs    = "+",
      type     = str,
      help     = "cpptraj output file from which to load dataset; " +
                 "may be plain text or gzip")
    parser.add_argument(
      "outfile",
      type     = str,
      help     = "HDF5 file to which to dataset will be output")

    # Parse arguments
    kwargs  = vars(parser.parse_args())
    if kwargs["kind"] == "dihedral":
        dtype = np.float32
        scaleoffset = 4
        address = "dihedral"
    elif kwargs["kind"] == "hbond":
        dtype = np.uint8
        scaleoffset = 1
        address = "hbond"
    elif kwargs["kind"] == "jcoupling":
        dtype = np.float32
        scaleoffset = 3
        address = "jcoupling"
    elif kwargs["kind"] == "natcon":
        dtype = np.float32
        scaleoffset = 4
        address = "natcon"
    elif kwargs["kind"] == "perresrmsd":
        dtype = np.float32
        scaleoffset = 4
        address = "perresrmsd"
    elif kwargs["kind"] == "saxs":
        from sys import exit
        saxs(**kwargs)
        exit()
    elif kwargs["kind"] == "secstruct":
        dtype = np.uint8
        scaleoffset = 3
        address = "secstruct"
    infile  = kwargs["infile"][0]
    outfile = kwargs["outfile"]

    if infile.endswith("gnu"):
        raw = np.genfromtxt(infile, skip_header=13,
                  skip_footer=2,invalid_raise=False,
                  dtype=np.float64)
        dim_1 = int(np.max(raw[:,0]))
        dim_2 = int(np.max(raw[:,1]))
        data = np.zeros((int(np.max(raw[:,0])), int(np.max(raw[:,1]))), dtype)
        for frame in raw:
            data[int(frame[0]) - 1, int(frame[1]) - 1] = dtype(frame[2])

        # Open hdf5 file
        with h5py.File(outfile) as hdf5_file:
            hdf5_file.create_dataset(address, data=data, dtype=dtype,
              chunks=True, compression="gzip", scaleoffset=scaleoffset)

    else:
        # Determine dimensions of dataset
        with open(devnull, "w") as fnull:

            # Determine field names and prepare dtype
            command = "head -n 1 {0}".format(infile)
            process = Popen(command, stdout=PIPE, stderr=fnull, shell=True)
            fields = process.stdout.read().split()
            if fields.pop(0) != "#Frame":
                raise TypeError("cpptraj2hdf5 currently only supports input " +
                                "in the form of '#Frame field_1 field_2 ...'")
            n_fields = len(fields)

        # Open hdf5 file
        with h5py.File(outfile) as hdf5_file:
            def iter_func():
                with open(infile, "r") as open_file:
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
            hdf5_file[address].attrs["fields"] = list(fields)
