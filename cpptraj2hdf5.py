#!/usr/bin/python
# -*- coding: utf-8 -*-
#   moldynplot.cpptraj2hdf5.py
#
#   Copyright (C) 2012-2016 Karl T Debiec
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
import h5py
################################## FUNCTIONS ##################################
def process_cpptraj(infile, outfile, address, dtype, scaleoffset,
    verbose=1, debug=0, **kwargs):
    """
    .. todo:
        - Support multiple infiles
    """
    from os import devnull
    from subprocess import Popen, PIPE

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
            print(data.shape)
            print(np.mean(data, axis=0))
            hdf5_file.create_dataset(address, data=data, dtype=dtype,
              chunks=True, compression="gzip", scaleoffset=scaleoffset)
            hdf5_file[address].attrs["fields"] = list(fields)

def process_saxs(package, infiles, outfile, address, dtype, scaleoffset,
    verbose=1, debug=0, **kwargs):
    """
    """
    from os.path import expandvars
    from glob import glob

    # Process arguments
    processed_infiles = []
    for infile in infiles:
        processed_infiles.extend(sorted(glob(expandvars(infile))))
    if len(processed_infiles) == 0:
        print("No infiles found matching '{0}', exiting".format(infiles))
        return
    infiles = processed_infiles
    if verbose >= 1:
        print("Loading SAXS data from {0} infiles, ".format(len(infiles)) +
              "starting with '{0}'".format(infiles[0]))
    outfile = expandvars(outfile)

    if package == "saxs_md":
        q = None
        data = None
        for i, infile in enumerate(infiles):
            if verbose >= 2:
                print("Loading SAXS data from {0}".format(infile))
            datum = np.loadtxt(infile, skiprows=3)
            if q is None:
                q = datum[:,0]
                if verbose >= 1:
                    print("q contains {0} points ".format(q.size) +
                          "ranging from {0} to {1} Å^-1".format(q[0], q[-1]))
                if verbose >= 2:
                    print("q:\n{0}".format(q))
            else:
                if not (datum[:,0] == q).all():
                    raise()
            if data is None:
                data = np.zeros((len(infiles), q.size))
            data[i,:] = datum[:,1]

    elif package == "foxs" or package == "crysol":
        q = None
        data = None
        j = 0
        for infile in infiles:
            if verbose >= 2:
                print("Loading SAXS data from {0}".format(infile))
            datum = np.loadtxt(infile, comments=["#","D"])
            if q is None:
                q = np.unique(datum[:,0])
                if verbose >= 1:
                    print("q contains {0} points ".format(q.size) +
                          "ranging from {0} to {1} Å^-1".format(q[0], q[-1]))
                if verbose >= 2:
                    print("q:\n{0}".format(q))
            else:
                if not (np.unique(datum[:,0]) == q).all():
                    raise()
            n_frames = datum.shape[0] / q.size
            if data is None:
                data = np.zeros((n_frames*len(infiles), q.size))
            intensity = np.reshape(datum[:,1], (n_frames, q.size))
            data[j:j+n_frames,:] = intensity
            j += n_frames

    if verbose >= 1:
        print("Loaded {0} intensity datasets".format(data.shape[0]))
    if verbose >= 2:
        print("Intensity:\n{0}".format(data))
        print(data.shape)

    # Open hdf5 file
    with h5py.File(outfile) as hdf5_file:
        if verbose >= 1:
            print("Writing q to '{0}[{1}/q]'".format(outfile, address))
        hdf5_file.create_dataset(address + "/q", data=q, dtype=dtype,
          chunks=True, compression="gzip", scaleoffset=scaleoffset)
        if verbose >= 1:
            print("Writing intensity to '{0}[{1}/intensity]'".format(outfile,
              address))
        hdf5_file.create_dataset(address + "/intensity", data=data,
          dtype=dtype, chunks=True, compression="gzip",
          scaleoffset=scaleoffset)

#################################### MAIN #####################################
if __name__ == "__main__":
    import argparse

    # Prepare argument parser
    parser            = argparse.ArgumentParser(
      description     = __doc__,
      formatter_class = argparse.RawTextHelpFormatter)

    kind_subparser = parser.add_subparsers(
      title           = "kind")

    # SAXS
    saxs_parser = kind_subparser.add_parser(
      name  = "saxs",
      help  = "small-angle X-ray scattering from 'saxs_md' or 'FOXS'")
    saxs_package = saxs_parser.add_mutually_exclusive_group()
    saxs_package.add_argument(
      "-s", "--saxs_md",
      action   = "store_const",
      const    = "saxs_md",
      default  = "saxs_md",
      dest     = "package",
      help     = "parse output from AmberTools' 'saxs_md' program (default)")
    saxs_package.add_argument(
      "-c", "--crysol",
      action   = "store_const",
      const    = "crysol",
      default  = "saxs_md",
      dest     = "package",
      help     = "parse output from 'crysol' program")
    saxs_package.add_argument(
      "-f", "--foxs",
      action   = "store_const",
      const    = "foxs",
      default  = "saxs_md",
      dest     = "package",
      help     = "parse output from 'FOXS' program")
    saxs_parser.add_argument(
      "infiles",
      nargs    = "+",
      type     = str,
      metavar  = "infile",
      help     = "file(s) from which to load calculated SAXS curves")
    saxs_parser.add_argument(
      "outfile",
      type     = str,
      help     = "HDF5 file to which to dataset will be output")
    saxs_parser.set_defaults(
      function    = process_saxs,
      address     = "saxs",
      dtype       = np.float32,
      scaleoffset = 3)

    # Cpptraj
    cpptraj_parser = argparse.ArgumentParser(add_help=False)
    cpptraj_parser.add_argument(
      "infile",
      type     = str,
      metavar  = "infile",
      help     = "file from which to load cpptraj output; "
                 "may be plain text or gzip")
    cpptraj_parser.add_argument(
      "outfile",
      type     = str,
      help     = "HDF5 file to which to dataset will be output")
    cpptraj_parser.set_defaults(
      function = process_cpptraj)

    kind_subparser.add_parser(
      name        = "dihedral",
      parents     = [cpptraj_parser],
      help        = "cpptraj's 'dihedral' or 'multidihedral'").set_defaults(
      address     = "dihedral",
      dtype       = np.float32,
      scaleoffset = 4)
    kind_subparser.add_parser(
      name        = "hbond",
      parents     = [cpptraj_parser],
      help        = "cpptraj's 'hbond'").set_defaults(
      address     = "hbond",
      dtype       = np.uint8,
      scaleoffset = 1)
    kind_subparser.add_parser(
      name        = "jcoupling",
      parents     = [cpptraj_parser],
      help        = "cpptraj's 'jcoupling'").set_defaults(
      address     = "jcoupling",
      dtype       = np.float32,
      scaleoffset = 3)
    kind_subparser.add_parser(
      name        = "natcon",
      parents     = [cpptraj_parser],
      help        = "native contacts, NOTE: actually inter-residue minimum "
                    "distances, not cpptraj's 'nativecontacts' command, which "
                    "does not appear to work").set_defaults(
      address     = "natcon",
      dtype       = np.float32,
      scaleoffset = 4)
    kind_subparser.add_parser(
      name        = "perresrmsd",
      parents     = [cpptraj_parser],
      help        = "cpptraj's 'rms' with 'perres' option").set_defaults(
      address     = "perresrmsd",
      dtype       = np.float32,
      scaleoffset = 4)
    kind_subparser.add_parser(
      name        = "secstruct",
      parents     = [cpptraj_parser],
      help        = "cpptraj's 'secstruct'").set_defaults(
      address     = "secstruct",
      dtype       = np.uint8,
      scaleoffset = 3)

    # Verbosity
    for p in kind_subparser.choices.values():
        p.add_argument(
          "-address",
          type     = str,
          help     = "Address within HDF5 file at which to output dataset: "
                     "(default: /%(default)s)")
        verbosity = p.add_mutually_exclusive_group()
        verbosity.add_argument(
          "-v", "--verbose",
          action   = "count",
          default  = 1,
          help     = "Enable verbose output, may be specified more than once")
        verbosity.add_argument(
          "-q", "--quiet",
          action   = "store_const",
          const    = 0,
          default  = 1,
          dest     = "verbose",
          help     = "Disable verbose output")

    # Parse arguments
    kwargs = vars(parser.parse_args())
    kwargs.pop("function")(**kwargs)
