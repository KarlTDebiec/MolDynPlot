#!/usr/bin/python
# -*- coding: utf-8 -*-
#   moldynplot.map2pdb.py
#
#   Copyright (C) 2012-2017 Karl T Debiec
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
################################## FUNCTIONS ##################################
def run(input_pdb, input_data, output_pdb, field, verbose=1, **kwargs):
    from os import devnull
    from os.path import expandvars
    from subprocess import Popen, PIPE
    import numpy as np
    import pandas as pd

    # Process arguments
    input_pdb  = expandvars(input_pdb)
    input_data = expandvars(input_data)
    output_pdb = expandvars(output_pdb)

    def iter_func():
        with open(input_pdb, "r") as open_file:
            for line in open_file:
                if line.startswith("TER") or line.startswith("ENDMDL"):
                    break
                if not line.startswith("ATOM"):
                    continue
                yield(line[6:11])           # Atom number
                yield(line[11:16])          # Atom name
                yield(line[16:20])          # Residue name
                yield(line[20:21])          # Conformation
                yield(line[21:22])          # Chain
                yield(line[22:26])          # Residue number
                yield(line[26:38])          # X
                yield(line[38:46])          # Y
                yield(line[46:54])          # Z
                yield(line[54:60])          # Beta
                yield(line[60:66])          # Occupancy
                yield(line[66:].strip())    # Element

    # Load input pdb
    if verbose >= 1:
        print("Loading input pdb file '{0}'".format(input_pdb))
    pdb_np = np.fromiter(iter_func(), dtype="S12").reshape((-1, 12))
    pdb = pd.DataFrame(
      data=np.array(pdb_np[:,1:]),
      index=np.array(pdb_np[:,0], dtype=np.int),
      columns=["atom_name", "residue_name", "conformation", "chain",
      "residue_number", "x", "y", "z", "beta", "occupancy", "element"])
    pdb.index.name = "atom_number"
    if verbose >= 2:
        print(pdb)

    # Load input data
    if verbose >= 1:
        print("Loading input data file '{0}'".format(input_data))
    data = pd.read_csv(input_data, index_col=0, delimiter=r"\s\s+")
    if verbose >= 2:
        print(data)

    # Store data in pdb field
    if verbose >= 1:
        print("Storing data in '{0}'".format(field))
    pdb[field] = pd.to_numeric(pdb[field])
    pdb[field] = 0
    for residue, row in data.iterrows():
        residue = int(residue.split(":")[1])
        pdb[field][pdb["residue_number"].astype(int) == residue] = \
        0.5 - row["delta cs"]/0.6
#        row["pre"]/2
    pdb[field][pdb[field] >= 0.5] = 0.5
    pdb[field][pdb[field] <= 0.0] = 0.0
    pdb[field][np.isnan(pdb[field])] = 0.75

    if verbose >= 1:
        print("Writing output pdb file '{0}'".format(output_pdb))
    with open(output_pdb, "w") as outfile:
        for index, row in pdb.iterrows():
            outfile.write("ATOM  " +
              "{0:5d}".format(index) +
              "{0:5s}".format(row["atom_name"]) + 
              "{0:4s}".format(row["residue_name"]) + 
              "{0:1s}".format(row["conformation"]) + 
              "{0:1s}".format(row["chain"]) + 
              "{0:1s}".format(row["residue_number"]) +
              "{0:12.3f}".format(float(row["x"])) +
              "{0:8.3f}".format(float(row["y"])) +
              "{0:8.3f}".format(float(row["z"])) +
              "{0:6.3f}".format(float(row["beta"])) +
              "{0:6.3f}".format(float(row["occupancy"])) +
              "{0:s}\n".format(row["element"]))

#################################### MAIN #####################################
if __name__ == "__main__":
    import argparse

    # Prepare argument parser
    parser            = argparse.ArgumentParser(
      description     = __doc__,
      formatter_class = argparse.RawTextHelpFormatter)
    parser.add_argument(
      "-b", "--beta",
      action   = "store_const",
      const    = "beta",
      default  = "beta",
      dest     = "field",
      help     = "Store value in beta field")
    parser.add_argument(
      "-o", "--occupancy",
      action   = "store_const",
      const    = "occupancy",
      default  = "beta",
      dest     = "field",
      help     = "Store value in occupancy field")
    parser.add_argument(
      "input_pdb",
      type     = str,
      help     = "input pdb file")
    parser.add_argument(
      "input_data",
      type     = str,
      help     = "input data to map onto pdb")
    parser.add_argument(
      "output_pdb",
      type     = str,
      help     = "output pdb file")

    verbosity = parser.add_mutually_exclusive_group()
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

    run(**kwargs)
