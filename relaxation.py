#!/usr/bin/python
# -*- coding: utf-8 -*-
#   moldynplot.relaxation.py
#
#   Copyright (C) 2012-2016 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Processes NMR relaxation and related data

.. todo:
  - document
  - Support rotational diffusion from simulation, including error from
    either block averaging of consecutive segments or standard error of
    independent simulations
  - Support rotational diffusion from experiment, including error from
    jackknife
"""
################################### MODULES ###################################
from __future__ import absolute_import,division,print_function,unicode_literals
if __name__ == "__main__":
    __package__ = str("moldynplot")
    import moldynplot
import h5py
import pandas as pd
import numpy as np
################################## FUNCTIONS ##################################
def spawn(function):
    """
    Wraps functions for use with :func:`multiprocess_map`

    Arguments:
      function (function): Function to run

    Returns:
      wrapped_function (function): wrapped function that accepts
        arguments from *queue_in* and outputs results of *function* to
        *queue_out*
    """
    def wrapped_function(queue_in, queue_out):
        while True:
            i, argument = queue_in.get()
            if i is None: break         # 'None' signals that queue is empty
            queue_out.put((i, function(argument)))

    return wrapped_function

def multiprocess_map(function, arguments, n_processes=1):
    """
    Runs *function* with *arguments* using *n_processes*.

    Replacement for multiproccessing.Pool.imap_unordered, which can only
    accept module-level functions.

    Arguments:
      function (function): Function to run
      arguments (iterable): Iterable of arguments to pass to function
      n_processes (int): Number of processes to use

    Returns:
      results (list): Results returned from *function*
    """
    from multiprocessing import Queue, Process

    # Initialize queues
    queue_in   = Queue(1)
    queue_out  = Queue()

    # Initialize processes and link to input and output queues
    processes = [Process(target = spawn(function),
      args = (queue_in, queue_out)) for i in range(n_processes)]
    for p in processes:
        p.daemon = True
        p.start()

    # Construct input queue, including 'None' signals to terminate
    input = [queue_in.put((i, argument))
      for i, argument in enumerate(arguments)]
    for i in range(n_processes):
        queue_in.put((None, None))

    # Retrieve output queue
    output = [queue_out.get() for i in range(len(input))]

    # Rejoin processes and return results
    for p in processes: p.join()
    return [x for i, x in sorted(output)]

def process_ired(infiles, outfile, indexfile=None, mode="mean",
    **kwargs):
    """
    Processes relaxation data calculated from MD simulations using the
    iRED method.

    Arguments:
      infiles (list): Path(s) to input file(s), may contain environment
        variables and wildcards
      outfile (str): Path to output text or hdf5 file
      indexfile (str): Path to index file used to map vector numbers
        listed in *infiles* to amino acid residues. Should contain one
        amino acid per line in the form 'XAA:#'
      mode (str): Relationship between input files; may be 'mean' if
        each infile represents an independent simulation, or
        'timeseries' if each infile represents a consecutive
        (potentially overlapping) excerpt of a longer simulation. If
        *mode* is 'mean', processed results will be the average across
        the simulations, including standard errors calculated using the
        standard deviation, if *mode* is 'timeseries', processed results
        will be a timeseries
      verbose (int): Level of verbose output
      kwargs (dict): Additional keyword arguments
    """
    from glob import glob
    from os import devnull
    from os.path import expandvars
    import re
    from subprocess import Popen, PIPE

#    # Process arguments
#    verbose = kwargs.get("verbose", 1)
#    processed_infiles = []
#    for infile in infiles:
#        print(infile)
#        processed_infiles.extend(sorted(glob(expandvars(infile))))
#    if len(processed_infiles) == 0:
#        print("No infiles found matching '{0}', exiting".format(infiles))
#        return
#    infiles = processed_infiles
#    if verbose >= 1:
#        print("Loading iRED data from {0} infiles, ".format(len(infiles)) +
#              "starting with '{0}'".format(infiles[0]))
#    outfile = expandvars(outfile)
#
#
#    # Load residue index
#    if indexfile is not None:
#        if verbose >= 1:
#            print("Loading residue indexes from '{0}'".format( indexfile))
#        res_index = np.loadtxt(expandvars(indexfile), dtype=np.str).flatten()
#
#    # Load data
#    r1r2noe_datasets = []
#    s2_datasets      = []
#    re_t1t2noe = re.compile(
#      "^#Vec\s+[\w_]+\[T1\]\s+[\w_]+\[\T2\]\s+[\w_]+\[NOE\]$")
#    re_s2 = re.compile("^#Vec\s+[\w_]+\[S2\]$")
#    for i, infile in enumerate(infiles):
#
#        # Determine if infile contains relaxation or order parameters
#        with open(devnull, "w") as fnull:
#            fields = Popen("head -n 1 {0}".format(infile),
#              stdout=PIPE, stderr=fnull, shell=True).stdout.read().strip()
#
#        # Parse relaxation
#        if re.match(re_t1t2noe, fields):
#            if verbose >= 1:
#                print("Loading iRED relaxation data from '{0}'".format(
#                  infile))
#            raw_data = np.loadtxt(infile, dtype=np.float32)
#            read_csv_kw = kwargs.get("read_csv_kw", dict(delim_whitespace=True,
#              header=0, index_col=0, names=["r1","r2","noe"]))
#            raw_data = pd.read_csv(infile, **read_csv_kw)
#            raw_data["r1"] = 1 / raw_data["r1"]
#            raw_data["r2"] = 1 / raw_data["r2"]
#            r1r2noe_datasets.append(raw_data)
#
#        # Parse order parameters
#        elif re.match(re_s2, fields):
#            if verbose >= 1:
#                print("Loading iRED order parameter data from '{0}'".format(
#                  infile))
#            raw_data = np.loadtxt(infile, dtype=np.float32)
#            read_csv_kw = kwargs.get("read_csv_kw", dict(delim_whitespace=True,
#              header=0, index_col=0, names=["s2"]))
#            raw_data = pd.read_csv(infile, **read_csv_kw)
#            s2_datasets.append(raw_data)
#
#        # Input file not understood
#        else:
#            raise Exception()
#
##    # Process index
##    items = []
##    fmt = []
##    if indexfile is not None:
##        items.append(("residue", residue))
##        fmt.append("%12s")
##    else:
##        fmt.append("%12d")
#
#    # Process relaxation
#    if mode == "mean":
#        if len(r1r2noe_datasets) == 1:
#            if verbose >= 1:
#                print("One one relaxation infile provided; skipping "
#                      "error calculation")
##            items.extend(
##                [("r1",  r1r2noe_datasets[0]["r1"]),
##                 ("r2",  r1r2noe_datasets[0]["r2"]),
##                 ("noe", r1r2noe_datasets[0]["noe"])])
##            fmt.extend(["%11.4f", "%11.4f", "%11.4f"])
#        elif len(r1r2noe_datasets) >= 2:
#            if verbose >= 1:
#                print("Calculating mean and standard error of "
#                      "{0} relaxation infiles".format(len(r1r2noe_datasets)))
##            r1r2noe_mean = pd.concat(r1r2noe_datasets).groupby(level=0).mean()
##            r1r2noe_se   = pd.concat(r1r2noe_datasets).groupby(
##                             level=0).std() / np.sqrt(len(r1r2noe_datasets))
##            items.extend(
##                [("r1",  r1r2noe_mean["r1"]),  ("r1 se",  r1r2noe_se["r1"]),
##                 ("r2",  r1r2noe_mean["r2"]),  ("r2 se",  r1r2noe_se["r2"]),
##                 ("noe", r1r2noe_mean["noe"]), ("noe se", r1r2noe_se["noe"])])
##            fmt.extend(["%11.4f", "%11.4f", "%11.4f",
##                        "%11.4f", "%11.4f", "%11.4f"])
#    elif mode == "timeseries":
#        if verbose >= 1:
#            print("Constructing timeseries from "
#                  "{0} relaxation infiles".format(len(r1r2noe_datasets)))
#            r1r2noe_ts = pd.concat([ds.stack() for ds in r1r2noe_datasets],
#                           axis=1).transpose()
#
#    # Process order parameters
#    if mode == "mean":
#        if len(s2_datasets) == 1:
#            if verbose >= 1:
#                print("One one order parameter infile provided; skipping "
#                      "error calculation")
##            items.extend([("s2",  s2_datasets[0]["s2"])])
##            fmt.extend(["%11.4f"])
#        elif len(s2_datasets) >= 2:
#            if verbose >= 1:
#                print("Calculating mean and standard error of "
#                      "{0} order parameter infiles".format(len(s2_datasets)))
##            s2_mean = pd.concat(s2_datasets).groupby(level=0).mean()
##            s2_se = pd.concat(s2_datasets).groupby(
##                      level=0).std() / np.sqrt(len(s2_datasets))
##            items.extend([("s2", s2_mean["s2"]), ("s2 se", s2_se["s2"])])
##            fmt.extend(["%11.4f", "%11.4f"])
#    elif mode == "timeseries":
#        if verbose >= 1:
#            print("Constructing timeseries from "
#                  "{0} order parameter infiles".format(len(r1r2noe_datasets)))
#            s2_ts = pd.concat([ds.stack() for ds in s2_datasets],
#                      axis=1).transpose()
#
#    # Organize data
#
#
#    # Organize data
#    if mode == "mean":
#        pass
#    elif mode == "timeseries":
#        timeseries = pd.merge(r1r2noe_ts, s2_ts, how="outer", left_index=True,
#                       right_index=True)
#
#        timeseries = timeseries[sorted(timeseries.columns.tolist(),
#                       key=lambda x: x[0])]
#        if res_index is not None:
#            timeseries.columns = timeseries.columns.set_levels(res_index,
#                                   level=0)
#        if verbose >= 1:
#            print("Writing timeseries to '{0}'".format(outfile))
#        if outfile.endswith(".h5") or outfile.endswith(".hdf5"):
#            with h5py.File(outfile) as hdf5_file:
#                hdf5_file.create_dataset("relax", data=timeseries,
#                  dtype=np.float32, chunks=True, compression="gzip",
#                  scaleoffset=5)
#                hdf5_file["relax"].attrs["columns"] = str(
#                  timeseries.columns.tolist())
#        else:
#            with open(outfile, "w") as text_file:
#                text_file.write(timeseries.to_string(col_space=12,
#                  sparsify=False))
#
#    data = pd.DataFrame.from_items(items)
#    if indexfile is not None:
#        data.set_index("residue", inplace=True)
#    else:
#        data.index.name = "vector"
#    columns = [data.index.name] + list(data.columns.values)
#    header = "{0:<11s}".format(columns.pop(0))
#    for column in columns:
#        header += "{0:>12s}".format(column)
#
#    # Output data
#    np.savetxt(outfile, np.column_stack((data.index.values, data.values)),
#      fmt=fmt, header=header, comments='#')

def process_error(sim_infiles, exp_infiles, outfile, **kwargs):
    """
    """
    from os import devnull

    if len(sim_infiles) != len(exp_infiles):
        raise ValueError("Number of simulation input files must match number "+
                         "of experimental input files, as they are treated "+
                         "pairwise. {0} simulation ".format(len(sim_infiles))+
                         "input file(s) and {0} ".format(len(exp_infiles))+
                         "experiment input file(s) provided.")

    # Work through each pair of infiles
    errs = []
    final_index = None
    for sim_infile, exp_infile in zip(sim_infiles, exp_infiles):
        print("Comparing simulation infile '{0}' ".format(sim_infile) +
              "with experimental infile '{0}':".format(exp_infile))

        # Load infiles and select shared indexes and columns
        sim = pd.read_csv(sim_infile, delim_whitespace=True, index_col=0)
        exp = pd.read_csv(exp_infile, delim_whitespace=True, index_col=0)
        overlap = sim.index.intersection(exp.index)
        if final_index is None:
            final_index = exp.index
        final_index = final_index.union(overlap)
        sim = sim.loc[overlap]
        exp = exp.loc[overlap]
        err_cols = [c for c in sim.columns.values
                      if not c.endswith(" se")
                      and c in exp.columns.values]
        err_se_cols = [c + " se" for c in err_cols
                       if  c + " se" in sim.columns.values
                       and c + " se" in exp.columns.values]
        print("   Files share fields {0} and {1} for {2} residues".format(
            str(map(str, err_cols)).replace("'", ""),
            str(map(str, err_se_cols)).replace("'", ""),
            len(overlap)))

        # Calculate error of available fields
        err = pd.DataFrame(0, index=overlap, columns=
          [x for t in zip(err_cols, err_se_cols) for x in t])
        err[err_cols] = (np.abs(exp[err_cols] - sim[err_cols])
                         / np.abs(exp[err_cols]))

        # Calculate uncertainty of error of available fields
        if len(err_se_cols) != 0:
            err[err_se_cols] = 0
            err[err_se_cols] = np.sqrt((err[err_cols].values) ** 2 *
                       ((np.sqrt(exp[err_se_cols].values ** 2
                       + sim[err_se_cols].values ** 2)
                         / (exp[err_cols].values - sim[err_cols].values)) ** 2
                        + (exp[err_se_cols].values / exp[err_cols].values) ** 2))
        errs.append(err)

    # Determine final columns and indexes
    final_cols = []
    final_index = sorted(final_index, key=lambda x: int(x.split(":")[1]))
    for err in errs:
        for col in err.columns.values:
            if not col in final_cols:
                final_cols.append(col)

    # Sum the columns
    final = pd.DataFrame(0.0, index=final_index, columns=final_cols)
    counts = pd.DataFrame(0, index=final_index, columns=final_cols)
    for err in errs:
        for col in err.columns.values:
            if not col.endswith(" se"):
                final[col].loc[err.index] += err[col].loc[err.index]
            else:
                final[col].loc[err.index] += err[col].loc[err.index] ** 2
            counts[col].loc[err.index] += 1

    # Average the columns
    print("Averaging fields:")
    for col in final_cols:
        if not col.endswith(" se"):
            print("    Averaging field '{0}'".format(col))
            final[col] /= counts[col]
        else:
            print("    Progagating uncertainty for field '{0}'".format(col))
            final[col] = np.sqrt(final[col]) / counts[col]

    # Write outfile
    print("Writing outfile '{0}' with fields ".format(outfile) +
          "{0} for ".format(str(map(str, final_cols)).replace("'", "")) +
          "{0} residues".format(len(final_index)))
    header = "residue    "
    for col in final_cols:
        header += "{0:>12s}".format(col)
    fmt = ["%12s"] + ["%11.5f"] * len(final_cols)
    np.savetxt(outfile, np.column_stack((final.index.values,
      final.values)), fmt=fmt, header=header, comments='#')

def process_relax(relax_type, peaklist, infiles, delays, error_method,
    n_synth_datasets, outfile, verbose=1, **kwargs):
    """
    Processes experimental relaxation data.

    Arguments:
      relax_type (str):
      peaklist (str):
      delays (list):
      error_
      infiles (list): Path(s) to input file(s), may contain environment
        variables and wildcards
      outfile (str): Path to output text file

    """
    from glob import glob
    from os.path import expandvars
    import nmrglue
    from scipy.interpolate import RectBivariateSpline
    from scipy.optimize import curve_fit

    # Process arguments
    processed_infiles = []
    for infile in infiles:
        processed_infiles += glob(expandvars(infile))
    infiles = processed_infiles
    if len(delays) != len(infiles):
        raise()
    peaklist = expandvars(peaklist)
    outfile = expandvars(outfile)

    # Load peaklist
    if verbose >= 1:
        print("Loading peaklist from '{0}'".format(peaklist))
    def convert_name(name):
        return "{0}:{1}".format(name[-4:-1].upper(), name[2:-4])
    relax = pd.read_csv(peaklist, sep="\t", usecols=[2,3,4], index_col=2,
      converters={4:convert_name}, names=["1H", "15N", "residue"], skiprows=1)

    # Load peak intensities from spectra
    for infile, delay in zip(infiles, delays):
        if verbose >= 1:
            print("Loading intensities from '{0}'".format(infile))
        parameters, intensity = nmrglue.pipe.read(infile)
        hydrogen = nmrglue.pipe.make_uc(parameters, intensity,
                     dim=1).ppm_scale()
        nitrogen = nmrglue.pipe.make_uc(parameters, intensity,
                     dim=0).ppm_scale()

        def calc_intensity(peak, **kwargs):
            H_index = np.argmin((hydrogen - peak["1H"])  ** 2)
            N_index = np.argmin((nitrogen - peak["15N"]) ** 2)
            return intensity[N_index, H_index]

        relax["{0} ms".format(delay)] = relax.apply(calc_intensity, axis=1)

    # Calculate relaxation rates
    delays = np.array(delays, np.float64) / 1000
    def calc_relax(peak, **kwargs):
        if verbose >= 1:
            print("Calculating relaxation for {0}".format(peak.name))

        def model_function(delay, intensity, relaxation):
            return intensity * np.exp(-1 * delay * relaxation)

        I = np.array(peak.filter(regex=(".*ms")).values, np.float64)
        I0, R = curve_fit(model_function, delays, I, p0=(I[0], 1.0))[0]

        # Calculate error
        if error_method == "rmse":
            error = np.sqrt(np.mean((I - model_function(delays, I0, R)) ** 2))
        elif error_method == "mae":
            error = np.mean(np.sqrt((I - model_function(delays, I0, R)) ** 2))

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
                    print("Unable to calculate standard error for {0}".format(
                      peak.name))
                return np.nan

        # Calculate standard error
        synth_Rs = multiprocess_map(synth_fit_decay, synth_datasets, 16)
        R_se = np.std(synth_Rs)

        return pd.Series([I0, R, R_se])

    # Calculate relaxation rates and standard errors
    fit = relax.apply(calc_relax, axis=1)
    fit.columns = ["I0", relax_type, relax_type + " se"]
    relax = relax.join(fit)

    # Write outfile
    if verbose >= 1:
        print("Writing outfile '{0}'".format(outfile))
    columns = [relax.index.name] + list(relax.columns.values)
    header = "{0:<11s}".format(columns.pop(0))
    for column in columns:
        header += "{0:>12s}".format(column)
    fmt = ["%12s", "%11.4f", "%11.4f"] + ["%11d"] * len(delays) + \
         ["%11d", "%11.4f", "%11.4f"]
    np.savetxt(outfile, np.column_stack((relax.index.values,
      relax.values)), fmt=fmt, header=header, comments='#')

def process_hetnoe(peaklist, infiles, outfile, verbose=1, **kwargs):
    """
    """
    from glob import glob
    from os.path import expandvars
    import nmrglue

    # Process arguments
    processed_infiles = []
    for infile in infiles:
        processed_infiles += glob(expandvars(infile))
    infiles = processed_infiles
    if len(infiles) != 2:
        raise()
    peaklist = expandvars(peaklist)
    outfile = expandvars(outfile)

    # Load peaklist
    if verbose >= 1:
        print("Loading peaklist from '{0}'".format(peaklist))
    def convert_name(name):
        return "{0}:{1}".format(name[-4:-1].upper(), name[2:-4])
    relax = pd.read_csv(peaklist, sep="\t", usecols=[2,3,4], index_col=2,
      converters={4:convert_name}, names=["1H", "15N", "residue"], skiprows=1)

    # Load peak intensities from spectra
    def calc_intensity(peak, **kwargs):
        H_index = np.argmin((hydrogen - peak["1H"])  ** 2)
        N_index = np.argmin((nitrogen - peak["15N"]) ** 2)
        return intensity[N_index, H_index]

    if verbose >= 1:
        print("Loading intensities from '{0}'".format(infiles[0]))
    parameters, intensity = nmrglue.pipe.read(infiles[0])
    hydrogen = nmrglue.pipe.make_uc(parameters, intensity,
                 dim=1).ppm_scale()
    nitrogen = nmrglue.pipe.make_uc(parameters, intensity,
                 dim=0).ppm_scale()
    hydrogen += 0.0612858
    nitrogen += 0.08399

    relax["sat"] = relax.apply(calc_intensity, axis=1)
    sat_se = intensity[np.logical_and(intensity > -intensity.std(),
      intensity < intensity.std())].std()

    if verbose >= 1:
        print("Loading intensities from '{0}'".format(infiles[1]))
    parameters, intensity = nmrglue.pipe.read(infiles[1])
    relax["nosat"] = relax.apply(calc_intensity, axis=1)
    nosat_se = intensity[np.logical_and(intensity > -intensity.std(),
      intensity < intensity.std())].std()

    relax["noe"] = relax["sat"] / relax["nosat"]
    relax["noe_se"] = np.sqrt((sat_se / relax["sat"]) ** 2 +
      (nosat_se / relax["nosat"]) ** 2) * relax["noe"]

    # Write outfile
    if verbose >= 1:
        print("Writing outfile '{0}'".format(outfile))
    columns = [relax.index.name] + list(relax.columns.values)
    header = "{0:<11s}".format(columns.pop(0))
    for column in columns:
        header += "{0:>12s}".format(column)
    fmt = ["%12s", "%11.4f", "%11.4f"] + ["%11d"] * 2 + \
         ["%11.4f", "%11.4f"]
    np.savetxt(outfile, np.column_stack((relax.index.values,
      relax.values)), fmt=fmt, header=header, comments='#')

def process_pre(dia_infile, para_infile, outfile, verbose=1, **kwargs):
    """
    """
    from glob import glob
    from os.path import expandvars
    import nmrglue

    # Process arguments
    dia_infile  = glob(expandvars(dia_infile))[0]
    para_infile = glob(expandvars(para_infile))[0]

    if verbose >= 1:
        print("Loading diamagnetic relaxation rates from '{0}'".format(
          dia_infile))
    dia_relax = pd.read_csv(dia_infile, index_col=0, delimiter=r"\s\s+")
    dia_relax.index.name = "residue"
    if verbose >= 1:
        print("Loading paramagnetic relaxation rates from '{0}'".format(
          para_infile))
    para_relax = pd.read_csv(para_infile, index_col=0, delimiter=r"\s\s+")
    para_relax.index.name = "residue"

    relax = dia_relax[["1H", "15N", "dia", "dia_se"]]
    relax = pd.concat((relax, para_relax[["para", "para_se"]]), axis=1)

    relax["pre"] = relax["dia"] / relax["para"]
    relax["pre_se"] = np.sqrt((relax["dia_se"] / relax["dia"]) ** 2 +
      (relax["para_se"] / relax["para"]) ** 2) * relax["pre"]
    relax["pre"][np.isinf(relax["pre"])] = 0
    relax["pre_se"][np.isnan(relax["pre_se"])] = 0
    print(relax)

    # Write outfile
    if verbose >= 1:
        print("Writing outfile '{0}'".format(outfile))
    columns = [relax.index.name] + list(relax.columns.values)
    header = "{0:<11s}".format(columns.pop(0))
    for column in columns:
        header += "{0:>12s}".format(column)
    fmt = ["%12s"] + ["%11.4f"]*8
    np.savetxt(outfile, np.column_stack((relax.index.values,
      relax.values)), fmt=fmt, header=header, comments='#')

#################################### MAIN #####################################
if __name__ == "__main__":
    import argparse

    # Prepare argument parser
    epilog = """
      Computer systems can fail at any time—crashes, bugs, power outages,
      whatever—so file systems need to anticipate and recover from these
      scenarios. The old-old-old school method is to plod along and then have a
      special utility to check and repair the file system during boot (fsck,
      short for file system check). More modern systems labor to achieve an
      always-consistent format, or only narrow windows of inconsistency,
      obviating the need for the full, expensive fsck. ZFS, for example, builds
      up new state on disk and then atomically transitions from the previous
      state to the new one with a single atomic operation."""
    parser = argparse.ArgumentParser(
      description     = """Processes datasets""",
      epilog          = epilog)
    subparsers = parser.add_subparsers(
                   dest        = "mode",
                   description = "")

    from .Dataset import IREDDataset, IREDTimeSeriesDataset

    IREDDataset.add_subparser(subparsers)
    IREDTimeSeriesDataset.add_subparser(subparsers)

#    # Prepare iRED subparser
#    ired_subparser = subparsers.add_parser(
#      name     = "ired",
#      help     = "Process simulated iRED relaxation data")
#    ired_subparser.set_defaults(
#      function   = process_ired)
#    input_group  = ired_subparser.add_argument_group("input")
#    action_group = ired_subparser.add_argument_group("action")
#    output_group = ired_subparser.add_argument_group("output")
#    input_group.add_argument(
#      "-infiles",
#      required = True,
#      dest     = "infiles",
#      metavar  = "INFILE",
#      nargs    = "+",
#      type     = str,
#      help     = "cpptraj iRED output file(s) from which to load datasets; " +
#                 "may be plain text or compressed, and may contain "
#                 "environment variables and wildcards")
#    input_group.add_argument(
#      "-indexfile",
#      required = False,
#      type     = str,
#      help     = "text file from which to load residue names; should list "
#                 "amino acids in the form 'XAA:#' separated by whitespace; if "
#                 "omitted will be taken from rows of first infile; may "
#                 "contain environment variables")
#
#    mode = input_group.add_mutually_exclusive_group()
#    mode.add_argument(
#      "--mean",
#      action   = "store_const",
#      const    = "mean",
#      default  = "mean",
#      dest     = "mode",
#      help     = "treat infiles as independent simulations; processed results "
#                 "will be the average across the simulations, including "
#                 "standard errors calculated using the standard deviation")
#    mode.add_argument(
#      "--timeseries",
#      action   = "store_const",
#      const    = "timeseries",
#      default  = "mean",
#      dest     = "mode",
#      help     = "treat infiles as consecutive (potentially overlapping) "
#                 "excerpts of a longer simulation; processed results will be "
#                 "a timeseries")
#    output_group.add_argument(
#      "-outfile",
#      required = True,
#      type     = str,
#      help     = "text or hdf5 file to which processed results will be "
#                 "output; may contain environment variables")

#    # Prepare error subparser
#    error_subparser  = subparsers.add_parser(
#      name        = "error",
#      help        = "Calculate error of simulated relaxation relative to " +
#                    "experiment",
#      description = "Calculates error of simulated relaxation relative to " +
#                    "experiment. The intended use case is to break down " +
#                    "errors relative to experimental data collected at " +
#                    "multiple magnetic fields or by multiple groups, " +
#                    "error(residue, measurement, magnet/group), into a form " +
#                    "that is easier to visualize and communicate, " +
#                    "error(residue, measurement). " +
#                    "Reads in a series of input files containing simulated " +
#                    "data and a series of files containing corresponding " +
#                    "experimental data. These files are treated in pairs " +
#                    "and the error between all data points present in both " +
#                    "(e.g. row 'GLN:2', column 'r1') calculated. " +
#                    "Columns ending in '_se' are treated as uncertainties, " +
#                    "and are propogated into uncertainties in the resulting " +
#                    "errors rather than being averaged. Take caution when " +
#                    "processing datasets that omit uncertainties alongside " +
#                    "those that do (experimental uncertainties are not " +
#                    "always reported), as the resulting uncertainties in " +
#                    "the residuals will be incorrect.")
#    error_subparser.set_defaults(
#      function   = process_error)
#    input_group  = error_subparser.add_argument_group("input")
#    action_group = error_subparser.add_argument_group("action")
#    output_group = error_subparser.add_argument_group("output")
#    input_group.add_argument(
#      "-sim_infile",
#      required = True,
#      dest     = "sim_infiles",
#      nargs    = "+",
#      type     = str,
#      help     = "input file(s) from which to load simulation datasets")
#    input_group.add_argument(
#      "-exp_infile",
#      required = True,
#      dest     = "exp_infiles",
#      nargs    = "+",
#      type     = str,
#      help     = "input file(s) from which to load experimental datasets")
#    output_group.add_argument(
#      "-outfile",
#      required = True,
#      type     = str,
#      help     = "Text file to which processed data will be output")
#
#    # Prepare relax subparser
#    relax_subparser  = subparsers.add_parser(
#      name     = "relax",
#      help     = "Process experimental R1 or R2 relaxation data")
#    relax_subparser.set_defaults(
#      function   = process_relax)
#    input_group  = relax_subparser.add_argument_group("input")
#    action_group = relax_subparser.add_argument_group("action")
#    output_group = relax_subparser.add_argument_group("output")
#    relax_type = input_group.add_mutually_exclusive_group()
#    relax_type.add_argument(
#      "--r1",
#      action   = "store_const",
#      const    = "r1",
#      default  = "r1",
#      dest     = "relax_type",
#      help     = "process R1 relaxation data")
#    relax_type.add_argument(
#      "--r2",
#      action   = "store_const",
#      const    = "r2",
#      default  = "r1",
#      dest     = "relax_type",
#      help     = "process R2 relaxation data")
#    relax_type.add_argument(
#      "--pre-dia",
#      action   = "store_const",
#      const    = "dia",
#      default  = "r1",
#      dest     = "relax_type",
#      help     = "process PRE diamagnetic relaxation data")
#    relax_type.add_argument(
#      "--pre-para",
#      action   = "store_const",
#      const    = "para",
#      default  = "r1",
#      dest     = "relax_type",
#      help     = "process PRE paramagnetic relaxation data")
#    input_group.add_argument(
#      "-peaklist",
#      required = True,
#      type     = str,
#      help     = "peak list (exported from ccpnmr)")
#    input_group.add_argument(
#      "-infile",
#      required = True,
#      dest     = "infiles",
#      metavar  = "INFILE",
#      nargs    = "+",
#      type     = str,
#      help     = "NMR spectra (NMRPipe format)")
#    input_group.add_argument(
#      "-delay",
#      required = True,
#      dest     = "delays",
#      metavar  = "DELAY",
#      nargs    = "+",
#      type     = str,
#      help     = "delays (ms); number of delays must match number of infiles")
#    action_group.add_argument(
#      "-synthetics",
#      required = False,
#      dest     = "n_synth_datasets",
#      default  = 100,
#      type     = int,
#      help     = "number of synthetic datasets to use to calculate error")
#    error_method = action_group.add_mutually_exclusive_group()
#    error_method.add_argument(
#      "--rmse",
#      action   = "store_const",
#      const    = "rmse",
#      default  = "rmse",
#      dest     = "error_method",
#      help     = "use root mean square error to generate synthetic datasets")
#    error_method.add_argument(
#      "--mae",
#      action   = "store_const",
#      const    = "mae",
#      default  = "rmse",
#      dest     = "error_method",
#      help     = "use mean absolute error to generate synthetic datasets")
#    output_group.add_argument(
#      "-outfile",
#      required = True,
#      type     = str,
#      help     = "text file to which processed data will be output")
#
#    # Prepare hetnoe subparser
#    hetnoe_subparser  = subparsers.add_parser(
#      name     = "hetnoe",
#      help     = "Process experimental heteronuclear NOE relaxation data")
#    hetnoe_subparser.set_defaults(
#      function   = process_hetnoe)
#    input_group  = hetnoe_subparser.add_argument_group("input")
#    action_group = hetnoe_subparser.add_argument_group("action")
#    output_group = hetnoe_subparser.add_argument_group("output")
#    input_group.add_argument(
#      "-peaklist",
#      required = True,
#      type     = str,
#      help     = "peak list (exported from ccpnmr)")
#    input_group.add_argument(
#      "-infile",
#      required = True,
#      dest     = "infiles",
#      metavar  = "INFILE",
#      nargs    = 2,
#      type     = str,
#      help     = "NMR spectra (NMRPipe format)")
#    output_group.add_argument(
#      "-outfile",
#      required = True,
#      type     = str,
#      help     = "text file to which processed data will be output")
#
#    # Prepare pre subparser
#    pre_subparser  = subparsers.add_parser(
#      name     = "pre",
#      help     = "Process experimental heteronuclear NOE relaxation data")
#    pre_subparser.set_defaults(
#      function   = process_pre)
#    input_group  = pre_subparser.add_argument_group("input")
#    action_group = pre_subparser.add_argument_group("action")
#    output_group = pre_subparser.add_argument_group("output")
#    input_group.add_argument(
#      "-dia",
#      required = True,
#      dest     = "dia_infile",
#      metavar  = "DIA_INFILE",
#      type     = str,
#      help     = "Diamagnetic relaxation rates")
#    input_group.add_argument(
#      "-para",
#      required = True,
#      dest     = "para_infile",
#      metavar  = "PARA_INFILE",
#      type     = str,
#      help     = "Paramagnetic relaxation rates")
#    output_group.add_argument(
#      "-outfile",
#      required = True,
#      type     = str,
#      help     = "text file to which processed data will be output")
#
#    # Prepare pre subparser
#    format_subparser = subparsers.add_parser(
#      name     = "format",
#      help     = "Formats r1, r2, and heteronuclear NOE data")
#    format_subparser.set_defaults(
#      function   = process_format)
#    input_group  = format_subparser.add_argument_group("input")
#    action_group = format_subparser.add_argument_group("action")
#    output_group = format_subparser.add_argument_group("output")
#    input_group.add_argument(
#      "-r1",
#      required = True,
#      dest     = "r1_infile",
#      metavar  = "R1_INFILE",
#      type     = str,
#      help     = "R1 relaxation infile")
#    input_group.add_argument(
#      "-r2",
#      required = True,
#      dest     = "r2_infile",
#      metavar  = "R2_INFILE",
#      type     = str,
#      help     = "R2 relaxation infile")
#    input_group.add_argument(
#      "-hetnoe",
#      required = True,
#      dest     = "hetnoe_infile",
#      metavar  = "hetnoe_INFILE",
#      type     = str,
#      help     = "Heteronuclear NOE relaxation infile")
#    output_group.add_argument(
#      "-outfile",
#      required = True,
#      type     = str,
#      help     = "text file to which processed data will be output")
#
#    # Verbosity
#    for p in subparsers.choices.values():
#        verbosity = p.add_mutually_exclusive_group()
#        verbosity.add_argument(
#          "-v", "--verbose",
#          action   = "count",
#          default  = 1,
#          help     = "enable verbose output, may be specified more than once")
#        verbosity.add_argument(
#          "-q", "--quiet",
#          action   = "store_const",
#          const    = 0,
#          default  = 1,
#          dest     = "verbose",
#          help     = "disable verbose output")

    # Parse arguments and run selected function
    kwargs  = vars(parser.parse_args())
    kwargs.pop("function")(**kwargs)
