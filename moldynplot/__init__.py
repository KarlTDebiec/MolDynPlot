# -*- coding: utf-8 -*-
#   moldynplot.__init__.py
#
#   Copyright (C) 2015-2016 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
General functions.
"""
################################### MODULES ###################################
from __future__ import absolute_import,division,print_function,unicode_literals
################################## FUNCTIONS ##################################
def dssp_cmap(z=None, vmin=None, vmax=None):
    """
    Generates colormap for Dictionary of Secondary Structure Prediction
    (DSSP) figures.

    Arguments:
      z (float, optional): Position along color axis
      vmin (float, optional): Lower bound of color axis
      vmax (float, optional): Upper bound of color axis

    Returns:
      LinearSegmentedColormap: DSSP color map; if *z*, *vmin*, and
      *vmax* are provided, provides the color at position *z* along the
      axis between *vmin* and *vmax*
    """
    from matplotlib.colors import LinearSegmentedColormap

    palette = [(1.000, 1.000, 1.000),   # White
               (0.573, 0.367, 0.256),   # Brown
               (0.826, 0.504, 0.178),   # Orange
               (0.769, 0.306, 0.321),   # Red
               (0.298, 0.447, 0.690),   # Blue
               (0.506, 0.447, 0.698),   # Purple
               (0.333, 0.659, 0.408),   # Green
               (0.769, 0.678, 0.400)]   # Yellow
    cdict = {"red": [], "green": [], "blue": []}
    for i, (red, green, blue) in enumerate(palette):
        cdict["red"]   += [(i / (len(palette) - 1), red,   red)]
        cdict["green"] += [(i / (len(palette) - 1), green, green)]
        cdict["blue"]  += [(i / (len(palette) - 1), blue,  blue)]
    cdict["red"]   = tuple(cdict["red"])
    cdict["green"] = tuple(cdict["green"])
    cdict["blue"]  = tuple(cdict["blue"])

    cmap = LinearSegmentedColormap("dssp", cdict)
    if z is None or vmin is None or vmax is None:
        return cmap
    else:
        return cmap(z / (vmax - vmin))

def ff99SB_cmap(z=None, vmin=None, vmax=None):
    """
    Generates purple->yellow->black colormap.

    Generates colormap in style of:
      Hornak, Viktor, Abel, Robert, Okur, Asim, Stockbine, Bentley,
      Roitberg, Adrian, Simmerling, Carlos, Comparison of Multiple Amber
      Force Fields and Development of Improved Protein Backbone
      Parameters. Proteins: Structure, Function, and Bioinformatics.
      2006. 65. 712-725.

    Arguments:
      z (float, optional): Position along color axis
      vmin (float, optional): Lower bound of color axis
      vmax (float, optional): Upper bound of color axis

    Returns:
      LinearSegmentedColormap: ff99SB-style color map; if *z*, *vmin*,
      and *vmax* are provided, provides the color at position *z* along
      the axis between *vmin* and *vmax*
    """
    from matplotlib.colors import LinearSegmentedColormap

    cdict = {"red":   [(0, 1, 1)],
             "green": [(0, 1, 1)],
             "blue":  [(0, 0, 0)]}
    for i in range(1, 193, 1):
        red   = 1.0
        green = 1.0 - ((i - 1) / 193)
        blue  = 0.0
        cdict["red"]   += [(i / 384, red,   red)]
        cdict["green"] += [(i / 384, green, green)]
        cdict["blue"]  += [(i / 384, blue,  blue)]
    for i in range(193, 385, 1):
        red   = (1.0 - ((i - 193) / 192)) * ((384 - i) / 192) ** 0.3
        green =  0.0
        blue  = (0.0 + ((i - 193) / 192)) * ((384 - i) / 192) ** 0.2
        cdict["red"]   += [(i / 384, red,   red)]
        cdict["green"] += [(i / 384, green, green)]
        cdict["blue"]  += [(i / 384, blue,  blue)]
    cdict["red"]   = tuple(cdict["red"])
    cdict["green"] = tuple(cdict["green"])
    cdict["blue"]  = tuple(cdict["blue"])

    cmap = LinearSegmentedColormap("ff99SB", cdict)
    if z is None or vmin is None or vmax is None:
        return cmap
    else:
        return cmap(z / (vmax - vmin))

def three_one(three):
    """
    Converts three-letter amino acid codes to one-letter.

    Arguments:
      three (str): Three letter amino acid code; AMBER and CHARMM
        nomenclature for alternative protonation states is supported,
        but lost by conversion.

    Returns:
      str: Corresponding one-letter amino acid code
    """
    return {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "ASH": "D",
            "CYS": "C", "CYM": "C", "CYX": "C", "GLN": "Q", "GLU": "E",
            "GLH": "E", "GLY": "G", "HIS": "H", "HID": "H", "HIE": "H",
            "HIP": "H", "HSD": "H", "HSE": "H", "HSP": "H", "ILE": "I",
            "LEU": "L", "LYS": "K", "LYN": "K", "MET": "M", "PHE": "F",
            "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y",
            "VAL": "V" }[three.upper()]

def multiprocess_map(function, arguments, n_processes=1):
    """
    Runs a function with arguments using n_processes.

    Meant as a replacement for :func:`multiproccessing.Pool.imap_unordered`,
    which can only accept module-level functions.

    Arguments:
      function (function): Function to run
      arguments (list): Iterable of arguments to pass to function
      n_processes (int): Number of processes to use

    Returns:
      list: results returned from function
    """
    from multiprocessing import Queue, Process

    def spawn(function):
        """
        Arguments:
          function: Function to run

        Returns:
          function: New Function that accepts arguments from queue_in
          and outputs results of function to queue_out
        """
        def run_function(queue_in, queue_out):
            while True:
                i, argument = queue_in.get()
                if i is None: break         # 'None' signals that queue is empty
                queue_out.put((i, function(argument)))
        return run_function

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

