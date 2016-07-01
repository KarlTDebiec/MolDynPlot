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
################################### CLASSES ###################################
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

    if z is None or vmin is None or vmax is None:
        return LinearSegmentedColormap("test", cdict)
    else:
        return LinearSegmentedColormap("test", cdict)(z / (vmax - vmin))

def three_one(three):
    return {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "ASH": "D",
            "CYS": "C", "CYM": "C", "CYX": "C", "GLN": "Q", "GLU": "E",
            "GLH": "E", "GLY": "G", "HIS": "H", "HID": "H", "HIE": "H",
            "HIP": "H", "HSD": "H", "HSE": "H", "HSP": "H", "ILE": "I",
            "LEU": "L", "LYS": "K", "LYN": "K", "MET": "M", "PHE": "F",
            "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y",
            "VAL": "V" }[three.upper()]
