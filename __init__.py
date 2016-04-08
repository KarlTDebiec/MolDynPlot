#!/usr/bin/python
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
def dssp_color_palette():
    """
    Generates colormap for dssp

    Returns:
      cmap (LinearSegmentedColormap): colormap
    """
    from matplotlib.colors import LinearSegmentedColormap

#    palette = [(1.000, 1.000, 1.000),   # White
#               (0.573, 0.367, 0.256),   # Brown
#               (0.826, 0.504, 0.178),   # Orange
#               (0.755, 0.241, 0.246),   # Red
#               (0.304, 0.485, 0.632),   # Blue
#               (0.555, 0.365, 0.575),   # Purple
#               (0.369, 0.615, 0.361),   # Green
#               (0.854, 0.850, 0.338)]   # Yellow
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

    return LinearSegmentedColormap("test", cdict)
