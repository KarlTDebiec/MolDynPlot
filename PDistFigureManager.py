#!/usr/bin/python
# -*- coding: utf-8 -*-
#   myplotspec_sim.PDistFigureManager.py
#
#   Copyright (C) 2015 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Generates one or more probability distribution figures to specifications
in a YAML file.
"""
################################### MODULES ###################################
from __future__ import absolute_import,division,print_function,unicode_literals
if __name__ == "__main__":
    __package__ = str("myplotspec_sim")
    import myplotspec_sim
from .myplotspec.FigureManager import FigureManager
################################### CLASSES ###################################
class PDistFigureManager(FigureManager):
    """
    Manages the generation of pdist figures using matplotlib.
    """

    from .myplotspec.manage_defaults_presets import manage_defaults_presets
    from .myplotspec.manage_kwargs import manage_kwargs

    defaults = """
        draw_figure:
          subplot_kw:
            autoscale_on: False
        draw_subplot:
          xticks: [0,1,2,3,4,5,6,7,8,9,10]
          xlabel: Distance (Å)
          yticks: [-4.5,-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5]
          ylabel: Potential of Mean Force (kcal/mol)
    """
    presets  = """
      count:
        draw_subplot:
          ylabel: Count
          yticks: [0,100,200,300,400,500,600,700,800,900,1000]
        draw_dataset:
          ykey: count
      probability:
        draw_subplot:
          ylabel: Probability
          yticks: [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
        draw_dataset:
          ykey: probability
      free_energy:
        draw_subplot:
          ylabel: Free Energy $(k_B T)$
          yticks: [0,1,2,3,4,5,6,7,8,9,10]
        draw_dataset:
          ykey: free energy
      pmf:
        draw_subplot:
          ylabel: Potential of Mean Force (kcal/mol)
          yticks: [-4.5,-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5]
        draw_dataset:
          ykey: pmf
          zero_line: True
      analogue:
        draw_subplot:
          xticks: [2,3,4,5,6,7,8]
          yticks: [-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5]
      notebook:
        help: Single plot for notebook (width ≤ 6.5", height ≤ 9")
        inherits: notebook
        draw_figure:
          left:       0.60
          sub_width:  4.70
          wspace:     0.30
          right:      0.20
          bottom:     0.50
          sub_height: 2.40
          top:        0.30
        draw_subplot:
          legend: False
          y2label_kw:
            labelpad: 6
    """

    @manage_defaults_presets()
    @manage_kwargs()
    def draw_dataset(self, subplot, xkey="center", ykey="pmf", label="",
        xoffset=0, yoffset=0, handles=None, debug=False, **kwargs):
        """
        Draws a dataset.

        Arguments:
          subplot (Axes): Axes on which to draw
          xkey (str): Name of field in dataset from which to load x
          ykey (str): Name of field in dataset from which to load y
          label (str, optional): Dataset label
          color (str, list, ndarray, float, optional): Dataset color
          plot_kw (dict, optional): Additional keyword arguments passed
            to subplot.plot()
          handles (OrderedDict, optional): Nascent OrderedDict of
            [labels]: handles on subplot
          zero_line (bool): Draw a horizontal line at y=0; really a
            subplot setting, but placed here to avoid needing to
            override draw_subplot
          kwargs (dict): Additional keyword arguments
        """
        from copy import copy
        from .myplotspec import get_color
        from .H5Dataset import H5Dataset

        plot_kw = copy(kwargs.get("plot_kw", {}))
        if "color" in plot_kw:
            plot_kw["color"] = get_color(plot_kw.pop("color"))
        elif "color" in kwargs:
            plot_kw["color"] = get_color(kwargs.pop("color"))

        # Load dataset, x, and y
        dataset = H5Dataset(default_address = "pdist", default_key = "pdist",
                    **kwargs)
        xmin, xmax = subplot.get_xbound()
        x = copy(dataset.datasets["pdist"][xkey])
        y = copy(dataset.datasets["pdist"][ykey])
        x = x[x >= xmin]
        y = y[-1 * x.size:]
        x = x[x <= xmax]
        y = y[:x.size]
        import numpy as np
        y -= np.min(y)
        dydx = (y[1:] - y[:-1]) / (x[1] - x[0])
        min_index  = (np.abs(x - 2.8)).argmin()
        max_index  = (np.abs(x - 3.2)).argmin()
        poi_index  = np.where(dydx == np.nanmax(
          dydx[min_index:max_index]))[0][0]
        min_index  = (np.abs(x - x[poi_index] + 0.10)).argmin()
        max_index  = (np.abs(x - x[poi_index] - 0.10)).argmin()
        def func(x, a, b, c):
            return a * (x - b) ** 2 + c
        from scipy.optimize import curve_fit
        a, b, c = curve_fit(func, x[min_index:max_index],
          dydx[min_index:max_index], p0 = [-100, 3, 4], maxfev = 10000)[0]
        subplot.axvline(b, linewidth = 1.0, color = plot_kw["color"])
        c = (np.abs(dataset.datasets["pdist"][xkey] - b)).argmin()
        print("POI = {0:4.2f}".format(b))
        # y=0 line; should really be in draw_subplot; must figure out
        #   how to derive from method without rewriting completely
        #   or disturbing argument collection
#        subplot.plot([xmin, xmax], [0, 0], color = "black")

        # Plot
        handle = subplot.plot(x, y, label = label, **plot_kw)[0]
        if handles is not None:
            handles[label] = handle

#################################### MAIN #####################################
if __name__ == "__main__":
    PDistFigureManager().main()
########################## AWAITING REIMPLEMENTATION ##########################
#def plot_three(content, **kwargs):
#    figure, subplots = get_figure_subplots(
#      nrows      = 2,    ncols  = 2,
#      fig_width  = 10.0, left   = 1.0, sub_width  = 4.00, wspace = 0.5,
#      fig_height =  7.5, bottom = 1.0, sub_height = 2.40, hspace = 0.5)
#def poi_subplot(subplot, datasets, **kwargs):
#        dPMFdx     = (y[1:] - y[:-1]) / (x[1] - x[0])
#        dPMFdx     = np.mean(
#          np.reshape(dPMFdx[:(dPMFdx.size - dPMFdx.size % downsample)],
#          (dPMFdx.size   / downsample, downsample)), axis = 1)
#        xp         = (x[:-1]  + x[1:])  / 2
#        min_index  = (np.abs(xp - 2.8)).argmin()
#        max_index  = (np.abs(xp - 3.2)).argmin()
#        poi_index  = np.where(dPMFdx == np.nanmax(
#          dPMFdx[min_index:max_index]))[0][0]
#        min_index  = (np.abs(xp - xp[poi_index] + 0.10)).argmin()
#        max_index  = (np.abs(xp - xp[poi_index] - 0.10)).argmin()
#        def func(x, a, b, c): return a * (x - b) ** 2 + c
#        a, b, c = curve_fit(func, xp[min_index:max_index],
#          dPMFdx[min_index:max_index], p0 = [-100, 3, 4], maxfev = 10000)[0]
#        subplot.plot(xp, dPMFdx, linewidth = 2.0, color = color)
#        subplot.axvline(b, linewidth = 1.0, color = color)
#        set_inset(subplot, "POI = {0:4.2f}".format(np.round(b,2)), fp = "14r",
#          xpro = 0.95, ypro = 0.95, ha = "right", va = "top")
#def plot_poi(datasets, **kwargs):
#    figure, subplots = get_figure_subplots(
#      nrows  = 2,    ncols      = len(datasets),
#      left   = 1.00, sub_width  = 2.65, hspace = 0.0, right = 0.50,
#      bottom = 0.75, sub_height = 2.00, wspace = 0.0, top   = 1.25)
