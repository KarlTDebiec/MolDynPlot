#!/usr/bin/python
# -*- coding: utf-8 -*-
#   moldynplot.PDistFigureManager.py
#
#   Copyright (C) 2015-2016 Karl T Debiec
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
    __package__ = str("moldynplot")
    import moldynplot
from .myplotspec.FigureManager import FigureManager
################################### CLASSES ###################################
class PDistFigureManager(FigureManager):
    """
    Manages the generation of probability distribution figures.
    """

    from .myplotspec.manage_defaults_presets import manage_defaults_presets
    from .myplotspec.manage_kwargs import manage_kwargs

    defaults = """
        draw_figure:
          subplot_kw:
            autoscale_on: False
        draw_subplot:
          xlabel: Distance (Å)
          xticks: [2,3,4,5,6,7,8,9,10]
          ylabel: Potential of Mean Force (kcal/mol)
          yticks: [-4.5,-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5]
        draw_dataset:
          xkey: center
          ykey: pmf
    """
    available_presets = """
      count:
        class: content
        draw_subplot:
          ylabel: Count
          yticks: [0,100,200,300,400,500,600,700,800,900,1000]
        draw_dataset:
          ykey: count
      probability:
        class: content
        draw_subplot:
          ylabel: Probability
          yticks: [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
        draw_dataset:
          ykey: probability
      free_energy:
        class: content
        draw_subplot:
          ylabel: Free Energy $(k_B T)$
          yticks: [0,1,2,3,4,5,6,7,8,9,10]
        draw_dataset:
          ykey: free energy
      pmf:
        class: content
        draw_subplot:
          ylabel: Potential of Mean Force (kcal/mol)
          yticks: [-4.5,-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5]
        draw_dataset:
          ykey: pmf
          zero_line: True
      notebook:
        class: target
        inherits: notebook
        draw_figure:
          left:       0.50
          sub_width:  5.00
          right:      0.25
          bottom:     1.00
          sub_height: 2.00
          hspace:     0.10
          top:        0.25
          title_kw:
            top: -0.1
          shared_legend: True
          shared_legend_kw:
            left:       0.50
            sub_width:  5.00
            right:      0.25
            bottom:     0.00
            sub_height: 0.50
            legend_kw:
              loc: 9
              ncol: 4
        draw_subplot:
          y2label_kw:
            labelpad: 12
    """

    @manage_defaults_presets()
    @manage_kwargs()
    def draw_dataset(self, subplot, xkey="x", ykey="pmf", label="", xoffset=0,
        yoffset=0, handles=None, verbose=1, debug=0, **kwargs):
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
        from .myplotspec import get_color, multi_get_copy
        from .H5Dataset import H5Dataset

        if kwargs.get("infile") is None:
            return

        plot_kw = multi_get_copy("plot_kw", kwargs, {})
        if "color" in plot_kw:
            plot_kw["color"] = get_color(plot_kw.pop("color"))
        elif "color" in kwargs:
            plot_kw["color"] = get_color(kwargs.pop("color"))

        # Load dataset, x, and y
        dataset = H5Dataset(default_address="pdist", default_key="pdist",
                    **kwargs)
        xmin, xmax = subplot.get_xbound()
        x = dataset.datasets["pdist"][xkey] + xoffset
        y = dataset.datasets["pdist"][ykey] + yoffset
        poi = kwargs.get("poi", False)
        if poi:
            import numpy as np
            from scipy.optimize import curve_fit
            xc = x[x >= 3.0]
            yc = y[-1 * xc.size:]
            xc = xc[xc <= 3.8]
            yc = yc[:xc.size]
            dPMFdx = (yc[1:] - yc[:-1]) / (xc[1] - xc[0])
            poi_index  = np.where(dPMFdx == np.nanmax(dPMFdx))[0][0]
            print(xc[poi_index])
            subplot.axvline(xc[poi_index], color=plot_kw.get("color"))

        # y=0 line; should really be in draw_subplot; must figure out
        #   how to derive from method without rewriting completely
        #   or disturbing argument collection
        if (kwargs.get("zero_line", False)
        and not hasattr(subplot, "_mps_zero_line")):
            subplot._mps_zero_line = subplot.axhline(0, color="k")

        # Plot
        handle = subplot.plot(x, y, label=label, **plot_kw)[0]
        if handles is not None:
            handles[label] = handle

#################################### MAIN #####################################
if __name__ == "__main__":
    PDistFigureManager().main()
########################## AWAITING REIMPLEMENTATION ##########################
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
