#!/usr/bin/python
# -*- coding: utf-8 -*-
#   moldynplot.StateProbFigureManager.py
#
#   Copyright (C) 2015-2017 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Generates one or more state probability figures to specifications in a
YAML file.
"""
################################### MODULES ###################################
from __future__ import (absolute_import, division, print_function,
    unicode_literals)

if __name__ == "__main__":
    __package__ = str("moldynplot")
    import moldynplot

from .myplotspec.FigureManager import FigureManager
from .myplotspec.manage_defaults_presets import manage_defaults_presets
from .myplotspec.manage_kwargs import manage_kwargs


################################### CLASSES ###################################
class StateProbFigureManager(FigureManager):
    """
    Class to manage the generation of probability distribution figures
    """

    defaults = """
        draw_figure:
          subplot_kw:
            autoscale_on: False
            axisbg: none
          multi_tick_params:
            bottom: off
            left: on
            right: off
            top: off
          shared_legend: True
          shared_legend_kw:
            spines: False
            handle_kw:
              ls: none
              marker: s
              mec: black
            legend_kw:
              frameon: False
              loc: 9
              numpoints: 1
              handletextpad: 0.0
              borderaxespad: 0
        draw_subplot:
          tick_params:
            bottom: off
            direction: out
            left: on
            right: off
            top: off
          title_kw:
            verticalalignment: bottom
          grid: True
          grid_kw:
            axis: y
            b: True
            color: [0.7,0.7,0.7]
            linestyle: '-'
          xticklabels:  []
          yticks: [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
        draw_dataset:
          bar_kw:
            align: center
            width: 0.6
            ecolor: black
            zorder: 3
            error_kw:
              zorder: 4
          handle_kw:
            ls: none
            marker: s
            mec: black
    """

    available_presets = """
      pbound:
        help: Probability of bound state
        draw_subplot:
          ylabel: $P_{bound}$
          yticklabels:  [0.0,"",0.2,"",0.4,"",0.6,"",0.8,"",1.0]
      presentation:
        class: target
        inherits: presentation
        draw_figure:
          left:       6.70
          sub_width:  3.20
          sub_height: 3.20
          bottom:     2.40
          shared_legend: True
          shared_legend_kw:
            left:        6.80
            sub_width:   3.00
            sub_height:  2.38
            bottom:      0.00
            legend_kw:
              columnspacing: 1
              labelspacing: 0.8
              legend_fp: 12r
              loc: 9
              ncol: 2
        draw_subplot:
          title_fp: 16b
          label_fp: 16b
          tick_fp:  12r
          tick_params:
            width: 2
        draw_dataset:
          bar_kw:
            lw: 2
            capsize: 4
            error_kw:
              capsize: 5
              capthick: 2
              elinewidth: 2
          handle_kw:
            ms: 10
            mew: 2
      manuscript:
        class: target
        inherits: manuscript
        draw_figure:
          left:       0.40
          sub_width:  1.50
          wspace:     0.10
          right:      0.10
          top:        0.35
          sub_height: 1.50
          bottom:     0.55
          shared_legend: True
          shared_legend_kw:
            left:        0.40
            sub_width:   1.50
            sub_height:  0.50
            bottom:      0.00
            legend_kw:
              columnspacing: 0.5
              labelspacing: 0.5
              legend_fp: 8r
              ncol: 4
        draw_dataset:
          bar_kw:
            error_kw:
              capsize: 2
              capthick: 0.5
              elinewidth: 0.5
            linewidth: 0.5
            edgecolor: [0,0,0]
          handle_kw:
            markeredgewidth: 0.5
            markersize: 5
    """

    @manage_defaults_presets()
    @manage_kwargs()
    def draw_dataset(self, subplot, experiment=None, x=None, label="",
            handles=None, bar=True, verbose=1, debug=0, **kwargs):
        """
        Draws a dataset.

        Arguments:
          subplot (Axes): Axes on which to draw
          x (float): X coordinate of bar
          label (str, optional): Dataset label
          color (str, list, ndarray, float, optional): Dataset color
          bar_kw (dict, optional): Additional keyword arguments passed
            to subplot.plot()
          handles (OrderedDict, optional): Nascent OrderedDict of
            [labels]: handles on subplot
          kwargs (dict): Additional keyword arguments
        """
        from .myplotspec import get_color, multi_get_copy
        from .dataset import H5Dataset

        # Handle missing input gracefully
        handle_kw = multi_get_copy("handle_kw", kwargs, {})
        if experiment is not None:
            subplot.axhspan(experiment[0], experiment[1], lw=2,
                color=[0.6, 0.6, 0.6])
            handles["Experiment"] = \
            subplot.plot([-10, -10], [-10, -10], mfc=[0.6, 0.6, 0.6],
                **handle_kw)[0]
            return
        if "infile" not in kwargs:
            if "P unbound" in kwargs and "P unbound se" in kwargs:
                y = 1.0 - kwargs.pop("P unbound")
                yerr = kwargs.pop("P unbound se") * 1.96
            elif "y" in kwargs and "y se" in kwargs:
                y = kwargs.pop("y")
                yerr = kwargs.pop("y se") * 1.96
            else:
                return
        else:
            dataset = H5Dataset(default_address="assign/stateprobs",
                default_key="pbound", **kwargs)
            y = 1.0 - dataset.datasets["pbound"]["P unbound"][0]
            yerr = dataset.datasets["pbound"]["P unbound se"][0] * 1.96

        # Configure plot settings

        # Plot
        if bar:
            bar_kw = multi_get_copy("bar_kw", kwargs, {})
            for color_key in ["color", "edgecolor", "facecolor", "ecolor"]:
                if color_key in bar_kw:
                    bar_kw[color_key] = get_color(bar_kw[color_key])
            barplot = subplot.bar(x, y, yerr=yerr, **bar_kw)

            handle_kw = multi_get_copy("handle_kw", kwargs, {})
            handle_kw["mfc"] = barplot.patches[0].get_facecolor()
            handle = subplot.plot([-10, -10], [-10, -10], **handle_kw)[0]
        if handles is not None and label is not None:
            handles[label] = handle


#################################### MAIN #####################################
if __name__ == "__main__":
    StateProbFigureManager().main()
