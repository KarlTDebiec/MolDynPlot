#!/usr/bin/python
# -*- coding: utf-8 -*-
#   moldynplot.StateProbFigureManager.py
#
#   Copyright (C) 2015 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Generates one or more state probability figures to specifications in a
YAML file.
"""
################################### MODULES ###################################
from __future__ import absolute_import,division,print_function,unicode_literals
if __name__ == "__main__":
    __package__ = str("moldynplot")
    import moldynplot
from .myplotspec.FigureManager import FigureManager
################################### CLASSES ###################################
class StateProbFigureManager(FigureManager):
    """
    Class to manage the generation of probability distribution figures
    """

    from .myplotspec.manage_defaults_presets import manage_defaults_presets
    from .myplotspec.manage_kwargs import manage_kwargs

    defaults = """
        draw_figure:
          subplot_kw:
            autoscale_on: False
        draw_subplot:
          xticklabels:  []
          yticks: [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
          tick_params:
            bottom: off
            top: off
          title_kw:
            verticalalignment: bottom
        draw_dataset:
          plot_kw:
            lw: 2
            align: center
            width: 0.6
            ecolor: black
            capsize: 4
            error_kw:
              elinewidth: 1
              capthick:   1
              capsize:    3
    """

    available_presets = """
      pbound:
        help: Probability of bound state
        draw_subplot:
          ylabel: $P_{bound}$
          yticklabels:  [0.0,"",0.2,"",0.4,"",0.6,"",0.8,"",1.0]
      notebook:
        class: target
        inherits: notebook
        draw_figure:
          left:       0.60
          sub_width:  1.80
          wspace:     0.10
          right:      0.20
          top:        0.50
          sub_height: 1.80
          bottom:     0.65
          subplots:
            1:
              ylabel: ""
              yticklabels: []
            2:
              ylabel: ""
              yticklabels: []
          shared_legend: True
          shared_legend_kw:
            left:        0.60
            sub_width:   1.80
            sub_height:  0.60
            bottom:      0.00
            legend_kw:
              frameon: False
              labelspacing: 0.5
              legend_fp: 8r
              loc: 9
              ncol: 4
    """

    @manage_defaults_presets()
    @manage_kwargs()
    def draw_dataset(self, subplot, experiment=None, x=None, label="",
        handles=None, verbose=1, debug=0, **kwargs):
        """
        Draws a dataset.

        Arguments:
          subplot (Axes): Axes on which to draw
          x (float): X coordinate of bar
          label (str, optional): Dataset label
          color (str, list, ndarray, float, optional): Dataset color
          plot_kw (dict, optional): Additional keyword arguments passed
            to subplot.plot()
          handles (OrderedDict, optional): Nascent OrderedDict of
            [labels]: handles on subplot
          kwargs (dict): Additional keyword arguments
        """
        from .myplotspec import get_color, multi_get_copy
        from .H5Dataset import H5Dataset

        # Handle missing input gracefully
        print(label, kwargs.keys())
        if experiment is not None:
            subplot.axhspan(experiment[0], experiment[1], lw=0,
              color=[0.7, 0.7, 0.7])
            handles["Experiment"] = subplot.plot([-10, -10], [-10, -10],
              color=[0.7, 0.7, 0.7], lw=5)[0]
            return
        if "infile" not in kwargs:
            if "P unbound" not in kwargs or "P unbound se" not in kwargs:
                return
            else:
                y    = 1.0 - kwargs.pop("P unbound")
                yerr = kwargs.pop("P unbound se") * 1.96
        else:
            dataset = H5Dataset(
                        default_address="assign/stateprobs",
                        default_key="pbound", **kwargs)
            y    = 1.0 - dataset.datasets["pbound"]["P unbound"][0]
            yerr = dataset.datasets["pbound"]["P unbound se"][0] * 1.96

        # Configure plot settings
        plot_kw = multi_get_copy("plot_kw", kwargs, {})
        if "color" in plot_kw:
            plot_kw["color"] = get_color(plot_kw["color"])
        elif "color" in kwargs:
            plot_kw["color"] = get_color(kwargs.pop("color"))
        if "ecolor" in plot_kw:
            plot_kw["ecolor"] = get_color(plot_kw["ecolor"])
        elif "color" in kwargs:
            plot_kw["ecolor"] = get_color(kwargs.pop("ecolor"))
        elif "color" in plot_kw:
            plot_kw["ecolor"] = plot_kw["color"]

        # Load dataset, x, and y
        if x is None:
            if handles is not None:
                x = 1 + len(handles)
            else:
                x = 1
        print(x, y, yerr)

        # Plot
        bar    = subplot.bar(x, y, yerr=yerr, **plot_kw)
        color  = bar.patches[0].get_facecolor()
        handle = subplot.plot([-10, -10], [-10, -10], color=color)[0]
        if handles is not None:
            handles[label] = handle

#################################### MAIN #####################################
if __name__ == "__main__":
    StateProbFigureManager().main()
