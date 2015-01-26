#!/usr/bin/python
# -*- coding: utf-8 -*-
#   myplotspec_sim.StateProbFigureManager.py
#
#   Copyright (C) 2015 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Class to manage the generation of binding probability figures
"""
################################### MODULES ###################################
from __future__ import absolute_import,division,print_function,unicode_literals
if __name__ == "__main__":
    __package__ = str("myplotspec_sim")
    import myplotspec_sim
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
            fig_width:     5.00
            left:          1.20
            sub_width:     3.25
            right:         0.55
            fig_height:    7.50
            top:           1.35
            sub_height:    3.25
            bottom:        2.90
            subplot_kw:
              autoscale_on: False
            shared_legend:
              left:       1.20
              sub_width:  3.25
              sub_height: 2.30
              bottom:     0.50
              legend_kw:
                frameon:      False
                labelspacing: 0.5
                loc:          9
              legend_lw:  5
          draw_subplot:
            ylabel:       $P_{bound}$
            xticklabels:  []
            yticks:       [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
            yticklabels:  [0.0,"",0.2,"",0.4,"",0.6,"",0.8,"",1.0]
            tick_params:
              bottom:     off
              top:        off
            lw:           2
          draw_dataset:
            plot_kw:
              lw:         2.0
              align:      center
              width:      0.6
              ecolor:     black
              capsize:    4
              error_kw:
                elinewidth: 2.0
                capthick:   2.0
    """
    presets  = """
        presentation:
          draw_figure:
            fig_width:     5.00
            left:          1.20
            sub_width:     3.25
            right:         0.55
            fig_height:    7.50
            top:           1.35
            sub_height:    3.25
            bottom:        2.90
            title_fp:     24r
            label_fp:     24r
            shared_legend:
              left:       1.20
              sub_width:  3.24
              sub_height: 2.30
              bottom:     0.50
              legend_fp:  14r
              legend_kw:
                frameon:      False
                labelspacing: 0.5
                loc:          9
          draw_subplot:
            title_fp:     24r
            label_fp:     24r
            tick_fp:      18r
            legend_fp:    18r
            lw:           2
          draw_dataset:
            plot_kw:
                lw:       2
        presentation_three:
          draw_figure:
            ncols:      3
            fig_width:  10.00
            left:        1.00
            sub_width:   2.50
            wspace:      0.50
            right:       0.50
            fig_height:  7.50
            bottom:      3.00
            sub_height:  2.50
            top:         2.00
            subplots:
              1:
                ylabel: ""
              2:
                ylabel: ""
            shared_legend:
              left:       6.80
              sub_width:  2.50
              sub_height: 2.50
              bottom:     0.50
              legend_fp:  14r
              legend_kw:
                loc:          2
                labelspacing: 0.5
          draw_subplot:
            title_fp:     16r
            label_fp:     24r
            tick_fp:      16r
            legend_fp:    18r
            lw:           2
          draw_dataset:
            plot_kw:
              lw:         2
    """

    @manage_defaults_presets()
    @manage_kwargs()
    def draw_subplot(self, subplot, experiment = None, title = None,
        legend = None, shared_handles = None, **in_kwargs):
        """
        Draws a subplot

        **Arguments:**
            :*subplot*:        <Axes> on which to act
            :*experiment*:     Confidence interval of experimental
                               y value; region shaded in gray
            :*title*:          Subplot's title
            :*legend*:         Subplot's legend
            :*shared_handles*: Nascent OrderedDict of handles and
                               labels shared among subplots of host
                               figure
        """
        from collections import OrderedDict
        from .myplotspec.axes import set_xaxis, set_yaxis
        from .myplotspec.legend import set_legend
        from .myplotspec.text import set_title

        # Format
        set_xaxis(subplot, **in_kwargs)
        set_yaxis(subplot, **in_kwargs)
        if str(title) != "None":
            set_title(subplot, title = title, **in_kwargs)

        # Configure and plot datasets
        handles         = OrderedDict()
        dataset_specs   = in_kwargs.pop("datasets", {})
        dataset_indexes = sorted([i for i in dataset_specs.keys()
                            if str(i).isdigit()])

        # Plot experimental data
        if experiment is not None:
            subplot.axhspan(experiment[0], experiment[1], lw = 0,
              color = [0.7, 0.7, 0.7])
            handles["Experiment"] = subplot.plot([-10, -10], [-10, -10],
              color = [0.7, 0.7, 0.7], lw = 5)[0]

        for i in dataset_indexes:
            out_kwargs              = {}
            out_kwargs["debug"]     = in_kwargs.get("debug",    False)
            out_kwargs["preset"]    = in_kwargs.get("preset",    [])
            out_kwargs["yaml_dict"] = in_kwargs.get("yaml_dict", {})
            out_kwargs["yaml_keys"] = [key
              for key2 in [[key3 + ["datasets", "all"],
                            key3 + ["datasets", i]]
              for key3 in in_kwargs.get("yaml_keys")]
              for key  in key2]
            out_kwargs["x"]         = i + 0.75
            self.draw_dataset(subplot = subplot, handles = handles,
              **out_kwargs)

        # Draw legend
        if str(legend) != "None":
            set_legend(subplot, handles = handles, **in_kwargs)
        if str(shared_handles) != "None":
            for label, handle in handles.items():
                if label not in shared_handles:
                    shared_handles[label] = handle

    @manage_defaults_presets()
    @manage_kwargs()
    def draw_dataset(self, subplot, x = None, label = "", handles = None,
        debug = False, **kwargs):
        """
        Draws a dataset

        **Arguments:**
            :*subplot*: <Axes> on which to act
            :*handles*: Nascent OrderedDict of [labels]:handles on
                        subplot
            :*x*:       X coordinate
            :*label*:   Dataset label
            :*color*:   Dataset color
            :*plot_kw*: Keyword arguments to pass to subplot.plot()
        """
        from .myplotspec import get_color
        from .H5Dataset import H5Dataset

        plot_kw = kwargs.get("plot_kw", {})
        if "color" in plot_kw:
            plot_kw["color"] = get_color(plot_kw["color"])
        elif "color" in kwargs:
            plot_kw["color"] = get_color(kwargs.pop("color"))

        # Load dataset, x, and y
        if x is None:
            if handles is not None:
                x = 1 + len(handles)
            else:
                x = 1
        if not ("P unbound" in kwargs and "P unbound se" in kwargs):
            dataset = H5Dataset(
                        default_address = "assign/stateprobs",
                        default_key = "pbound", **kwargs)
            y    = 1.0 - dataset.datasets["pbound"]["P unbound"][0]
            yerr = dataset.datasets["pbound"]["P unbound se"][0] * 1.96
        else:
            y    = 1.0 - kwargs.pop("P unbound")
            yerr = kwargs.pop("P unbound se") * 1.96

        # Plot
        bar    = subplot.bar(x, y, yerr = yerr, **plot_kw)
        color  = bar.patches[0].get_facecolor()
        handle = subplot.plot([-10, -10], [-10, -10], color = color)[0]
        if handles is not None:
            handles[label] = handle

#################################### MAIN #####################################
if __name__ == "__main__":
    StateProbFigureManager().main()
