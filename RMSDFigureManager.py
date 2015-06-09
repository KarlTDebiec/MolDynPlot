#!/usr/bin/python
# -*- coding: utf-8 -*-
#   myplotspec_sim.RMSDFigureManager.py
#
#   Copyright (C) 2015 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Generates one or more RMSD figures to specifications in a YAML file.
"""
################################### MODULES ###################################
from __future__ import absolute_import,division,print_function,unicode_literals
if __name__ == "__main__":
    __package__ = str("myplotspec_sim")
    import myplotspec_sim
from .myplotspec.FigureManager import FigureManager
################################### CLASSES ###################################
class RMSDFigureManager(FigureManager):
    """
    Class to manage the generation of RMSD figures
    """

    from .myplotspec.manage_defaults_presets import manage_defaults_presets
    from .myplotspec.manage_kwargs import manage_kwargs
    from .myplotspec.manage_output import manage_output

    defaults = """
        draw_figure:
          subplot_kw:
            autoscale_on: False
        draw_subplot:
          xticks:       [0,200,400,600,800,1000]
          xlabel:       Time (ns)
          yticks:       [0,1,2,3,4,5]
          ylabel:       RMSD ($\\AA$)
          ylabel_kw:
            va:         center
            rotation:   vertical
        draw_dataset:
          xscale:       0.001
    """

    presets  = """
      notebook:
        inherits: notebook
        draw_figure:
          left:          0.50
          sub_width:     4.40
          right:         1.60
          bottom:        0.90
          sub_height:    1.80
          top:           0.40
          subplot2:
            left:        5.00
            sub_width:   1.20
            bottom:      0.90
            sub_height:  1.80
          shared_legend:
            left:        0.50
            sub_width:   5.70
            sub_height:  0.50
            bottom:      0.00
            legend_lw:  3
            legend_kw:
              frameon:      False
              labelspacing: 0.5
              legend_fp:    8r
              loc:          9
              ncol:         2
        draw_subplot:
          ylabel_kw:
            labelpad:   10
          subplot2_kw:
            xticks: [0.00,0.04,0.08,0.12,0.16]
            xlabel: Probability
            yticks: [0,1,2,3,4,5]
            yticklabels: []
            title_fp:     10b
            label_fp:     10b
            tick_fp:      8r
            legend_fp:    8r
            tick_params:
              length:     2
              pad:        6
        draw_dataset:
          plot_kw:
            lw:         0.5
          plot2_kw:
            lw:         1.0
    """

    @manage_defaults_presets()
    @manage_kwargs()
    @manage_output()
    def draw_figure(self, title=None, shared_xlabel=None,
        shared_ylabel=None, shared_legend=None, subplot2=None,
        **in_kwargs):
        from collections import OrderedDict
        from .myplotspec import get_figure_subplots
        from .myplotspec.legend import set_shared_legend
        from .myplotspec.text import set_title, set_shared_xlabel
        from .myplotspec.text import set_shared_ylabel

        # Prepare figure and subplots with specified dimensions
        subplot_specs = in_kwargs.pop("subplots", {})
        subplot_indexes = sorted([i for i in subplot_specs.keys()
                            if str(i).isdigit()])
        figure, subplots = get_figure_subplots(**in_kwargs)

        # Format Figure
        if title is not None:
            set_title(figure, title = title, **in_kwargs)
        if shared_xlabel is not None:
            set_shared_xlabel(figure, xlabel=shared_xlabel, **in_kwargs)
        if shared_ylabel is not None:
            set_shared_ylabel(figure, ylabel=shared_ylabel, **in_kwargs)
        if shared_legend is not None:
            shared_handles = OrderedDict()

        if subplot2 is not None:
            figure, subplot2s = get_figure_subplots(figure=figure,
              **subplot2)

        # Configure and plot subplots
        for i in subplot_indexes:
            out_kwargs = subplot_specs[i].copy()
            out_kwargs["verbose"] = in_kwargs.get("verbose", False)
            out_kwargs["debug"] = in_kwargs.get("debug", False)
            out_kwargs["preset"] = in_kwargs.get("preset", [])[:]
            out_kwargs["yaml_dict"] = in_kwargs.get("yaml_dict", {})
            out_kwargs["yaml_keys"] = [key
              for key2 in [[key3 + ["subplots", "all"],
                            key3 + ["subplots", i]]
              for key3 in in_kwargs.get("yaml_keys")]
              for key  in key2]
            if shared_legend is not None:
                out_kwargs["shared_handles"] = shared_handles
            if i in subplots:
                if i in subplot2s:
                    self.draw_subplot(subplot=subplots[i],
                      subplot2=subplot2s[i], **out_kwargs)
                else:
                    self.draw_subplot(subplot=subplots[i], **out_kwargs)

        # Draw legend
        if shared_legend is not None:
            set_shared_legend(figure, subplots, handles = shared_handles,
              **shared_legend)

        # Return results
        return figure

    @manage_defaults_presets()
    @manage_kwargs()
    def draw_subplot(self, subplot, subplot2=None, title=None,
        legend=None, shared_handles=None, **in_kwargs):
        from collections import OrderedDict
        from .myplotspec.axes import set_xaxis, set_yaxis
        from .myplotspec.legend import set_legend
        from .myplotspec.text import set_title

        # Format
        set_xaxis(subplot, **in_kwargs)
        set_yaxis(subplot, **in_kwargs)
        if title is not None:
            set_title(subplot, title = title, **in_kwargs)
        if subplot2 is not None:
            subplot2_kw = in_kwargs.get("subplot2_kw", {})
            set_xaxis(subplot2, **subplot2_kw)
            set_yaxis(subplot2, **subplot2_kw)
            subplot2.set_autoscale_on(False)

        # Configure and plot datasets
        handles = OrderedDict()
        dataset_specs = in_kwargs.pop("datasets", {})
        dataset_indexes = sorted([i for i in dataset_specs.keys()
                            if str(i).isdigit()])
        for i in dataset_indexes:
            out_kwargs = dataset_specs[i].copy()
            out_kwargs["verbose"] = in_kwargs.get("verbose", False)
            out_kwargs["debug"] = in_kwargs.get("debug", False)
            out_kwargs["preset"] = in_kwargs.get("preset", [])[:]
            out_kwargs["yaml_dict"] = in_kwargs.get("yaml_dict", {})
            out_kwargs["yaml_keys"] = [key
              for key2 in [[key3 + ["datasets", "all"],
                            key3 + ["datasets", i]]
              for key3 in in_kwargs.get("yaml_keys")]
              for key  in key2]
            if subplot2 is not None:
                out_kwargs["subplot2"] = subplot2
            self.draw_dataset(subplot=subplot, handles=handles,
              **out_kwargs)

        # Draw legend
        if legend is not None and legend is not False:
            set_legend(subplot, handles = handles, **in_kwargs)
        if shared_handles is not None:
            for label, handle in handles.items():
                if label not in shared_handles:
                    shared_handles[label] = handle

    @manage_defaults_presets()
    @manage_kwargs()
    def draw_dataset(self, subplot, subplot2=None, infile=None, label="",
        xscale=1, handles=None, verbose=False, debug=False, **kwargs):
        from .myplotspec import get_color
        import pandas
        from os.path import expandvars

        # Load data
        if  ("infile" is None):
            return
        if verbose >= 1:
            print("Loading {0}".format(expandvars(infile)))
        data = pandas.read_csv(expandvars(infile), delim_whitespace=True,
          comment="@", names=["time", "RMSD"])
        data["time"] *= xscale

        # Configure plot settings
        plot_kw = kwargs.get("plot_kw", {}).copy()
        if "color" in plot_kw:
            plot_kw["color"] = get_color(plot_kw["color"])
        elif "color" in kwargs:
            plot_kw["color"] = get_color(kwargs.pop("color"))

        # Plot
        handle = subplot.plot(data["time"], data["RMSD"], **plot_kw)[0]
        if handles is not None:
            handles[label] = handle
        if subplot2 is not None:
            import numpy as np
            from sklearn.neighbors import KernelDensity
            plot2_kw = plot_kw.copy()
            plot2_kw.update(kwargs.get("plot2_kw", {}))

            bandwidth   = 0.1
            grid        = np.linspace(0,5,100)
            kde = KernelDensity(bandwidth=bandwidth)
            kde.fit(data["RMSD"][:, np.newaxis])
            pdf = np.exp(kde.score_samples(grid[:, np.newaxis]))
            pdf /= pdf.sum()
            subplot2.plot(pdf, grid, **plot2_kw)

#################################### MAIN #####################################
if __name__ == "__main__":
    RMSDFigureManager().main()
