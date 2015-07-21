#!/usr/bin/python
# -*- coding: utf-8 -*-
#   myplotspec_sim.TimeSeriesFigureManager.py
#
#   Copyright (C) 2015 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Generates one or more time series figures to specifications in a YAML
file.
"""
################################### MODULES ###################################
from __future__ import absolute_import,division,print_function,unicode_literals
if __name__ == "__main__":
    __package__ = str("myplotspec_sim")
    import myplotspec_sim
from .myplotspec.FigureManager import FigureManager
################################### CLASSES ###################################
class TimeSeriesFigureManager(FigureManager):
    """
    Class to manage the generation of time series figures
    """

    from .myplotspec.manage_defaults_presets import manage_defaults_presets
    from .myplotspec.manage_kwargs import manage_kwargs
    from .myplotspec.manage_output import manage_output

    defaults = """
        draw_figure:
          subplot_kw:
            autoscale_on: False
        draw_subplot:
          xticks: [0,200,400,600,800,1000]
          xlabel: Time (ns)
          yticks: [0,1,2,3,4,5,6,7,8,9,10]
          ylabel_kw:
            va: center
            rotation: vertical
    """

    presets  = """
      rmsd:
        help: Root Mean Standard Deviation (RMSD) vs. time
        draw_subplot:
          yticks: [0,1,2,3,4,5]
          ylabel: RMSD ($\\AA$)
          subplot2_kw:
            xticks: [0.00,0.05,0.10,0.15,0.20]
            xticklabels: ["0.00","0.05","0.10","0.15","0.20"]
            yticks: [0,1,2,3,4,5]
            y2ticks: [0,1,2,3,4,5]
        draw_dataset:
          load_kw:
            comment: "@"
            delim_whitespace: True
            names: ["time", "RMSD"]
          xscale: 0.001
          xkey: time
          ykey: RMSD
          grid: !!python/object/apply:numpy.linspace [0,5,100]
      radgyr:
        help: Radius of gyration vs. time
        draw_subplot:
          yticks: [10,11,12,13,14,15]
          ylabel: Radius of Gyration ($\\AA$)
          subplot2_kw:
            xticks: [0.00,0.05,0.10,0.15,0.20]
            xticklabels: ["0.00","0.05","0.10","0.15","0.20"]
            yticks: [10,11,12,13,14,15]
            y2ticks: [10,11,12,13,14,15]
        draw_dataset:
          load_kw:
            sequential_ds: True
            comment: "@"
            names: ["radgyr", "maxdist"]
          xscale: 0.1
          xkey: time
          ykey: radgyr
          grid: !!python/object/apply:numpy.linspace [10,15,100]
      notebook_pdist:
        help: Second plot including probability density
        extends: notebook
        draw_figure:
          right:         1.60
          subplot2:
            left:        5.00
            sub_width:   1.20
            bottom:      0.90
            sub_height:  1.80
          shared_legend:
            sub_width:   5.70
        draw_subplot:
          subplot2_kw:
            xticks: [0.0,0.2,0.4,0.6,0.8,1.0]
            xticklabels: ["0.0","0.2","0.4","0.6","0.8","1.0"]
            xlabel: Probability
            yticks: [0,1,2,3,4,5,6,7,8,9,10]
            yticklabels: []
            title_fp:  10b
            label_fp:  10b
            tick_fp:   8r
            legend_fp: 8r
            tick_params:
              length: 2
              pad: 6
        draw_dataset:
          plot2_kw:
            lw: 1.0
          bandwidth: 0.1
          grid: !!python/object/apply:numpy.linspace [0,10,1000]
      notebook:
        help: Single plot for notebook (width ≤ 6.5", height ≤ 9")
        inherits: notebook
        draw_figure:
          left:          0.50
          sub_width:     4.40
          right:         0.20
          bottom:        0.90
          sub_height:    1.80
          top:           0.40
          shared_legend:
            left:        0.50
            sub_width:   4.30
            sub_height:  0.50
            bottom:      0.00
            legend_lw: 3
            legend_kw:
              frameon: False
              labelspacing: 0.5
              legend_fp: 8r
              loc: 9
              ncol: 2
        draw_subplot:
          ylabel_kw:
            labelpad: 10
        draw_dataset:
          plot_kw:
            lw: 0.5
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
                if subplot2 is not None and i in subplot2s:
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
    def draw_dataset(self, subplot, xkey, ykey, subplot2=None,
        infile=None, label="", xscale=None, handles=None, verbose=False,
        debug=False, **kwargs):
        from .myplotspec import get_color
        import numpy as np
        import pandas
        from os.path import expandvars

        # Load data
        if  ("infile" is None):
            return
        if verbose >= 1:
            print("Loading {0}".format(expandvars(infile)))
        load_kw = kwargs.get("load_kw", {}).copy()
        if load_kw.get("sequential_ds", False):
            names = load_kw.pop("names")
            lines_per_block = load_kw.get("lines_per_block", 10000)
            col_i = 0
            data = pandas.DataFrame(columns=names,
              index=range(lines_per_block))
            with open(expandvars(infile)) as f:
                for line in f:
                    if line.startswith("@"):
                        continue
                    else:
                        line = np.array(line.split(), np.float)
                        index = int(line[0])
                        if index not in data.index:
                            data = pandas.concat((data,
                              pandas.DataFrame(columns=names, 
                              index= range(index, index+lines_per_block))))
                        if not np.isnan(data.loc[index][names[col_i]]):
                            col_i += 1
                        data.loc[index][names[col_i]] = line[1]
            data.insert(0, "time", data.index)
            data = data.dropna()
        else:
            data = pandas.read_csv(expandvars(infile), **load_kw)
        if xscale is not None:
            data[xkey] *= xscale

        # Configure plot settings
        plot_kw = kwargs.get("plot_kw", {}).copy()
        if "color" in plot_kw:
            plot_kw["color"] = get_color(plot_kw["color"])
        elif "color" in kwargs:
            plot_kw["color"] = get_color(kwargs.pop("color"))

        # Plot
        handle = subplot.plot(data[xkey], data[ykey], **plot_kw)[0]
        if handles is not None:
            handles[label] = handle
        if subplot2 is not None:
            from sklearn.neighbors import KernelDensity

            bandwidth = kwargs.get("bandwidth", 0.1)
            grid = kwargs.get("grid", np.linspace(0,10,100))
            kde = KernelDensity(bandwidth=bandwidth)
            kde.fit(data[ykey][:, np.newaxis])
            pdf = np.exp(kde.score_samples(grid[:, np.newaxis]))
            pdf /= pdf.sum()
            plot2_kw = plot_kw.copy()
            plot2_kw.update(kwargs.get("plot2_kw", {}))
            subplot2.plot(pdf, grid, **plot2_kw)

#################################### MAIN #####################################
if __name__ == "__main__":
    TimeSeriesFigureManager().main()
