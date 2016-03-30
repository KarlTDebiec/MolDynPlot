#!/usr/bin/python
# -*- coding: utf-8 -*-
#   moldynplot.SAXSFigureManager.py
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
    __package__ = str("moldynplot")
    import moldynplot
from .myplotspec.FigureManager import FigureManager
################################### CLASSES ###################################
class SAXSFigureManager(FigureManager):
    """
    Manages the generation of time series figures
    """

    from .myplotspec.manage_defaults_presets import manage_defaults_presets
    from .myplotspec.manage_kwargs import manage_kwargs

    defaults = """
        draw_figure:
          subplot_kw:
            autoscale_on: False
          multi_tick_params:
            left: on
            right: off
            bottom: on
            top: off
          shared_legend: True
          shared_legend_kw:
            legend_kw:
              frameon: False
              loc: 9
              numpoints: 1
              handletextpad: 0
        draw_subplot:
          title_kw:
            verticalalignment: bottom
          xlabel:      "q ($Å^{-1}$)"
          xticks:      [ 0.00,  0.05,  0.10,  0.15,  0.20,  0.25,  0.30,  0.35]
          xticklabels: ["0.00","0.05","0.10","0.15","0.20","0.25","0.30","0.35"]
          ylabel:      "Intensity"
          yticks:      [ 0.00,  0.02,  0.04,  0.06,  0.08,  0.10,  0.12]
          yticklabels: ["0.00","0.02","0.04","0.06","0.08","0.10","0.12"]
          tick_params:
            direction: out
            left: on
            right: off
            bottom: on
            top: off
            direction: out
          grid: True
          grid_kw:
            b: True
            color: [0.8,0.8,0.8]
            linestyle: '-'
        draw_dataset:
          plot_kw:
            zorder: 10
          handle_kw:
            ls: none
            marker: s
            mec: black
    """

    available_presets = """
      amber:
        class: content
        help: Data from sax_md
        draw_dataset:
          dataset_kw:
            cls: moldynplot.CpptrajDataset.SAXSDataset
            mean: True
      experiment:
        class: content
        help: Data from experiment
        draw_dataset:
          dataset_kw:
            cls: moldynplot.ExperimentDataset.SAXSDataset
            read_csv_kw:
              engine: python
              skiprows: 2
              skipfooter: 5
              names: [q, intensity, intensity_se]
              sep: " "
              skipinitialspace: True
              index_col: 0
      logx:
        class: appearance
        help: Plot x axis using base 10 logarithmic scale
        draw_subplot:
          xlabel: "$log_{10}$(q)"
        draw_dataset:
          logx: True
      logy:
        class: appearance
        help: Plot y axis using base 10 logarithmic scale
        draw_subplot:
          ylabel:      "$log_{10}$(Intensity)"
          yticks:      [-5,-4,-3,-2,-1,0]
          yticklabels: [-5,-4,-3,-2,-1,0]
        draw_dataset:
          logy: True
      presentation:
        class: target
        inherits: presentation
        draw_figure:
          left:       1.20
          sub_width:  7.00
          bottom:     2.90
          sub_height: 3.00
          shared_legend: True
          shared_legend_kw:
            left:       1.20
            sub_width:  7.00
            bottom:     1.70
            sub_height: 0.50
            legend_kw:
              labelspacing: 0.5
              ncol: 2
        draw_dataset:
          partner_kw:
            sub_width: 1.2
            title_fp: 18r
            xlabel_kw:
              labelpad: 20
            label_fp: 18r
            tick_fp: 14r
            xticks:
            lw: 2
            tick_params:
              length: 3
              pad: 6
              width: 2
      manuscript:
        class: target
        inherits: manuscript
        draw_figure:
          left:       0.50
          sub_width:  4.40
          right:      0.20
          bottom:     0.70
          sub_height: 1.80
          top:        0.25
          shared_legend_kw:
            left:       0.50
            sub_width:  4.40
            bottom:     0.00
            sub_height: 0.30
            legend_kw:
              labelspacing: 0.5
              legend_fp: 7r
              ncol: 5
        draw_subplot:
          xlabel_kw:
            labelpad: 3
          ylabel_kw:
            labelpad: 6
        draw_dataset:
          partner_kw:
            hspace:    0.05
            sub_width: 0.8
            title_fp: 8b
            xlabel_kw:
              labelpad: 12.5
            label_fp: 8b
            tick_fp: 6r
            tick_params:
              length: 2
              pad: 3
              width: 1
          plot_kw:
            lw: 1
          handle_kw:
            ms: 6
            mew: 1
      notebook:
        class: target
        inherits: notebook
        draw_figure:
          left:       0.60
          sub_width:  4.40
          right:      0.20
          bottom:     1.00
          sub_height: 1.80
          top:        0.30
          shared_legend: True
          shared_legend_kw:
            left:       0.60
            sub_width:  4.40
            right:      0.20
            bottom:     0.00
            sub_height: 0.50
            legend_kw:
              labelspacing: 0.5
              legend_fp: 8r
              ncol: 2
        draw_dataset:
          plot_kw:
            lw: 1.0
          partner_kw:
            sub_width: 0.8
            title_fp: 10b
            xlabel_kw:
              labelpad:  12.5
            label_fp: 10b
            tick_fp: 8r
            xticks:
            tick_params:
              length: 2
              pad: 6
              width: 1
    """

    @manage_defaults_presets()
    @manage_kwargs()
    def draw_dataset(self, subplot, label=None,
        handles=None, logx=False, logy=False,
        draw_fill_between=False, draw_plot=True,
        verbose=1, debug=0, **kwargs):
        import numpy as np
        import pandas as pd
        from .myplotspec import get_colors, multi_get_copy

        # Load data
        dataset_kw = multi_get_copy("dataset_kw", kwargs, {})
        if "infile" in kwargs:
            dataset_kw["infile"] = kwargs["infile"]
        dataset = self.load_dataset(verbose=verbose, debug=debug, **dataset_kw)
        dataframe = dataset.dataframe

        # Configure plot settings
        plot_kw = multi_get_copy("plot_kw", kwargs, {})
        get_colors(plot_kw, kwargs)

        # Plot fill_between
        if draw_fill_between:
            fill_between_kw = multi_get_copy("fill_between_kw", kwargs, {})
            get_colors(fill_between_kw, plot_kw)
            x = dataframe.index.values
            if logx:
                x = np.log10(x)
            y = dataframe["intensity"]
            y_se = dataframe["intensity_se"]
            if logy:
                y_se = (y_se / (y * np.log(10)))
                y = np.log10(y)
            y_lb = y - 1.96 * y_se
            y_ub = y + 1.96 * y_se
            subplot.fill_between(x, y_lb, y_ub, **fill_between_kw)

        # Plot series
        if draw_plot:
            x = dataframe.index.values
            if logx:
                x = np.log10(x)
            y = dataframe["intensity"]
            if logy:
                y = np.log10(y)
            print(y.values.min(), y.values.max())
            plot = subplot.plot(x, y, **plot_kw)[0]
#            handle_kw = multi_get_copy("handle_kw", kwargs, {})
#            handle_kw["mfc"] = plot.get_color()
#            handle = subplot.plot([-10, -10], [-10, -10], **handle_kw)[0]
#            if handles is not None and label is not None:
#                handles[label] = handle

#################################### MAIN #####################################
if __name__ == "__main__":
    SAXSFigureManager().main()
