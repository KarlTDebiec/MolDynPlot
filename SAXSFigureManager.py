#!/usr/bin/python
# -*- coding: utf-8 -*-
#   moldynplot.SAXSFigureManager.py
#
#   Copyright (C) 2015-2016 Karl T Debiec
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
          multi_xticklabels: ["0.00","0.05","0.10","0.15",
                              "0.20","0.25","0.30","0.35"]
          multi_yticklabels: ["0.00","0.02","0.04","0.06",
                              "0.08","0.10","0.12"]
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
          label_kw:
            zorder: 10
            horizontalalignment: right
            verticalalignment: top
        draw_dataset:
          plot_kw:
            zorder: 10
          handle_kw:
            ls: none
            marker: s
            mec: black
    """

    available_presets = """
      simulation:
        class: content
        help: Data from simulation (saxs_md, crysol, foxs)
        draw_dataset:
          dataset_kw:
            cls: moldynplot.Dataset.SAXSTimeSeriesDataset
            calc_mean: True
            calc_error: True
      experiment:
        class: content
        help: Data from experiment
        draw_dataset:
          dataset_kw:
            cls: moldynplot.Dataset.SAXSExperimentDataset
            read_csv_kw:
              engine: python
              skiprows: 2
              skipfooter: 5
              names: [q, intensity, intensity_se]
              sep: " "
              skipinitialspace: True
              index_col: 0
      envelope:
        class: content
        help: Data back-calculated from molecular envelope
        draw_dataset:
          dataset_kw:
            cls: moldynplot.Dataset.SAXSExperimentDataset
            read_csv_kw:
              engine: python
              skiprows: 1
              names: [q, intensity, something]
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
        draw_figure:
          multi_yticklabels: [-4,-3,-2,-1,0]
        draw_subplot:
          ylabel: "$log_{10}$(Intensity)"
          yticks: [-4,-3,-2,-1,0]
        draw_dataset:
          logy: True
      kratky:
        class: appearance
        help: Plot y axis using intensity*q^2
        draw_figure:
          subplot_kw:
            autoscale_on: False
          multi_yticklabels: null
          multi_yticklabels: [0, 5, 10, 15, 20, 25, 30]
        draw_subplot:
          ylabel: "Intensity*$q^2$ ($x10^{-5}$)"
          yticks: [0, 0.00005, 0.00010, 0.00015, 0.00020, 0.00025, 0.00030]
          yticklabels: [0, 5, 10, 15, 20, 25, 30]
        draw_dataset:
          kratky: True
      diffy:
        class: appearance
        help: Plot difference on y axis
        draw_figure:
          multi_yticklabels: ["-0.020","-0.015","-0.010","-0.005","0.000",
                               "0.005", "0.010", "0.015", "0.020"]
        draw_subplot:
          ylabel: "Δ Intensity"
          yticks:      [ -0.020,  -0.015,  -0.010,  -0.005,  0.000,
                          0.005,   0.010,   0.015,   0.020]
          yticklabels: ["-0.020","-0.015","-0.010","-0.005","0.000",
                         "0.005", "0.010", "0.015", "0.020"]
        draw_dataset:
          dataset_kw:
            cls: moldynplot.Dataset.SAXSDiffDataset
      presentation:
        class: target
        inherits: presentation
        draw_figure:
          left:       1.25
          sub_width:  7.00
          bottom:     2.90
          sub_height: 3.00
          shared_legend: True
          shared_legend_kw:
            left:       8.35
            sub_width:  2.00
            bottom:     2.90
            sub_height: 3.00
            legend_kw:
              labelspacing: 0.5
              loc: 2
              ncol: 1
        draw_dataset:
          handle_kw:
            ms: 12
            mew: 2
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
          title_kw:
            top: -0.125
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
          draw_label: True
          label_kw:
            border_lw: 1
            xabs: -0.05
            yabs: -0.05
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
    """

    @manage_defaults_presets()
    @manage_kwargs()
    def draw_dataset(self, subplot, label=None, handles=None, logx=False,
        logy=False, kratky=False, draw_fill_between=False, draw_plot=True,
        draw_handle=True, draw_label=True, verbose=1, debug=0, **kwargs):
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
                y    = np.log10(y)
            elif kratky:
                y_se = y_se*x*x
                y    = y*x*x
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
            elif kratky:
                y = y*x*x
            plot = subplot.plot(x, y, **plot_kw)[0]

        if draw_handle:
            handle_kw = multi_get_copy("handle_kw", kwargs, {})
            handle_kw["mfc"] = plot_kw.get("color", "blue")
            handle = subplot.plot([-10, -10], [-10, -10], **handle_kw)[0]
            if handles is not None and label is not None:
                handles[label] = handle

#################################### MAIN #####################################
if __name__ == "__main__":
    SAXSFigureManager().main()
