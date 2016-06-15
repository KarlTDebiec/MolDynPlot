#!/usr/bin/python
# -*- coding: utf-8 -*-
#   moldynplot.CorrFigureManager.py
#
#   Copyright (C) 2015-2016 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Generates one or more correlation figures to specifications in a YAML file.
"""
################################### MODULES ###################################
from __future__ import absolute_import,division,print_function,unicode_literals
if __name__ == "__main__":
    __package__ = str("moldynplot")
    import moldynplot
from .myplotspec.FigureManager import FigureManager
################################### CLASSES ###################################
class CorrFigureManager(FigureManager):
    """
    Manages the generation of correlation figures.
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
            handle_kw:
              marker: s
              ls: none
              mec: black
            legend_kw:
              frameon: False
              loc: 9
              numpoints: 1
              handletextpad: 0.0
        draw_subplot:
          title_kw:
            verticalalignment: bottom
          tick_params:
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
          dataset_kw:
            cls: moldynplot.CorrDataset.CorrDataset
            x_kw:
              read_csv_kw:
                delim_whitespace: True
                index_col: 0
            y_kw:
              read_csv_kw:
                delim_whitespace: True
                index_col: 0
          corr_line_kw:
            x: [0, 100]
            y: [0, 100]
            color: [0.5,0.5,0.5]
            zorder: 9
          errorbar_kw:
            zorder: 10
          plot_kw:
            zorder: 10
          handle_kw:
            marker: s
            ls: none
            mec: black
    """

    available_presets = """
      r1:
        class: content
        help: Format subplot for R1 relaxation
        draw_subplot:
          xlabel: "$R_1$ Experiment"
          xticks: [0,1,2,3,4,5]
          yticks: [0,1,2,3,4,5]
          ylabel: "$R_1$ Simulation"
        draw_dataset:
          xkey:   r1
          xsekey: r1_se
          ykey:   r1
          ysekey: r1_se
      r2:
        class: content
        help: Format subplot for R2 relaxation
        draw_subplot:
          xlabel: "$R_2$ Experiment"
          xticks: [0,2,4,6,8,10,12,14,16,18,20]
          ylabel: "$R_2$ Simulation"
          yticks: [0,2,4,6,8,10,12,14,16,18,20]
        draw_dataset:
          xkey:   r2
          xsekey: r2_se
          ykey:   r2
          ysekey: r2_se
      hetnoe:
        class: content
        help: Format subplot for Heteronuclear NOE relaxation
        draw_subplot:
          xlabel: "HetNOE Experiment"
          xticks: [0.0,0.2,0.4,0.6,0.8,1.0]
          ylabel: "HetNOE Simulation"
          yticks: [0.0,0.2,0.4,0.6,0.8,1.0]
        draw_dataset:
          xkey:   noe
          xsekey: noe_se
          ykey:   noe
          ysekey: noe_se
      s2:
        class: content
        help: Format subplot for S2 order parameter
        draw_subplot:
          xlabel: "$S^2$ Experiment"
          xticks: [0.0,0.2,0.4,0.6,0.8,1.0]
          ylabel: "$S^2$ Simulation"
          yticks: [0.0,0.2,0.4,0.6,0.8,1.0]
        draw_dataset:
          xkey:   s2
          xsekey: s2_se
          ykey:   s2
          ysekey: s2_se
      rotdif:
        class: content
        help: Format subplot for rotational diffusion
        draw_figure:
          multi_xticklabels:  [0,1,2,3,4,5,6,7,8]
          multi_yticklabels:  [0,1,2,3,4,5,6,7,8]
        draw_subplot:
          xlabel: $\\tau_c$ Experiment
          xticks: [0,1,2,3,4,5,6,7,8]
          ylabel: $\\tau_c$ Simulation
          yticks: [0,1,2,3,4,5,6,7,8]
      relaxation_4:
        help: Four plots (2x2) including R1, R2, HetNOE, and S2
        draw_figure:
          nrows: 2
          ncols: 2
          subplots:
            0:
              preset: r1
            1:
              preset: r2
            2:
              preset: hetnoe
            3:
              preset: s2
      manuscript:
        class: target
        inherits: manuscript
        draw_figure:
          left:       0.40
          sub_width:  2.50
          wspace:     0.10
          right:      0.20
          bottom:     0.70
          sub_height: 2.50
          hspace:     0.10
          top:        0.20
          shared_legend: True
          shared_legend_kw:
            left:       0.40
            sub_width:  2.50
            sub_height: 0.30
            bottom:     0.00
            legend_kw:
              labelspacing: 0.5
              legend_fp: 7r
              ncol: 5
        draw_dataset:
          errorbar_kw:
            capsize: 2
            capthick: 1
            elinewidth: 1
          handle_kw:
            ms: 6
            mew: 1
          label_kw:
            fp: 7r
            text_kw:
              ha: center
              va: center
            border_lw: 1
      presentation:
        class: target
        inherits: presentation
        draw_figure:
          fig_width:  5.12
          left:       2.00
          sub_width:  1.944444
          bottom:     2.20
          sub_height: 3.50
          shared_legend: True
          shared_legend_kw:
            left:       2.00
            sub_width:  1.944444
            sub_height: 1.5
            bottom:     0.00
            legend_kw:
              labelspacing: 0.5
              ncol: 5
        draw_dataset:
          errorbar_kw:
            capsize: 6
            capthick: 2
            elinewidth: 2
          handle_kw:
            ms: 12
            mew: 2
          label_kw:
            fp: 14r
            text_kw:
              ha: center
              va: center
            border_lw: 2
    """

    @manage_defaults_presets()
    @manage_kwargs()
    def draw_dataset(self, subplot,
        xkey=None, ykey=None, xsekey=None, ysekey=None, label=None,
        handles=None,
        draw_corr_line=True, draw_errorbar=True, draw_plot=True,
        draw_handle=False, draw_label=False,
        verbose=1, debug=0, **kwargs):
        from os.path import expandvars
        from warnings import warn
        import numpy as np
        import pandas as pd
        from .myplotspec import get_colors, multi_get_copy

        # Load data
        dataset_kw = multi_get_copy("dataset_kw", kwargs, {})
        if "cls" in dataset_kw and dataset_kw["cls"] is not None:
            dataset = self.load_dataset(verbose=verbose, debug=debug,
              **dataset_kw)
            dataframe = dataset.dataframe
        else:
            dataset = dataframe = None

        # Configure plot settings
        plot_kw = multi_get_copy("plot_kw", kwargs, {})
        get_colors(plot_kw, kwargs)

        # Plot expected correlation line
        if draw_corr_line and not hasattr(subplot, "_mps_corr_line"):
            corr_line_kw = multi_get_copy("corr_line_kw", kwargs, {})
            get_colors(corr_line_kw)
            # Must pass 'x' and 'y' as positional arguments
            subplot._mps_corr_line = subplot.plot(
              corr_line_kw.pop("x", [0,10]), corr_line_kw.pop("y", [0,10]),
              **corr_line_kw)

        # Plot error bar
        if draw_errorbar:
            errorbar_kw = multi_get_copy("errorbar_kw", kwargs, {})
            get_colors(errorbar_kw, plot_kw)

            if "x" in errorbar_kw:
                eb_x = errorbar_kw.pop("x")
            elif dataframe is not None:
                eb_x = dataframe[xkey, "x"]
            if "y" in errorbar_kw:
                eb_y = errorbar_kw.pop("y")
            elif dataframe is not None:
                eb_y = dataframe[ykey, "y"]
            if "xerr" not in errorbar_kw:
                if "xse" in errorbar_kw:
                    errorbar_kw["xerr"] = errorbar_kw.pop("xse") * 1.96
                elif "xse" in kwargs:
                    errorbar_kw["xerr"] = kwargs.pop("xse") * 1.96
                elif xsekey is not None and (xsekey, "x") in dataframe:
                    errorbar_kw["xerr"] = errorbar_kw.get("xse",
                      dataframe[xsekey, "x"]) * 1.96
            if "yerr" not in errorbar_kw:
                if "yse" in errorbar_kw:
                    errorbar_kw["yerr"] = errorbar_kw.pop("yse") * 1.96
                elif "yse" in kwargs:
                    errorbar_kw["yerr"] = kwargs.pop("yse") * 1.96
                elif ysekey is not None and (ysekey, "y") in dataframe:
                    errorbar_kw["yerr"] = errorbar_kw.get("yse",
                      dataframe[ysekey, "y"]) * 1.96
            subplot.errorbar(eb_x, eb_y, **errorbar_kw)

        # Plot legend handles
        if draw_handle:
            if label is None:
                warn("'draw_handle' is enabled but a value for argument"
                     "'label' was not provided; skipping")
            else:
                handle_kw = multi_get_copy("handle_kw", kwargs, {})
                get_colors(handle_kw, plot_kw)
                handle = subplot.plot(
                  handle_kw.pop("x", [-10,-10]), handle_kw.pop("y", [-10,-10]),
                  **handle_kw)[0]
                if handles is not None:
                    handles[label] = handle

        # Draw label
        if draw_label:
            from .myplotspec.text import set_text

            label_kw = multi_get_copy("label_kw", kwargs, {})
            set_text(subplot, **label_kw)

#################################### MAIN #####################################
if __name__ == "__main__":
    CorrFigureManager().main()
