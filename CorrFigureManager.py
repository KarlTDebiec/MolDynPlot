#!/usr/bin/python
# -*- coding: utf-8 -*-
#   moldynplot.CorrFigureManager.py
#
#   Copyright (C) 2015 Karl T Debiec
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
          corr_line_kw:
            x: [0, 10]
            y: [0, 10]
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
      R1:
        class: content
        help: Format subplot for R1 relaxation
        draw_subplot:
          xlabel:      $R_1$ Experiment
          xticks:      [0.0,0.5,1.0,1.5,2.0,2.5,3.0]
          xticklabels: [0.0,0.5,1.0,1.5,2.0,2.5,3.0]
          ylabel:      $R_1$ Simulation
          yticks:      [0.0,0.5,1.0,1.5,2.0,2.5,3.0]
          yticklabels: [0.0,0.5,1.0,1.5,2.0,2.5,3.0]
        draw_dataset:
          y_key:    R1
          yerr_key: R1 se
      R2:
        class: content
        help: Format subplot for R2 relaxation
        draw_subplot:
          xlabel:      $R_2$ Experiment
          xticks:      [0.0,3.0,6.0,9.0,12.0,15.0]
          xticklabels: [0.0,3.0,6.0,9.0,12.0,15.0]
          ylabel:      $R_2$ Simulation
          yticks:      [0.0,3.0,6.0,9.0,12.0,15.0]
          yticklabels: [0.0,3.0,6.0,9.0,12.0,15.0]
        draw_dataset:
          y_key:    R2
          yerr_key: R2 se
      HetNOE:
        class: content
        help: Format subplot for Heteronuclear NOE relaxation
        draw_subplot:
          xlabel:      "HetNOE Experiment"
          xticks:      [0.0,0.2,0.4,0.6,0.8,1.0]
          xticklabels: [0.0,0.2,0.4,0.6,0.8,1.0]
          ylabel:      "HetNOE Simulation"
          yticks:      [0.0,0.2,0.4,0.6,0.8,1.0]
          yticklabels: [0.0,0.2,0.4,0.6,0.8,1.0]
        draw_dataset:
          y_key:    "NOE"
          yerr_key: "NOE se"
      order:
        class: content
        help: Format subplot for S2 order parameter
        draw_subplot:
          xlabel:      $S^2$ Experiment
          xticks:      [0.0,0.2,0.4,0.6,0.8,1.0]
          xticklabels: [0.0,0.2,0.4,0.6,0.8,1.0]
          ylabel:      $S^2$ Simulation
          yticks:      [0.0,0.2,0.4,0.6,0.8,1.0]
          yticklabels: [0.0,0.2,0.4,0.6,0.8,1.0]
        draw_dataset:
          y_key:    "S2"
          yerr_key: "S2 se"
      rotdif:
        class: content
        help: Format subplot for rotational diffusion
        draw_figure:
          multi_xticklabels:  [0,1,2,3,4,5,6,7,8]
          multi_yticklabels:  [0,1,2,3,4,5,6,7,8]
        draw_subplot:
          xlabel:      $\\tau_c$ Experiment
          xticks:      [0,1,2,3,4,5,6,7,8]
          ylabel:      $\\tau_c$ Simulation
          yticks:      [0,1,2,3,4,5,6,7,8]
      relaxation_4:
        help: Four plots (2x2) including R1, R2, HetNOE, and S2
        draw_figure:
          nrows: 2
          ncols: 2
          subplots:
            0:
              preset: R1
            1:
              preset: R2
            2:
              preset: HetNOE
            3:
              preset: order
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
    """

    @manage_defaults_presets()
    @manage_kwargs()
    def draw_dataset(self, subplot,
        label=None,
        handles=None,
        draw_corr_line=True, draw_errorbar=True, draw_plot=True,
        draw_handle=False, draw_label=False,
        verbose=1, debug=0, **kwargs):
        from os.path import expandvars
        from warnings import warn
        import numpy as np
        import pandas as pd
        from .myplotspec import get_colors, multi_get_copy

        # Configure plot settings
        plot_kw = multi_get_copy("plot_kw", kwargs, {})
        get_colors(plot_kw, kwargs)

        # Plot x=y correlation line
        if draw_corr_line and not hasattr(subplot, "_mps_corr_line"):
            corr_line_kw = multi_get_copy("corr_line_kw", kwargs, {})
            get_colors(corr_line_kw)
            subplot._mps_corr_line = subplot.plot(
              corr_line_kw.pop("x", [0,10]), corr_line_kw.pop("y", [0,10]),
              **corr_line_kw)

        # Plot error bar
        if draw_errorbar:
            errorbar_kw = multi_get_copy("errorbar_kw", kwargs, {})
            get_colors(errorbar_kw, plot_kw)
            subplot.errorbar(
              errorbar_kw.pop("x", [0,10]), errorbar_kw.pop("y", [0,10]),
              yerr=errorbar_kw.pop("yerr")*1.96,
              **errorbar_kw)

        # Plot manual legend handles
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
            print(label_kw)
            set_text(subplot, **label_kw)

#################################### MAIN #####################################
if __name__ == "__main__":
    CorrFigureManager().main()
