#!/usr/bin/python
# -*- coding: utf-8 -*-
#   moldynplot.SequenceFigureManager.py
#
#   Copyright (C) 2015 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Generates one or more sequence figures to specifications in a YAML file.
"""
################################### MODULES ###################################
from __future__ import absolute_import,division,print_function,unicode_literals
import matplotlib
if __name__ == "__main__":
    __package__ = str("moldynplot")
    import moldynplot
from .myplotspec.FigureManager import FigureManager
################################### CLASSES ###################################
class SequenceFigureManager(FigureManager):
    """
    Manages the generation of sequence figures.
    """

    from .myplotspec.manage_defaults_presets import manage_defaults_presets
    from .myplotspec.manage_kwargs import manage_kwargs

    defaults = """
        draw_figure:
          subplot_kw:
            autoscale_on: True
          multi_tick_params:
            left: on
            right: off
            bottom: on
            top: off
        draw_subplot:
          title_kw:
            verticalalignment: bottom
          xlabel: Residue
          ylabel_kw:
            va: center
          tick_params:
            bottom: on
            top: off
            right: off
            left: on
            direction: out
          grid: True
          grid_kw:
            b: True
            linestyle: '-'
            color: [0.8,0.8,0.8]
        draw_dataset:
          dataset_kw:
            read_csv_kw:
              delim_whitespace: True
              index_col: 0
          errorbar_kw:
            ls: None
            zorder: 10
    """
    available_presets = """
      secstruct:
        class: content
        help: Plot secondary structure as a function of sequence
        draw_subplot:
          ylabel: '% α Helix'
          yticks: [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
          yticklabels: [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
        draw_dataset:
          kind: ss
          ykey: Alpha
      R1:
        class: content
        help: Format subplot for R1 relaxation
        draw_subplot:
          xticklabels: []
          ylabel:      "$R_1$"
          yticks:      [0,1,2,3,4,5]
        draw_dataset:
          ykey:   r1
          ysekey: r1_se
      R2:
        class: content
        help: Format subplot for R2 relaxation
        draw_subplot:
          xticklabels: []
          ylabel:      "$R_2$"
          yticks:      [0,2,4,6,8,10,12,14,16,18,20]
        draw_dataset:
          ykey:   r2
          ysekey: r2_se
      HetNOE:
        class: content
        help: Format subplot for Heteronuclear NOE relaxation
        draw_subplot:
          xlabel:      Residue
          ylabel:      "Het\n\nNOE"
          yticks:      [0.0,0.2,0.4,0.6,0.8,1.0]
        draw_dataset:
          ykey:   noe
          ysekey: noe_se
      order:
        class: content
        help: Format subplot for S2 order parameter
        draw_subplot:
          xlabel:      Residue
          ylabel:      "$S^2$"
          yticks:      [0.0,0.2,0.4,0.6,0.8,1.0]
          yticklabels: [0.0,0.2,0.4,0.6,0.8,1.0]
        draw_dataset:
          y_key:    "S2"
          yerr_key: "S2 se"
      relaxation_3:
        class: content
        help: Three stacked plots including R1, R2, and HetNOE
        draw_figure:
          multiplot: True
          nrows: 3
          subplots:
            0:
              preset: R1
            1:
              preset: R2
            2:
              preset: HetNOE
      relaxation_4:
        class: content
        help: Four stacked plots including R1, R2, HetNOE, and S2
        draw_figure:
          nrows: 4
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
          left:       0.60
          sub_width:  6.00
          wspace:     0.05
          right:      0.12
          bottom:     0.35
          sub_height: 1.00
          hspace:     0.00
          top:        0.25
          multiplot: True
          title_kw:
            top: -0.1
        draw_subplot:
          legend: False
          xlabel_kw:
            labelpad: 3
          ylabel_kw:
            rotation: horizontal
            labelpad: 3
          grid_kw:
            alpha: 0.3
        draw_dataset:
          errorbar_kw:
            lw: 1
            capthick: 1
      notebook:
        help: Notebook (width ≤ 6.5, height ≤ 9)
        inherits: notebook
        draw_figure:
          left:       0.70
          sub_width:  5.60
          right:      0.20
          bottom:     0.40
          sub_height: 1.00
          hspace:     0.00
          top:        0.30
        draw_subplot:
          ylabel_kw:
            labelpad: 15
        draw_dataset:
          bar_kw:
            lw: 0.5
      presentation:
        help: Three stacked plots for presentation (width = 10.24,
              height = 7.68)
        inherits: presentation
        draw_figure:
          left:       1.50
          sub_width:  7.50
          bottom:     3.60
          sub_height: 2.60
        draw_subplot:
          legend: False
    """

    @manage_defaults_presets()
    @manage_kwargs()
    def draw_dataset(self, subplot,
        ykey=None, ysekey=None, label=None,
        handles=None,
        draw_fill_between=False, draw_errorbar=True,
        verbose=1, debug=0, **kwargs):
        import numpy as np
        from .myplotspec import get_colors, multi_get_copy
        from .myplotspec.Dataset import Dataset

        # Load data
        dataset_kw = multi_get_copy("dataset_kw", kwargs, {})
        if "infile" in kwargs:
            dataset_kw["infile"] = kwargs["infile"]
        dataframe= self.load_dataset(Dataset, verbose=verbose, debug=debug,
          **dataset_kw).dataframe
        x = np.array([filter(lambda x: x in '0123456789.', s)
              for s in dataframe.index.values], np.int)

        # Configure plot settings
        plot_kw = multi_get_copy("plot_kw", kwargs, {})
        get_colors(plot_kw, kwargs)

        # Plot fill_between
        if draw_fill_between:
            fill_between_kw = multi_get_copy("fill_between_kw", kwargs, {})
            get_colors(fill_between_kw, plot_kw)
            fb_x, fb_lb, fb_ub = [], [], []
            yse_min = kwargs.get("yse_min")
            for residue in range(x.min(), x.max()+1):
                if residue in x:
                    fb_x.extend([residue-0.5, residue+0.5])
                    index = np.argmax(x==residue)
                    if ysekey is not None:
                        yse = dataframe[ysekey][index]
                        if yse_min is not None and yse < yse_min:
                            yse = yse_min
                    elif yse_min is not None:
                        yse = yse_min
                    fb_lb.extend([
                      dataframe[ykey][index] - yse * 1.96,
                      dataframe[ykey][index] - yse * 1.96])
                    fb_ub.extend([
                      dataframe[ykey][index] + yse * 1.96,
                      dataframe[ykey][index] + yse * 1.96])
                else:
                    fb_x.append(None)
                    fb_lb.append(None)
                    fb_ub.append(None)
#                fb_x.append(None)
#                fb_lb.append(None)
#                fb_ub.append(None)
            fb_x = np.array(fb_x, np.float)
            fb_lb = np.array(fb_lb, np.float)
            fb_ub = np.array(fb_ub, np.float)
                
            subplot.fill_between(fb_x, fb_lb, fb_ub, **fill_between_kw)

        # Plot error bar
        if draw_errorbar:
            errorbar_kw = multi_get_copy("errorbar_kw", kwargs, {})
            get_colors(errorbar_kw, plot_kw)
            subplot.errorbar(x, y=dataframe[ykey], yerr=dataframe[ysekey]*1.96,
              **errorbar_kw)

#################################### MAIN #####################################
if __name__ == "__main__":
    SequenceFigureManager().main()
