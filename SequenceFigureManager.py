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
            autoscale_on: False
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
          yticks:      [0.0,0.5,1.0,1.5,2.0,2.5]
          yticklabels: [0.0,0.5,1.0,1.5,2.0,2.5]
        draw_dataset:
          y_key:    "R1"
          yerr_key: "R1 se"
      R2:
        class: content
        help: Format subplot for R2 relaxation
        draw_subplot:
          xticklabels: []
          ylabel:      "$R_2$"
          yticks:      [0.0,3.0,6.0,9.0,12.0,15.0]
          yticklabels: [0.0,3.0,6.0,9.0,12.0,15.0]
        draw_dataset:
          y_key:    "R2"
          yerr_key: "R2 se"
      HetNOE:
        class: content
        help: Format subplot for Heteronuclear NOE relaxation
        draw_subplot:
          xlabel:      Residue
          ylabel:      "Het\n\nNOE"
          yticks:      [0.0,0.2,0.4,0.6,0.8,1.0]
          yticklabels: [0.0,0.2,0.4,0.6,0.8,1.0]
        draw_dataset:
          y_key:    "NOE"
          yerr_key: "NOE se"
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
          nrows: 3
          subplots:
            0:
              preset: R1
              yticklabels:  [0.0,0.5,1.0,1.5,2.0,2.5]
            1:
              preset: R2
              yticklabels:  [0.0,3.0,6.0,9.0,12.0]
            2:
              preset: HetNOE
              yticklabels:  [0.0,0.2,0.4,0.6,0.8]
      relaxation_4:
        class: content
        help: Four stacked plots including R1, R2, HetNOE, and S2
        draw_figure:
          nrows: 4
          subplots:
            0:
              preset: R1
              yticklabels: [0.0,0.5,1.0,1.5,2.0,2.5]
            1:
              preset: R2
              yticklabels: [0.0,3.0,6.0,9.0,12.0]
            2:
              preset: HetNOE
              yticklabels: [0.0,0.2,0.4,0.6,0.8]
            3:
              preset: order
              yticklabels: [0.0,0.2,0.4,0.6,0.8]
      manuscript:
        class: target
        inherits: manuscript
        draw_figure:
          left:       0.51
          sub_width:  6.31
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
    def draw_dataset(self, subplot, kind, ykey=None, ysekey=None,
        label=None, handles=None, verbose=1, debug=0, **kwargs):
        from .myplotspec import get_color, multi_get_copy
        from .myplotspec.Dataset import Dataset

        # Load data
        dataset_kw = multi_get_copy("dataset_kw", kwargs, {})
        if "infile" in kwargs:
            dataset_kw["infile"] = kwargs["infile"]
        dataframe= self.load_dataset(Dataset, verbose=verbose, debug=debug,
          **dataset_kw).data
        print(dataframe)

        # Configure plot settings
        errorbar_kw = multi_get_copy("errorbar_kw", kwargs, {})
        if "color" in errorbar_kw:
            errorbar_kw["color"] = get_color(errorbar_kw["color"])
        elif "color" in kwargs:
            errorbar_kw["color"] = get_color(kwargs.pop("color"))

#        # Plot data
        x = [int(c.split(":")[1]) for c in list(dataframe.columns.values)]
        print(ykey)
        subplot.errorbar(x, y=dataframe.loc[ykey], yerr=dataframe.loc[ysekey],
          **errorbar_kw)

#################################### MAIN #####################################
if __name__ == "__main__":
    SequenceFigureManager().main()
