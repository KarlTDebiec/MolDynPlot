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
#        draw_figure:
#          subplot_kw:
#            autoscale_on: False
         draw_subplot:
          xlabel: Residue
         draw_dataset:
           bar_kw:
             align: center
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
        help: Three stacked plots including R1, R2, and HetNOE
        draw_figure:
          nrows:        3
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
        from .SSDataset import SSDataset
        dataset_classes = {"ss": SSDataset}

        # Load data
        kind = kind.lower()
        dataset_kw = multi_get_copy("dataset_kw", kwargs, {})
        if "infile" in kwargs:
            dataset_kw["infile"] = kwargs["infile"]
        try:
            dataset = self.load_dataset(dataset_classes[kind],
                        dataset_classes=dataset_classes,
                        verbose=verbose, debug=debug, **dataset_kw).data
        except TypeError:
            from warnings import warn
            warn("{0} has raised an ".format(dataset_classes[kind].__name__) +
              "error; skipping this dataset.")
            return

        # Configure plot settings
        bar_kw = multi_get_copy("bar_kw", kwargs, {})
        for color_key in ["color", "edgecolor", "facecolor", "ecolor"]:
            if color_key in bar_kw:
                bar_kw[color_key] = get_color(bar_kw[color_key])
        if label is not None:
            bar_kw["label"] = label

        # Plot data
        subplot.bar(dataset.index - 1, dataset[ykey], **bar_kw)

#################################### MAIN #####################################
if __name__ == "__main__":
    SequenceFigureManager().main()
