#!/usr/bin/python
# -*- coding: utf-8 -*-
#   moldynplot.SSFigureManager.py
#
#   Copyright (C) 2015 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Generates one or more secondary structure figures to specifications in a
YAML file.
"""
################################### MODULES ###################################
from __future__ import absolute_import,division,print_function,unicode_literals
if __name__ == "__main__":
    __package__ = str("moldynplot")
    import moldynplot
from .myplotspec.FigureManager import FigureManager
################################### CLASSES ###################################
class SSFigureManager(FigureManager):
    """
    Manages the generation of secondary structure figures.
    """

    from .myplotspec.manage_defaults_presets import manage_defaults_presets
    from .myplotspec.manage_kwargs import manage_kwargs

    defaults = """
        draw_figure:
          subplot_kw:
            autoscale_on: False
        draw_subplot:
          xlabel: Time (µs)
          ylabel_kw:
            va: center
          tick_params:
            left: off
            right: off
            bottom: off
            top: off
          y2ticks: []
          y2label_kw:
            rotation: 270
          grid: True
          grid_kw:
            b: True
            linestyle: '-'
            axis: x
        draw_dataset:
          heatmap_kw:
            cmap: !!python/object/apply:moldynplot.dssp_color_palette []
            edgecolors: none
            rasterized: True
    """

    available_presets = """
      dssp:
        class: content
        help: Dynamic secondary structure of proteins calculated by cpptraj
        draw_dataset:
          labels:
            0: None
            1: Parallel β Sheet
            2: Antiparallel β Sheet
            3: $3_{10}$ Helix
            4: α Helix
            5: π Helix
            6: Turn
            7: Bend
          heatmap_kw:
            vmin: 0
            vmax: 7
      notebook:
        class: target
        inherits: notebook
        draw_figure:
          left:       0.50
          sub_width:  5.00
          right:      0.25
          bottom:     1.00
          sub_height: 2.00
          hspace:     0.10
          top:        0.25
          title_kw:
            top: -0.1
          shared_legend: True
          shared_legend_kw:
            left:       0.50
            sub_width:  5.00
            right:      0.25
            bottom:     0.00
            sub_height: 0.50
            legend_kw:
              loc: 9
              numpoints: 1
              handletextpad: 0
              ncol: 4
        draw_subplot:
          y2label_kw:
            labelpad: 12
    """

    @manage_defaults_presets()
    @manage_kwargs()
    def draw_dataset(self, subplot, downsample=None, label=None, handles=None,
        verbose=1, debug=0, **kwargs):
        """
        """
        from os.path import expandvars
        import numpy as np
        import pandas as pd
        from .myplotspec import multi_get_copy
        from sys import exit

        # Load data
        infile = expandvars(kwargs.get("infile"))
        dataset = pd.read_csv(infile, delim_whitespace=True,
          index_col=0)

        # Scale:
        dt = kwargs.get("dt", 0.001)
        dataset.index *= dt
        dataset.index.names = ["time"]

        # Downsample
        if downsample is not None:
            full_size = dataset.shape[0]
            reduced_size = int(full_size / downsample)
            reduced = pd.DataFrame(0.0, index=range(0, reduced_size),
              columns=dataset.columns, dtype=np.int64)
            for i in range(1, reduced_size):
                reduced.loc[i] = dataset[
                  i*downsample:(i+1)*downsample].mode().loc[0]
            reduced.index *= (dt * downsample) / 1000
            reduced.index.names = ["time"]
            dataset = reduced

        heatmap = True
        if heatmap:
            heatmap_kw = multi_get_copy("heatmap_kw", kwargs, {})
            pcolormesh = subplot.pcolor(
              dataset.index.values,
              np.arange(len(dataset.columns)),
              dataset.T.values,
            **heatmap_kw)

        if handles is not None:
            cmap = pcolormesh.cmap
            labels = kwargs.get("labels")
            h_kw=dict(ls="none", marker="s", ms=5, mec="k")

            for i in sorted(labels.keys()):
                label = labels[i]
                handles[label] = subplot.plot((0),(0),
                  mfc=cmap(i / (len(labels) - 1)),
                  **h_kw)[0]

#################################### MAIN #####################################
if __name__ == "__main__":
    SSFigureManager().main()
