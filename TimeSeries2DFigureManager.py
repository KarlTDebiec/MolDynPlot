#!/usr/bin/python
# -*- coding: utf-8 -*-
#   moldynplot.TimeSeries2DFigureManager.py
#
#   Copyright (C) 2015 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Generates one or more 2D time series figures to specifications in a YAML
file.
"""
################################### MODULES ###################################
from __future__ import absolute_import,division,print_function,unicode_literals
if __name__ == "__main__":
    __package__ = str("moldynplot")
    import moldynplot
from .myplotspec.FigureManager import FigureManager
################################### CLASSES ###################################
class TimeSeries2DFigureManager(FigureManager):
    """
    Manages the generation of 2D time series figures.
    """

    from .myplotspec.manage_defaults_presets import manage_defaults_presets
    from .myplotspec.manage_kwargs import manage_kwargs

    defaults = """
        draw_figure:
          subplot_kw:
            autoscale_on: False
          shared_legend_kw:
            legend_kw:
              frameon: False
              loc: 9
              numpoints: 1
              handletextpad: 0
        draw_subplot:
          xlabel: Time (µs)
          ylabel: Residue
          ylabel_kw:
            va: center
          tick_params:
            left: on
            right: off
            bottom: on
            top: off
            direction: out
          y2ticks: []
          y2label_kw:
            rotation: 270
          grid: True
          grid_kw:
            b: True
            linestyle: '-'
            axis: both
            color: [0.2,0.2,0.2]
        draw_dataset:
          heatmap: True
          heatmap_kw:
            edgecolors: none
            rasterized: True
    """

    available_presets = """
      dssp:
        class: content
        help: Dynamic secondary structure of proteins calculated by cpptraj
        draw_subplot:
          ylabel: Residue
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
            cmap: !!python/object/apply:moldynplot.dssp_color_palette []
            vmin: 0
            vmax: 7
      perres_rmsd:
        class: content
        help: Per-residue RMSD calculated by cpptraj
        draw_dataset:
          heatmap_kw:
            cmap: afmhot
            vmin: 0
            vmax: 10
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
              ncol: 4
        draw_subplot:
          y2label_kw:
            labelpad: 12
    """

    @manage_defaults_presets()
    @manage_kwargs()
    def draw_dataset(self, subplot, downsample=None, label=None, heatmap=True,
        handles=None, verbose=1, debug=0, **kwargs):
        """
        """
        from os.path import expandvars
        import numpy as np
        import pandas as pd
        from .myplotspec import multi_get_copy

        # Load data
        print("loading")
        infile = expandvars(kwargs.get("infile"))
        dataset = pd.read_csv(infile, delim_whitespace=True, index_col=0)
        print("loaded")

        # Scale:
        dt = kwargs.get("dt", 0.001)
        dataset.index *= dt
        dataset.index.names = ["time"]

        # Downsample
#        if downsample is not None:
#            full_size = dataset.shape[0]
#            reduced_size = int(full_size / downsample)
#            reduced = pd.DataFrame(0.0, index=range(0, reduced_size),
#              columns=dataset.columns, dtype=np.int64)
#            for i in range(1, reduced_size):
#                reduced.iloc[i] = dataset.iloc[
#                  i*downsample:(i+1)*downsample].mode().iloc[0]
#            reduced.index *= (dt * downsample) / 1000
#            reduced.index.names = ["time"]
#            dataset = reduced

        if heatmap:
            heatmap_kw = multi_get_copy("heatmap_kw", kwargs, {})
            pcolormesh = subplot.pcolor(
              dataset.index.values,
              np.arange(0, len(dataset.columns) + 1, 1),
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
    TimeSeries2DFigureManager().main()
