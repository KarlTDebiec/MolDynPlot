#!/usr/bin/python
# -*- coding: utf-8 -*-
#   moldynplot.StateProbFigureManager.py
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
        draw_subplot:
          ylabel_kw:
            va: center
          y2label_kw:
            rotation: 270
          grid: True
          grid_kw:
            b: True
            linestyle: '-'
        draw_dataset:
          heatmap_kw:
            cmap: cubehelix_r
            edgecolors: none
            rasterized: True
            vmin: 0
            vmax: 7
    """

    available_presets = """
      notebook:
        class: target
        inherits: notebook
        draw_figure:
          left:       0.50
          sub_width:  5.00
          right:      0.25
          bottom:     0.70
          sub_height: 2.00
          hspace:     0.10
          top:        0.25
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
    """

    @manage_defaults_presets()
    @manage_kwargs()
    def draw_dataset(self, subplot, label=None, handles=None, verbose=1,
        debug=0, **kwargs):
        from os.path import expandvars
        import numpy as np
        import pandas
        from .myplotspec import multi_get_copy
        from sys import exit

        infile = expandvars(kwargs.get("infile"))
        dataset = pandas.read_csv(infile, delim_whitespace=True,
          index_col=0, dtype=np.int64)

        # Downsample
        downsample = kwargs.get("downsample", 1000)
        if downsample is not None:
            full_size = dataset.shape[0]
            reduced_size = int(full_size / downsample)
            reduced = pandas.DataFrame(0.0, index=range(1, reduced_size + 1),
              columns = dataset.columns, dtype=np.int64)
            for i in range(1, reduced_size):
                reduced.loc[i] = dataset[
                  i*downsample:(i+1)*downsample].mode().loc[0]
        dataset = reduced

        heatmap = True
        if heatmap:
            heatmap_kw = multi_get_copy("heatmap_kw", kwargs, {})
            pcolormesh = subplot.pcolor(dataset.transpose(),
            **heatmap_kw)
        if handles is not None:
            h_ar=((-10,-10), (-10,-10))
            h_kw=dict(ls="none", marker="s", ms=5, mec="k", mfc="r")

            handles["None"] = subplot.plot(*h_ar, **h_kw)[0]
            handles["Parallel β Sheet"] = subplot.plot(*h_ar, **h_kw)[0]
            handles["Antiparallel β Sheet"] = subplot.plot(*h_ar, **h_kw)[0]
            handles["$3_{10}$ Helix"] = subplot.plot(*h_ar, **h_kw)[0]
            handles["α Helix"] = subplot.plot(*h_ar, **h_kw)[0]
            handles["π Helix"] = subplot.plot(*h_ar, **h_kw)[0]
            handles["Turn"] = subplot.plot(*h_ar, **h_kw)[0]
            handles["Bend"] = subplot.plot(*h_ar, **h_kw)[0]
#################################### MAIN #####################################
if __name__ == "__main__":
    SSFigureManager().main()
