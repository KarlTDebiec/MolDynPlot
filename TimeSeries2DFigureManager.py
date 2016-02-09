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
          title_kw:
            verticalalignment: bottom
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
        draw_figure:
          shared_legend: True
        draw_subplot:
          ylabel: Residue
        draw_dataset:
          downsample_mode: mode
          heatmap_kw:
            cmap: !!python/object/apply:moldynplot.dssp_color_palette []
            vmin: 0
            vmax: 7
          legend: True
          labels:
            0: None
            1: Parallel β Sheet
            2: Antiparallel β Sheet
            3: $3_{10}$ Helix
            4: α Helix
            5: π Helix
            6: Turn
            7: Bend
      perres_rmsd:
        class: content
        help: Per-residue RMSD calculated by cpptraj
        draw_dataset:
          downsample_mode: mean
          heatmap_kw:
            cmap: afmhot_r
            vmin: 0
            vmax: 10
          colorbar: True
          colorbar_kw:
            zticks: [0,1,2,3,4,5,6,7,8,9,10]
            zlabel: Per-Residue Backbone RMSD (Å)
      notebook:
        class: target
        inherits: notebook
        draw_figure:
          left:       0.60
          sub_width:  5.00
          right:      0.25
          bottom:     1.10
          sub_height: 2.00
          hspace:     0.10
          top:        0.35
          title_kw:
            top: -0.1
          shared_legend_kw:
            left:       0.50
            sub_width:  5.00
            right:      0.25
            bottom:     0.10
            sub_height: 0.50
            top:        0.60
            legend_kw:
              ncol: 4
        draw_subplot:
          ylabel_kw:
            labelpad: 12
        draw_dataset:
          partner_kw:
            position: bottom
            left:       2.10
            sub_width:  2.00
            right:      1.75
            bottom:     0.40
            sub_height: 0.10
            top:        0.50
          colorbar_kw:
            ztick_fp: 8r
            zlabel_fp: 10b
            zlabel_kw:
              labelpad: 3
    """

    @manage_defaults_presets()
    @manage_kwargs()
    def draw_dataset(self, subplot, downsample=None, label=None, heatmap=True,
        colorbar=False, legend=False, handles=None, labels=None, verbose=1,
        debug=0, **kwargs):
        """
        """
        from os.path import expandvars
        import h5py
        import numpy as np
        from .myplotspec import multi_get_copy

        # Load data
        with h5py.File(expandvars(kwargs.get("infile"))) as h5_file:
            key = kwargs.get("key", h5_file.keys()[0])
            dataset = np.array(h5_file[key])
        if "usecols" in kwargs:
            dataset = dataset[:,kwargs.get("usecols")]

        # Scale:
        dt = kwargs.get("dt", 0.001)
        time = np.arange(dataset.shape[0]) * dt

        # Downsample
        if downsample is not None:
            from scipy.stats.mstats import mode
            full_size = dataset.shape[0] - (dataset.shape[0] % downsample)
            dataset = np.reshape(dataset[:full_size],
              (int(full_size / downsample), downsample, dataset.shape[1]))
            downsample_mode = kwargs.get("downsample_mode", "mean")
            if downsample_mode == "mean":
                dataset = np.mean(dataset, axis=1)
            elif downsample_mode == "mode":
                dataset = np.squeeze(mode(dataset, axis=1)[0])
            time = np.arange(dataset.shape[0]) * (dt * downsample)

        if heatmap:
            heatmap_kw = multi_get_copy("heatmap_kw", kwargs, {})
            y = kwargs.get("y", np.arange(1, dataset.shape[1]+2))
            pcolormesh = subplot.pcolor(time, y, dataset.T, **heatmap_kw)

            if colorbar:
                from .myplotspec.axes import set_colorbar
                if not hasattr(subplot, "_mps_partner_subplot"):
                    from .myplotspec.axes import add_partner_subplot
                    add_partner_subplot(subplot, verbose=verbose,
                      debug=debug, **kwargs)
                set_colorbar(subplot, pcolormesh, **kwargs)

            if legend:
                if handles is not None and labels is not None:
                    cmap = pcolormesh.cmap
                    h_kw=dict(ls="none", marker="s", ms=5, mec="k")

                    for i in sorted(labels.keys()):
                        label = labels[i]
                        handles[label] = subplot.plot((-1),(-1),
                          mfc=cmap(i / (len(labels) - 1)), **h_kw)[0]

#################################### MAIN #####################################
if __name__ == "__main__":
    TimeSeries2DFigureManager().main()
