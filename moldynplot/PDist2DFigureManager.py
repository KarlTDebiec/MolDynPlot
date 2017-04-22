#!/usr/bin/python
# -*- coding: utf-8 -*-
#   moldynplot.PDist2DFigureManager.py
#
#   Copyright (C) 2015-2017 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Generates 2D probability distribution figures to specifications
"""
################################### MODULES ###################################
from __future__ import (absolute_import, division, print_function,
    unicode_literals)
if __name__ == "__main__":
    __package__ = str("moldynplot")
    import moldynplot
from IPython import embed
from .myplotspec.FigureManager import FigureManager
from .myplotspec.manage_defaults_presets import manage_defaults_presets
from .myplotspec.manage_kwargs import manage_kwargs
################################### CLASSES ###################################
class PDist2DFigureManager(FigureManager):
    """
    Manages the generation of 2D probability distribution figures.
    """

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
          ylabel: Distance
          tick_params:
            direction: out
            left: on
            right: off
            bottom: on
            top: off
          grid: True
          grid_kw:
            b: True
            color: [0.8,0.8,0.8]
            linestyle: '-'
          label_kw:
            zorder: 10
            horizontalalignment: left
            verticalalignment: top
        draw_dataset:
          draw_heatmap: True
          heatmap_kw:
            edgecolors: none
            rasterized: True
            zorder: 0.1
    """
    available_presets = """
      distance:
        class: content
        help: Plot inter-atomic distance in terms of probability
        draw_figure:
          multi_yticklabels: [0,10,20,30,40,50,60,70]
        draw_subplot:
          yticks: [0,10,20,30,40,50,60,70]
        draw_dataset:
          draw_heatmap: True
          min_cutoff: 0.000001
          heatmap_kw:
            cmap: afmhot
            vmin: 0.000
            vmax: 0.05
          draw_colorbar: True
          colorbar_kw:
            zticks: [0.00,0.01,0.02,0.03,0.04,0.05]
            zticklabels: []
            zlabel: Probability
      distance_free_energy:
        class: content
        help: Plot inter-atomic distance in terms of free energy
        extends: distance
        draw_dataset:
          logz: True
          heatmap_kw:
            vmin: 0
            vmax: 5
          colorbar_kw:
            zticks: [0,1,2,3,4,5]
            zlabel: ΔG (kcal/mol)
      manuscript:
        class: target
        inherits: manuscript
        draw_figure:
          left:       0.50
          sub_width:  2.50
          wspace:     0.10
          right:      0.20
          bottom:     0.70
          sub_height: 1.00
          hspace:     0.10
          top:        0.25
          title_kw:
            top: -0.1
          shared_legend_kw:
            left:       0.50
            sub_width:  2.50
            bottom:     0.00
            sub_height: 0.30
            handle_kw:
              ms: 5
            legend_kw:
              labelspacing: 0.5
              legend_fp: 7r
              ncol: 6
        draw_subplot:
          xlabel_kw:
            labelpad: 3
          ylabel_kw:
            labelpad: 6
          y2ticks: []
          y2label_kw:
            rotation: 270
            verticalalignment: bottom
          draw_label: True
          label_kw:
            border_lw: 1
            xabs:  0.020
            yabs: -0.025
        draw_dataset:
          partner_kw:
            position: right
            wspace:    0.05
            sub_width: 0.05
          colorbar_kw:
            ztick_fp:  6r
            zlabel_fp: 8b
            zlabel_kw:
              labelpad: 2
    """

    @manage_defaults_presets()
    @manage_kwargs()
    def draw_dataset(self, subplot, min_cutoff=None, logz=False,
        draw_heatmap=True, draw_colorbar=True, **kwargs):
        import numpy as np
        from .myplotspec import get_colors, multi_get_copy

        # Process arguments and load data
        verbose = kwargs.get("verbose", 1)
        dataset_kw = multi_get_copy("dataset_kw", kwargs, {})
        if "infile" in kwargs:
            dataset_kw["infile"] = kwargs["infile"]
        dataset = self.load_dataset(verbose=verbose, **dataset_kw)
        if dataset is not None and hasattr(dataset, "pdist_df"):
            pdist_df = dataset.pdist_df
        else:
            pdist_df = dataset.dataframe
        x = np.array( [filter(lambda x: x in "0123456789.", s)
              for s in pdist_df.columns.values], np.int)

        # Configure plot settings
        plot_kw = multi_get_copy("plot_kw", kwargs, {})
        get_colors(plot_kw, kwargs)

        # Draw heatmap
        if draw_heatmap:
            heatmap_kw = multi_get_copy("heatmap_kw", kwargs, {})
            hm_x = np.arange(x.min(), x.max()+2, 1, np.int)
            hm_y = pdist_df.index.values
            hm_z = np.zeros((hm_x.size, hm_y.size))
            for index_within_pdist, residue in enumerate(x):
                pdist_of_residue = pdist_df.values[:, index_within_pdist]
                index_within_hm = np.argmax(hm_x == residue)
                hm_z[index_within_hm] = pdist_of_residue
                print(
                  "Residue {0} is at x[{1}](={2}) and hm_x[{3}](={4})".format(
                  residue, index_within_pdist, x[index_within_pdist],
                  index_within_hm, hm_x[index_within_hm]))
            hm_z = hm_z.T
            if min_cutoff is not None:
                hm_z[hm_z < min_cutoff] = np.inf
            if logz:
                hm_z = -1 * np.log10(hm_z)
            pcolormesh = subplot.pcolor(hm_x - 0.5, hm_y, hm_z, **heatmap_kw)

            # Draw colorbar
            if draw_colorbar:
                from .myplotspec.axes import set_colorbar
                if not hasattr(subplot, "_mps_partner_subplot"):
                    from .myplotspec.axes import add_partner_subplot
                    add_partner_subplot(subplot, **kwargs)
                set_colorbar(subplot, pcolormesh, **kwargs)

#################################### MAIN #####################################
if __name__ == "__main__":
    PDist2DFigureManager().main()
