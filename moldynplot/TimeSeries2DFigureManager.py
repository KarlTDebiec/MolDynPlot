#!/usr/bin/python
# -*- coding: utf-8 -*-
#   moldynplot.TimeSeries2DFigureManager.py
#
#   Copyright (C) 2015-2017 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Generates 2D time series figures to specifications
"""
################################### MODULES ###################################
from __future__ import (absolute_import, division, print_function,
    unicode_literals)
if __name__ == "__main__":
    __package__ = str("moldynplot")
    import moldynplot
from .myplotspec.FigureManager import FigureManager
from .myplotspec.manage_defaults_presets import manage_defaults_presets
from .myplotspec.manage_kwargs import manage_kwargs
################################### CLASSES ###################################
class TimeSeries2DFigureManager(FigureManager):
    """
    Manages the generation of 2D time series figures.

    **Supported Presets:**

    Per-residue root-mean standard deviation (``perresrmsd``):

    .. image:: _static/p53/perresrmsd.png
        :scale: 35
        :align: center

    Dictionary of secondary structure prediction (``dssp``):

    .. image:: _static/p53/dssp.png
        :scale: 35
        :align: center
    """

    defaults = """
        draw_figure:
          subplot_kw:
            autoscale_on: False
          multi_tick_params:
            bottom: on
            left: on
            right: off
            top: off
          shared_legend: False
          shared_legend_kw:
            spines: False
            handle_kw:
              ls: none
              marker: s
              mec: black
            legend_kw:
              borderaxespad: 0
              frameon: False
              handletextpad: 0
              loc: 9
              numpoints: 1
        draw_subplot:
          grid: True
          grid_kw:
            color: [0.3, 0.3, 0.3]
          hline_kw:
            color: [0.0, 0.0, 0.0]
            zorder: 9
          label_kw:
            horizontalalignment: left
            verticalalignment: top
            zorder: 10
          tick_params:
            bottom: on
            direction: out
            top: off
            left: on
            right: off
          title_kw:
            verticalalignment: bottom
          xlabel: Time
        draw_dataset:
          colorbar_kw:
            ztick_params:
              bottom: off
              left: off
              right: off
              top: off
          draw_heatmap: True
          heatmap_kw:
            edgecolors: none
            rasterized: True
            zorder: 0.1
    """

    available_presets = """
      dssp:
        class: content
        help: Dynamic secondary structure of proteins calculated by cpptraj
        draw_figure:
          shared_legend: True
          shared_legend_kw:
            handles:
              - ["None",
                 {color: !!python/object/apply:moldynplot.dssp_cmap [0,0,7]}]
              - ["Parallel β Sheet",
                 {color: !!python/object/apply:moldynplot.dssp_cmap [1,0,7]}]
              - ["Antiparallel β Sheet",
                 {color: !!python/object/apply:moldynplot.dssp_cmap [2,0,7]}]
              - ["$3_{10}$ Helix",
                 {color: !!python/object/apply:moldynplot.dssp_cmap [3,0,7]}]
              - ["α Helix",
                 {color: !!python/object/apply:moldynplot.dssp_cmap [4,0,7]}]
              - ["π Helix",
                 {color: !!python/object/apply:moldynplot.dssp_cmap [5,0,7]}]
              - ["Turn",
                 {color: !!python/object/apply:moldynplot.dssp_cmap [6,0,7]}]
              - ["Bend",
                 {color: !!python/object/apply:moldynplot.dssp_cmap [7,0,7]}]
            legend_kw:
              title: Secondary Structure
        draw_subplot:
          ylabel: Residue
        draw_dataset:
          dataset_kw:
            cls: moldynplot.dataset.TimeSeriesDataset.TimeSeriesDataset
            downsample_mode: mode
          heatmap_kw:
            cmap: !!python/object/apply:moldynplot.dssp_cmap []
            vmin: 0
            vmax: 7
      perresrmsd:
        class: content
        help: Per-residue RMSD calculated by cpptraj
        draw_figure:
          shared_legend: False
        draw_subplot:
          ylabel: Residue
        draw_dataset:
          dataset_kw:
            cls: moldynplot.dataset.TimeSeriesDataset.TimeSeriesDataset
            downsample_mode: mean
          draw_colorbar: True
          heatmap_kw:
            cmap: afmhot_r
            vmax: 10
            vmin: 0
          colorbar_kw:
            zlabel: Per-Residue Backbone RMSD (Å)
            zticks: [0,1,2,3,4,5,6,7,8,9,10]
      saxs:
        class: content
        help: Small-angle X-ray scattering calculated by saxs_md
        draw_figure:
          shared_legend: False
        draw_subplot:
          ylabel:      "q ($Å^{-1}$)"
          yticks:      [0.00,0.05,0.10,0.15,0.20,0.25,0.30,0.35]
          yticklabels: ["0.00","0.05","0.10","0.15","0.20","0.25","0.30",
          "0.35"]
          grid_kw:
            color: [0.5,0.5,0.5]
        draw_dataset:
          dataset_kw:
            cls: moldynplot.dataset.SAXSTimeSeriesDataset.SAXSTimeSeriesDataset
          heatmap_kw:
            cmap: bone
            vmin: 0
            vmax: 10000000
          draw_colorbar: True
          colorbar_kw:
            zlabel: Intensity
            zticks: [0,5000000,10000000]
      saxs_log:
        class: content
        help: Intensity on log scale
        extends: saxs
        draw_dataset:
          logz: True
          heatmap_kw:
            vmin: -4
            vmax: 0
          colorbar_kw:
            zlabel: $log_{10}$(Intensity)
            zticks: [-4,-3,-2,-1,0]
          contour_kw:
            levels: [-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5]
      manuscript:
        class: target
        inherits: manuscript
        draw_figure:
          bottom: 0.35
          hspace: 0.10
          left: 0.50
          right: 0.10
          shared_legend_kw:
            bottom: 0.00
            handle_kw:
              mew: 0.5
              ms: 5
            left: 0.50
            legend_kw:
              ncol: 4
            sub_height: 0.40
            sub_width: 3.80
          shared_ylabel_kw:
            left: -0.34
          sub_height: 1.00
          sub_width:  3.80
          title_kw:
            top: -0.1
          top:        0.10
          wspace:     0.10
        draw_subplot:
          draw_label: True
          label_kw:
            border_lw: 1
            xabs:  0.020
            yabs: -0.025
          xlabel_kw:
            labelpad: 3
          ylabel_kw:
            labelpad: 6
        draw_dataset:
          colorbar_kw:
            zlabel_kw:
              labelpad: 0
            ztick_params:
              pad: 0
          partner_kw:
            bottom:     0.30
            hspace:     0.42
            left:       1.82
            position: bottom
            sub_height: 0.05
            sub_width:  1.76
      manuscript_stacked_dssp:
        class: target
        extends: manuscript
        help: Series of vertically stacked DSSP plots
        draw_figure:
          sub_width: 5.90
          bottom:    0.65
          hspace:    0.00
          shared_ylabel_kw:
            left: -0.37
          shared_legend_kw:
            sub_width:  5.90
            sub_height: 0.30
            legend_kw:
              loc: 10
              ncol: 8
              title: NULL
        draw_subplot:
          ylabel:
          ylabel_kw:
            labelpad: 4
            rotation: horizontal
            verticalalignment: center
          yticklabels: []
          ylabel_fp: 7r
      manuscript_dssp_legend_right:
        class: target
        extends: manuscript
        help: Stacked DSSP plot; per-residue RMSD plot to be generated
              separately
        draw_figure:
          sub_width:  5.15
          right:      1.35
          bottom:     0.05
          top:        0.25
          title_kw:
            top: -0.1
          shared_legend_kw:
            mode: partner
            partner_kw:
              position:   right
              sub_width:  1.25
              wspace:     0.05
              sub_height: 1.15
            legend_kw:
              title:  "Secondary Structure"
              ncol:   1
        draw_subplot:
          xlabel: null
          tick_params:
            bottom: off
      manuscript_dssp_perresrmsd:
        class: target
        extends: manuscript
        help: Stacked DSSP and per-residue RMSD plots
        draw_figure:
          nrows: 2
          right:      1.35
          sub_width:  5.15
          left:       0.50
          bottom:     0.40
          multiplot: True
          shared_legend: True
          shared_legend_kw:
            mode: "partner"
            partner_kw:
              position: right
              sub_width:  1.24
              wspace:     0.01
              sub_height: 1.15
            handles:
              - ["None",
                 {color: !!python/object/apply:moldynplot.dssp_cmap [0,0,7]}]
              - ["Parallel β Sheet",
                 {color: !!python/object/apply:moldynplot.dssp_cmap [1,0,7]}]
              - ["Antiparallel β Sheet",
                 {color: !!python/object/apply:moldynplot.dssp_cmap [2,0,7]}]
              - ["$3_{10}$ Helix",
                 {color: !!python/object/apply:moldynplot.dssp_cmap [3,0,7]}]
              - ["α Helix",
                 {color: !!python/object/apply:moldynplot.dssp_cmap [4,0,7]}]
              - ["π Helix",
                 {color: !!python/object/apply:moldynplot.dssp_cmap [5,0,7]}]
              - ["Turn",
                 {color: !!python/object/apply:moldynplot.dssp_cmap [6,0,7]}]
              - ["Bend",
                 {color: !!python/object/apply:moldynplot.dssp_cmap [7,0,7]}]
            legend_kw:
              title: Secondary Structure
              loc: 2
              ncol: 1
          subplots:
            0:
              preset: dssp
            1:
              preset: perresrmsd
        draw_subplot:
          datasets:
            0:
              colorbar_kw:
                zlabel: "Backbone RMSD (Å)"
        draw_dataset:
          partner_kw:
            position: right
            left:       null
            sub_width:  1.000
            wspace:     0.150
            bottom:     0.400
            sub_height: 0.063
          colorbar_kw:
            position: top
            zlabel: "Backbone RMSD (Å)"
            zlabel_kw:
              rotation: 0
              labelpad: 5
            zticklabels: [0,"",2,"",4,"",6,"",8,"",10]
    """

    @manage_defaults_presets()
    @manage_kwargs()
    def draw_dataset(self, subplot, label=None, handles=None, logz=False,
            draw_heatmap=False, draw_colorbar=False, draw_contour=False,
            draw_legend=False, draw_label=True, **kwargs):
        import numpy as np
        import six
        from .myplotspec import get_colors, multi_get_copy

        # Process arguments and load data
        verbose = kwargs.get("verbose", 1)
        dataset_kw = multi_get_copy("dataset_kw", kwargs, {})
        if "infile" in kwargs:
            dataset_kw["infile"] = kwargs["infile"]
        dataset = self.load_dataset(verbose=verbose, **dataset_kw)
        timeseries_df = dataset.timeseries_df

        # Draw heatmap
        if draw_heatmap:
            heatmap_kw = multi_get_copy("heatmap_kw", kwargs, {})
            hm_x = timeseries_df.index.values
            if "hm_y" in heatmap_kw:
                hm_y = heatmap_kw.pop("hm_y")
            else:
                try:
                    hm_y = np.array(timeseries_df.columns, np.float)
                except:
                    hm_y = np.array(range(1, timeseries_df.shape[1] + 2))
            hm_z = timeseries_df.values.T
            if logz:
                hm_z = np.log10(hm_z)
            pcolormesh = subplot.pcolor(hm_x, hm_y, hm_z, **heatmap_kw)

            # Draw colorbar
            if draw_colorbar:
                from .myplotspec.axes import set_colorbar
                if not hasattr(subplot, "_mps_partner_subplot"):
                    from .myplotspec.axes import add_partner_subplot
                    add_partner_subplot(subplot, **kwargs)
                set_colorbar(subplot, pcolormesh, **kwargs)

        # Draw contour
        if draw_contour:
            contour_kw = multi_get_copy("contour_kw", kwargs, {})
            get_colors(contour_kw)
            if "levels" not in contour_kw:
                contour_kw["levels"] = range(0,
                    int(np.ceil(np.nanmax(timeseries_df.values))))
            ct_x = timeseries_df.index.values
            try:
                ct_y = np.array(timeseries_df.columns, np.float)
            except:
                ct_y = np.array(range(1, timeseries_df.shape[1] + 2))
            ct_z = timeseries_df.values.T
            if logz:
                ct_z = np.log10(ct_z)
            contour = subplot.contour(ct_x, ct_y, ct_z, **contour_kw)

        # Draw label
        if draw_label and label is not None:
            from .myplotspec.text import set_text

            label_kw = multi_get_copy("label_kw", kwargs, {})
            if (isinstance(label, six.string_types) and isinstance(label_kw,
                dict)):
                set_text(subplot, s=label, **label_kw)
            elif (isinstance(label, list) and isinstance(label_kw,
                    list) and len(label) == len(label_kw)):
                for l, kw in zip(label, label_kw):
                    set_text(subplot, s=l, **kw)
            else:
                raise Exception("bad label arguments")


#################################### MAIN #####################################
if __name__ == "__main__":
    TimeSeries2DFigureManager().main()
