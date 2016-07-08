#!/usr/bin/python
# -*- coding: utf-8 -*-
#   moldynplot.TimeSeries2DFigureManager.py
#
#   Copyright (C) 2015-2016 Karl T Debiec
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
from .myplotspec.manage_defaults_presets import manage_defaults_presets
from .myplotspec.manage_kwargs import manage_kwargs
################################### CLASSES ###################################
class TimeSeries2DFigureManager(FigureManager):
    """
    Manages the generation of 2D time series figures.
    """

    defaults = """
        draw_figure:
          subplot_kw:
            autoscale_on: False
          multi_tick_params:
            left:   on
            right:  off
            bottom: on
            top:    off
          shared_legend_kw:
            spines: False
            handle_kw:
              ls: none
              marker: s
              mec: black
            legend_kw:
              frameon: False
              handletextpad: 0
              loc: 9
              numpoints: 1
        draw_subplot:
          title_kw:
            verticalalignment: bottom
          xlabel: Time
          tick_params:
            direction: out
            bottom: on
            top:    off
            right:  off
            left:   on
          grid: True
          grid_kw:
            b: True
            color: [0.2,0.2,0.2]
            linestyle: '-'
            lw: 0.5
          hline_kw:
            color: [0.2,0.2,0.2]
            linestyle: '-'
            lw: 0.5
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
          contour_kw:
            colors: '0.25'
            levels: [1, 2, 3, 4, 5]
            linestyles: solid
            zorder: 0.2
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
        draw_subplot:
          ylabel: Residue
        draw_dataset:
          dataset_kw:
            cls: moldynplot.Dataset.TimeSeriesDataset
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
            cls: moldynplot.Dataset.TimeSeriesDataset
            downsample_mode: mean
          heatmap_kw:
            cmap: afmhot_r
            vmin: 0
            vmax: 10
          draw_colorbar: True
          colorbar_kw:
            zticks: [0,1,2,3,4,5,6,7,8,9,10]
            zlabel: Per-Residue Backbone RMSD (Å)
      saxs:
        class: content
        help: Small-angle X-ray scattering calculated by saxs_md
        draw_figure:
          shared_legend: False
        draw_subplot:
          ylabel:      "q ($Å^{-1}$)"
          yticks:      [0.00,0.05,0.10,0.15,0.20,0.25,0.30,0.35]
          yticklabels: ["0.00","0.05","0.10","0.15","0.20","0.25","0.30","0.35"]
          grid_kw:
            color: [0.5,0.5,0.5]
        draw_dataset:
          dataset_kw:
            cls: moldynplot.Dataset.SAXSTimeSeriesDataset
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
          left:       0.50
          sub_width:  4.40
          wspace:     0.05
          right:      0.20
          bottom:     0.80
          sub_height: 1.00
          hspace:     0.05
          top:        0.25
          title_kw:
            top: -0.1
          shared_legend_kw:
            left:       0.50
            sub_width:  4.40
            bottom:     0.00
            sub_height: 0.40
            handle_kw:
              ms: 5
            legend_kw:
              ncol: 4
        draw_subplot:
          xlabel_kw:
            labelpad: 3
          ylabel_kw:
            labelpad: 6
          draw_label: True
          label_kw:
            border_lw: 1
            xabs:  0.020
            yabs: -0.025
        draw_dataset:
          partner_kw:
            position: bottom
            hspace:     0.42
            sub_height: 0.05
            left:       1.82
            sub_width:  1.76
            bottom:     0.30
          colorbar_kw:
            ztick_fp:  6r
            zlabel_fp: 8b
            zlabel_kw:
              labelpad: 2
      manuscript_stacked_dssp:
        class: target
        extends: manuscript
        help: Series of vertically stacked DSSP plots
        draw_figure:
          sub_width: 6.30
          bottom:    0.60
          hspace:    0.00
          shared_ylabel_kw:
            left: -0.37
          shared_legend_kw:
            sub_width:  6.30
            sub_height: 0.25
            legend_kw:
              ncol: 8
        draw_subplot:
          ylabel:
          ylabel_kw:
            labelpad: 4
            rotation: horizontal
            verticalalignment: center
          yticklabels: []
          ylabel_fp: 7r
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
      notebook:
        class: target
        inherits: notebook
        draw_figure:
          left:       0.60
          sub_width:  5.00
          wspace:     0.10
          right:      0.25
          bottom:     1.10
          sub_height: 2.00
          hspace:     0.10
          top:        0.35
          shared_legend_kw:
            left:       0.50
            sub_width:  5.00
            right:      0.25
            bottom:     0.10
            sub_height: 0.50
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
    def draw_dataset(self, subplot, label=None,
        handles=None, logz=False,
        draw_heatmap=False, draw_colorbar=False, draw_contour=False,
        draw_legend=False, draw_label=True, 
        verbose=1, debug=0, **kwargs):
        import numpy as np
        import six
        from .myplotspec import get_colors, multi_get_copy

        # Load data
        dataset_kw = multi_get_copy("dataset_kw", kwargs, {})
        if "infile" in kwargs:
            dataset_kw["infile"] = kwargs["infile"]
        dataset = self.load_dataset(verbose=verbose, debug=debug, **dataset_kw)
        timeseries_df = dataset.timeseries_df

        # Draw heatmap
        if draw_heatmap:
            heatmap_kw = multi_get_copy("heatmap_kw", kwargs, {})
            hm_x = timeseries_df.index.values
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

                    add_partner_subplot(subplot, verbose=verbose,
                      debug=debug, **kwargs)

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
            if (isinstance(label, six.string_types)
            and isinstance(label_kw, dict)):
                set_text(subplot, s=label, **label_kw)
            elif (isinstance(label, list) and isinstance(label_kw, list)
            and len(label) == len(label_kw)):
                for l, kw in zip(label, label_kw):
                    set_text(subplot, s=l, **kw)
            else:
                raise Exception("bad label arguments")

#################################### MAIN #####################################
if __name__ == "__main__":
    TimeSeries2DFigureManager().main()
