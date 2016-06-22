#!/usr/bin/python
# -*- coding: utf-8 -*-
#   moldynplot.PDistFigureManager.py
#
#   Copyright (C) 2015-2016 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Generates one or more probability distribution figures to specifications
in a YAML file.
"""
################################### MODULES ###################################
from __future__ import absolute_import,division,print_function,unicode_literals
if __name__ == "__main__":
    __package__ = str("moldynplot")
    import moldynplot
from .myplotspec.FigureManager import FigureManager
################################### CLASSES ###################################
class PDistFigureManager(FigureManager):
    """
    Manages the generation of probability distribution figures.
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
          shared_legend: True
          shared_legend_kw:
            spines: False
            handle_kw:
              ls: none
              marker: s
              mec: black
            legend_kw:
              frameon: False
              loc: 9
              numpoints: 1
              handletextpad: 0
        draw_subplot:
          title_kw:
            verticalalignment: bottom
          ylabel:      "Probability Distribution"
          yticklabels: []
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
          plot_kw:
            zorder: 10
          fill_between_kw:
            color: [0.7, 0.7, 0.7]
            lw: 0
            ylb: 0
            yub: 1
          handle_kw:
            ls: none
            marker: s
            mec: black
    """
    available_presets = """
      rg:
        class: content
        help: Radius of Gyration (Rg)
        draw_figure:
          multi_xticklabels: [0,5,10,15,20,25,30]
        draw_subplot:
          xlabel: $R_g$ (Å)
          xticks: [0,5,10,15,20,25,30]
        draw_dataset:
          ykey: rg
          dataset_kw:
            cls: moldynplot.Dataset.TimeSeriesDataset
            calc_pdist: True
            pdist_kw:
                bandwidth: 0.1
            read_csv_kw:
              delim_whitespace: True
              header: 0
              index_col: 0
              names: [frame, rg, rgmax]
      rmsd:
        class: content
        help: Root Mean Standard Deviation (RMSD)
        draw_figure:
          multi_xticklabels: [0,1,2,3,4,5]
        draw_subplot:
          xlabel: RMSD (Å)
          xticks: [0,1,2,3,4,5]
        draw_dataset:
          ykey: rmsd
          dataset_kw:
            cls: moldynplot.Dataset.TimeSeriesDataset
            calc_pdist: True
            pdist_kw:
                bandwidth: 0.1
                grid: !!python/object/apply:numpy.linspace [0,5,1000]
            read_csv_kw:
              delim_whitespace: True
              header: 0
              index_col: 0
              names: [frame, rmsd]
      r1:
        class: content
        help: Format subplot for R1 relaxation
        draw_subplot:
          xlabel:      "$R_1$"
          xticks:      [0.0,0.5,1.0,1.5,2.0,2.5,3.0]
      r2:
        class: content
        help: Format subplot for R2 relaxation
        draw_subplot:
          xlabel: "$R_2$"
          xticks: [0,2,4,6,8,10,12,14,16,18,20]
      hetnoe:
        class: content
        help: Format subplot for Heteronuclear NOE relaxation
        draw_subplot:
          xlabel: "Heteronuclear NOE"
          xticks: [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
      relaxation_3:
        class: content
        help: Three stacked plots including R1, R2, and HetNOE
        draw_figure:
          nrows: 3
          shared_ylabel: "Prbability Distribution"
          subplots:
            all:
              ylabel: null
            0:
              preset: r1
            1:
              preset: r2
            2:
              preset: hetnoe
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
            wspace:    0.10
            sub_width: 0.80
            title_fp: 8b
            xlabel_kw:
              labelpad: 8.5
            label_fp: 8b
            tick_fp: 6r
            tick_params:
              length: 2
              pad: 3
              width: 1
          plot_kw:
            lw: 1
          handle_kw:
            ms: 6
            mew: 1
    """

    @manage_defaults_presets()
    @manage_kwargs()
    def draw_dataset(self, subplot, label=None, ykey=None, handles=None,
        draw_pdist=True, draw_fill_between=False, verbose=1, debug=0,
        **kwargs):
        """
        """
        from warnings import warn
        from .myplotspec import get_colors, multi_get_copy

        # Load data
        dataset_kw = multi_get_copy("dataset_kw", kwargs, {})
        if "infile" in kwargs:
            dataset_kw["infile"] = kwargs["infile"]
        dataset = self.load_dataset(verbose=verbose, debug=debug, **dataset_kw)
        if dataset is not None and hasattr(dataset, "pdist"):
            pdist = dataset.pdist
        else:
            pdist = None

        # Configure plot settings
        plot_kw = multi_get_copy("plot_kw", kwargs, {})
        get_colors(plot_kw, kwargs)

        # Plot fill_between
        if draw_fill_between:
            fill_between_kw = multi_get_copy("fill_between_kw", kwargs, {})
            get_colors(fill_between_kw, plot_kw)

            if "x" in fill_between_kw:
                fb_x = fill_between_kw.pop("x")
            if "ylb" in fill_between_kw:
                fb_ylb = fill_between_kw.pop("ylb")
            if "yub" in fill_between_kw:
                fb_yub = fill_between_kw.pop("yub")
            subplot.fill_between(fb_x, fb_ylb, fb_yub, **fill_between_kw)

        # Plot pdist
        if draw_pdist:
            if not hasattr(dataset, "pdist"):
                warn("'draw_pdist' is enabled but dataset does not have the "
                     "necessary attribute 'pdist', skipping.")
            else:
                pdist = pdist[ykey]
                pdist_kw = plot_kw.copy()
                pdist_kw.update(kwargs.get("pdist_kw", {}))

                pd_x = pdist.index.values
                pd_y = pdist.values

                subplot.plot(pd_x, pd_y, **pdist_kw)
                pdist_rescale = True
                if pdist_rescale:
                    pdist_max = pd_y.max()
                    y_max = subplot.get_ybound()[1]
                    if (pdist_max > y_max / 1.25
                    or not hasattr(subplot, "_mps_rescaled")):
                        subplot.set_ybound(0, pdist_max*1.25)
                        yticks = [0, pdist_max*0.25, pdist_max*0.50,
                          pdist_max*0.75, pdist_max, pdist_max*1.25]
                        subplot.set_yticks(yticks)
                        subplot._mps_rescaled = True

            if draw_fill_between:
                subplot.fill_between(fb_x, fb_ylb, fb_yub, **fill_between_kw)

#################################### MAIN #####################################
if __name__ == "__main__":
    PDistFigureManager().main()
