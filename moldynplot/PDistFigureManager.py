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
from .myplotspec.manage_defaults_presets import manage_defaults_presets
from .myplotspec.manage_kwargs import manage_kwargs
################################### CLASSES ###################################
class PDistFigureManager(FigureManager):
    """
    Manages the generation of probability distribution figures.
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
              borderaxespad: 0
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
            zorder: 1
          handle_kw:
            ls: none
            marker: s
            mec: black
          mean_kw:
            ls: none
            marker: o
            mec: black
            zorder: 11
    """
    available_presets = """
      radgyr:
        class: content
        help: Radius of Gyration (Rg)
        draw_figure:
          multi_xticklabels: [0,5,10,15,20,25,30]
        draw_subplot:
          xlabel: $R_g$ (Å)
          xticks: [0,5,10,15,20,25,30]
        draw_dataset:
          column: rg
          dataset_kw:
            cls: moldynplot.dataset.TimeSeriesDataset.TimeSeriesDataset
            calc_pdist: True
            pdist_kw:
              bandwidth: 0.1
              grid: !!python/object/apply:numpy.linspace [0,30,1000]
            read_csv_kw:
              delim_whitespace: True
              header: 0
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
          column: rmsd
          dataset_kw:
            cls: moldynplot.dataset.TimeSeriesDataset.TimeSeriesDataset
            calc_pdist: True
            pdist_kw:
                bandwidth: 0.1
                grid: !!python/object/apply:numpy.linspace [0,5,1000]
            read_csv_kw:
              delim_whitespace: True
              header: 0
              names: [frame, rmsd]
      r1:
        class: content
        help: Format subplot for R1 relaxation
        draw_subplot:
          xlabel:      "$R_1$"
          xticks:      [0.0,0.5,1.0,1.5,2.0,2.5,3.0]
        draw_dataset:
          dataset_kw:
            pdist_kw:
              bandwidth: 0.02
          column: r1
      r2:
        class: content
        help: Format subplot for R2 relaxation
        draw_subplot:
          xlabel: "$R_2$"
          xticks: [0,2,4,6,8,10,12,14,16,18,20]
        draw_dataset:
          dataset_kw:
            pdist_kw:
              bandwidth: 0.3
          column: r2
      r2/r1:
        class: content
        help: Format subplot for R2/R1 relaxation
        draw_subplot:
          xlabel: "$R_2$/$R_1$"
          xticks: [3,4,5,6,7,8,9,10,11]
        draw_dataset:
          dataset_kw:
            pdist_kw:
              bandwidth:
                r2/r1: 0.1
          column: r2/r1
      hetnoe:
        class: content
        help: Format subplot for Heteronuclear NOE relaxation
        draw_subplot:
          xlabel: "Heteronuclear NOE"
          xticks: [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
        draw_dataset:
          column: noe
          dataset_kw:
            pdist_kw:
              bandwidth: 0.03
      rotdif:
        class: content
        help: Format subplot for rotational diffusion
        draw_subplot:
          xlabel: "$τ_c$ (ns)"
          xticks: [5,6,7,8,9,10,11,12,13,14]
        draw_dataset:
          column: rotdif
          dataset_kw:
            pdist_kw:
              bandwidth: 0.2
      relaxation_3:
        class: content
        help: Three stacked plots including R1, R2, and HetNOE
        draw_figure:
          nrows: 3
          shared_ylabel: "Probability Distribution"
          subplots:
            0:
              preset: r1
              ylabel: null
            1:
              preset: r2
              ylabel: null
            2:
              preset: hetnoe
              ylabel: null
      relaxation_4:
        class: content
        help: Four stacked plots including R1, R2, R2/R1, and HetNOE
        draw_figure:
          nrows: 4
          shared_ylabel: "Probability Distribution"
          subplots:
            0:
              preset: r1
              ylabel: null
            1:
              preset: r2
              ylabel: null
            2:
              preset: r2/r1
              ylabel: null
            3:
              preset: hetnoe
              ylabel: null
      rotdif_2:
        class: content
        help: Two stacked plots including R2/R1 rotdif
        draw_figure:
          nrows: 2
          shared_ylabel: "Probability Distribution"
          subplots:
            0:
              preset: r2/r1
              ylabel: null
            1:
              preset: rotdif
              ylabel: null
      rotdif_4:
        class: content
        help: Two stacked plots including R2/R1 rotdif
        draw_figure:
          nrows: 2
          ncols: 2
          shared_ylabel: "Probability Distribution"
          subplots:
            0:
              preset: r2/r1
              ylabel: null
            1:
              preset: r2/r1
              ylabel: null
            2:
              preset: rotdif
              ylabel: null
            3:
              preset: rotdif
              ylabel: null
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
          mean_kw:
            ms: 2
          handle_kw:
            ms: 6
            mew: 1
    """

    @manage_defaults_presets()
    @manage_kwargs()
    def draw_dataset(self, subplot, column=None,
        draw_pdist=True, draw_fill_between=False, draw_mean=False, **kwargs):
        """
        Loads a dataset and draws it on a subplot.

        Loaded dataset should have attribute `pdist_df`.

        Arguments:
          subplot (Axes): :class:`Axes<matplotlib.axes.Axes>` on
            which to draw
          dataset_kw (dict): Keyword arguments passed to
            :meth:`load_dataset
            <myplotspec.FigureManager.FigureManager.load_dataset>`
          plot_kw (dict): Keyword arguments passed to methods of
            :class:`Axes<matplotlib.axes.Axes>`
          column (str): Column within `pdist_df` to use
          draw_fill_between (bool): Fill between specified region
          fill_between_kw (dict): Keyword arguments used to configure
            call to
            :meth:`fill_between<matplotlib.axes.Axes.fill_between>`
          fill_between_kw[x] (list, ndarray): x values passed to
            :meth:`fill_between<matplotlib.axes.Axes.fill_between>`
          fill_between_kw[ylb] (list, ndarray): y lower bound values
            passed to
            :meth:`fill_between<matplotlib.axes.Axes.fill_between>`
          fill_between_kw[yub] (list, ndarray): y upper bound values
            passed to
            :meth:`fill_between<matplotlib.axes.Axes.fill_between>`
          draw_pdist (bool): Draw probability distribution
          pdist_kw (dict): Keyword arguments using to configure call to
            :meth:`plot<matplotlib.axes.Axes.plot>`
          draw_mean (bool): Draw point at mean value
          mean_kw (dict): Keyword arguments used to configure call to
            :meth:`plot<matplotlib.axes.Axes.plot>`
          verbose (int): Level of verbose output
          kwargs (dict): Additional keyword arguments
        """
        from warnings import warn
        import numpy as np
        from .myplotspec import get_colors, multi_get_copy

        # Process arguments
        verbose = kwargs.get("verbose", 1)
        dataset_kw = multi_get_copy("dataset_kw", kwargs, {})
        if "infile" in kwargs:
            dataset_kw["infile"] = kwargs["infile"]
        dataset = self.load_dataset(verbose=verbose, **dataset_kw)
        if dataset is not None and hasattr(dataset, "pdist_df"):
            pdist_df = dataset.pdist_df
        else:
            pdist_df = None


        # Configure plot settings
        plot_kw = multi_get_copy("plot_kw", kwargs, {})
        get_colors(plot_kw, kwargs)

        # Draw fill_between
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

        # Draw pdist
        if draw_pdist:
            if not hasattr(dataset, "pdist_df"):
                warn("'draw_pdist' is enabled but dataset does not have the "
                     "necessary attribute 'pdist_df', skipping.")
            else:
                pdist = pdist_df[column]
                pdist_kw = plot_kw.copy()
                pdist_kw.update(kwargs.get("pdist_kw", {}))

                pd_x = pdist.index.values
                pd_y = np.squeeze(pdist.values)

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
                if draw_mean:
                    mean_kw = plot_kw.copy()
                    mean_kw.update(kwargs.get("mean_kw", {}))
                    mean = np.sum(np.array(pd_x, np.float64)
                                 *np.array(pd_y, np.float64))
                    subplot.plot(mean, pd_y[np.abs(pd_x - mean).argmin()],
                      **mean_kw)

#################################### MAIN #####################################
if __name__ == "__main__":
    PDistFigureManager().main()
