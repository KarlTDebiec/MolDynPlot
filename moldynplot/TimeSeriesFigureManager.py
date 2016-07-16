#!/usr/bin/python
# -*- coding: utf-8 -*-
#   moldynplot.TimeSeriesFigureManager.py
#
#   Copyright (C) 2015-2016 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Generates one or more time series figures to specifications in a YAML
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
class TimeSeriesFigureManager(FigureManager):
    """
    Manages the generation of time series figures.

    .. image:: _static/p53/rmsd.png
    """

    defaults = """
        draw_figure:
          subplot_kw:
            autoscale_on: False
            axisbg: none
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
          xlabel: Time
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
          partner_kw:
            position: right
            y2label_kw:
              rotation: 270
              verticalalignment: bottom
            tick_params:
              direction: out
              bottom: on
              top: off
              right: off
              left: off
            grid: True
            grid_kw:
              b: True
              color: [0.8,0.8,0.8]
              linestyle: '-'
            xticks:
            tick_params:
              direction: out
              bottom: on
              top: off
              right: off
              left: off
          plot_kw:
            zorder: 10
          fill_between_kw:
            color: [0.7, 0.7, 0.7]
            lw: 0
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
      natcon:
        class: content
        help: "% Native contacts vs. time"
        draw_subplot:
          ylabel: "% Native Contacts"
          yticks:      [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
          yticklabels: [0,10,20,30,40,50,60,70,80,90,100]
        draw_dataset:
          column: percent_native_contacts
          dataset_kw:
            cls: moldynplot.Dataset.NatConTimeSeriesDataset
            downsample_mode: mean
          partner_kw:
            yticks: [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
          plot_kw:
            drawstyle: steps
      rmsd:
        class: content
        help: Root Mean Standard Deviation (RMSD)
        draw_subplot:
          ylabel: RMSD (Å)
          yticks: [0,2,4,6,8,10]
        draw_dataset:
          column: rmsd
          partner_kw:
            yticks: [0,2,4,6,8,10]
          dataset_kw:
            cls: moldynplot.Dataset.TimeSeriesDataset
            calc_pdist: True
            pdist_kw:
              bandwidth: 0.1
              grid: !!python/object/apply:numpy.linspace [0,10,1000]
            read_csv_kw:
              header: 0
              names: [frame, rmsd]
      rg:
        class: content
        help: Radius of gyration (Rg)
        draw_figure:
          multi_yticklabels: [0,5,10,15,20,25,30]
        draw_subplot:
          ylabel: $R_g$ (Å)
          yticks: [0,5,10,15,20,25,30]
        draw_dataset:
          column: rg
          partner_kw:
            yticks: [0,5,10,15,20,25,30]
          dataset_kw:
            cls: moldynplot.Dataset.TimeSeriesDataset
            calc_pdist: True
            pdist_kw:
              bandwidth: 0.1
              grid: !!python/object/apply:numpy.linspace [0,30,1000]
            read_csv_kw:
              delim_whitespace: True
              header: 0
              names: [frame, rg, rgmax]
      rotdif:
        class: content
        help: Rotational correlatio time (τc)
        draw_subplot:
          ylabel: $τ_c$ (ns)
          yticks: [0,5,10,15]
        draw_dataset:
          column: rotdif
          partner_kw:
            yticks: [0,5,10,15]
          dataset_kw:
            cls: moldynplot.Dataset.TimeSeriesDataset
            calc_pdist: True
            pdist_kw:
              bandwidth:
                rotdif: 0.2
      presentation:
        class: target
        inherits: presentation
        draw_figure:
          left:       1.20
          sub_width:  7.00
          bottom:     3.10
          sub_height: 3.00
          shared_legend: True
          shared_legend_kw:
            left:       1.20
            sub_width:  7.00
            bottom:     1.90
            sub_height: 0.50
            legend_kw:
              labelspacing: 0.0
              ncol: 5
        draw_dataset:
          partner_kw:
            sub_width: 1.2
            title_fp: 18r
            xlabel_kw:
              labelpad: 20
            label_fp: 18r
            tick_fp: 14r
            xticks:
            lw: 2
            tick_params:
              length: 3
              pad: 6
              width: 2
          handle_kw:
            ms: 12
            mew: 2
      manuscript:
        class: target
        inherits: manuscript
        draw_figure:
          left:       0.50
          sub_width:  4.40
          wspace:     0.10
          right:      0.20
          bottom:     0.70
          sub_height: 1.80
          hspace:     0.10
          top:        0.25
          title_kw:
            top: -0.1
          shared_legend_kw:
            left:       0.50
            sub_width:  4.40
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
          mean_kw:
            ms: 2
      notebook:
        class: target
        inherits: notebook
        draw_figure:
          left:       0.60
          sub_width:  4.40
          right:      0.20
          bottom:     1.00
          sub_height: 1.80
          top:        0.30
          shared_legend: True
          shared_legend_kw:
            left:       0.60
            sub_width:  4.40
            right:      0.20
            bottom:     0.00
            sub_height: 0.50
            legend_kw:
              labelspacing: 0.5
              legend_fp: 8r
              ncol: 2
        draw_dataset:
          plot_kw:
            lw: 1.0
          partner_kw:
            sub_width: 0.8
            title_fp: 10b
            xlabel_kw:
              labelpad:  12.5
            label_fp: 10b
            tick_fp: 8r
            xticks:
            tick_params:
              length: 2
              pad: 6
              width: 1
      pdist:
        class: appearance
        help: Draw probability distribution on right side of plot
        draw_dataset:
          draw_pdist: True
          partner_kw:
            xlabel:      Distribution
            xticks:      [0,0.000001]
            xticklabels: []
            yticklabels: []
    """

    @manage_defaults_presets()
    @manage_kwargs()
    def draw_dataset(self, subplot, label=None, column=None, handles=None,
        draw_pdist=False, draw_fill_between=False, draw_mean=False,
        draw_plot=True, **kwargs):
        """
        Draws a dataset on a subplot.

        Loaded dataset should have attribute `timeseries_df`

        Arguments:
          subplot (Axes): :class:`Axes<matplotlib.axes.Axes>` on
            which to draw
          dataset_kw (dict): Keyword arguments passed to
            :meth:`load_dataset
            <myplotspec.FigureManager.FigureManager.load_dataset>`
          plot_kw (dict): Keyword arguments passed to methods of
            :class:`Axes<matplotlib.axes.Axes>`
          draw_plot (bool): Draw plot
          draw_pdist (bool): Draw probability distribution
          draw_fill_between (bool): Fill between specified region for this
            dataset
          draw_mean (bool): Draw point at mean value value of
            probability distribution
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
        if dataset is not None and hasattr(dataset, "timeseries_df"):
            timeseries = dataset.timeseries_df
        else:
            timeseries = None

        # Configure plot settings
        plot_kw = multi_get_copy("plot_kw", kwargs, {})
        get_colors(plot_kw, kwargs)

        # Draw fill_between
        if draw_fill_between:
            fill_between_kw = multi_get_copy("fill_between_kw", kwargs, {})
            get_colors(fill_between_kw, plot_kw)

            if "x" in fill_between_kw:
                fb_x = fill_between_kw.pop("x")
            else:
                fb_x = timeseries.index.values
            if "ylb" in fill_between_kw:
                fb_ylb = fill_between_kw.pop("ylb")
            elif "fill_between_lb_key" in fill_between_kw:
                fill_between_lb_key = fill_between_kw.pop(
                  "fill_between_lb_key")
                fb_ylb = timeseries[fill_between_lb_key]
            elif "fill_between_lb_key" in kwargs:
                fill_between_lb_key = kwargs.get( "fill_between_lb_key")
                fb_ylb = timeseries[fill_between_lb_key]
            else:
                warn("inappropriate fill_between settings")
            if "yub" in fill_between_kw:
                fb_yub = fill_between_kw.pop("yub")
            elif "fill_between_ub_key" in fill_between_kw:
                fill_between_ub_key = fill_between_kw.pop(
                  "fill_between_ub_key")
                fb_yub = timeseries[fill_between_ub_key]
            elif "fill_between_ub_key" in kwargs:
                fill_between_ub_key = kwargs.get( "fill_between_ub_key")
                fb_yub = timeseries[fill_between_ub_key]
            else:
                warn("inappropriate fill_between settings")
            subplot.fill_between(fb_x, fb_ylb, fb_yub, **fill_between_kw)

        # Draw plot
        if draw_plot:
            if verbose >= 2:
                print("mean  {0}: {1:6.3f}".format(column,
                  timeseries[column].mean()))
                print("stdev {0}: {1:6.3f}".format(column,
                  timeseries[column].std()))
            plot = subplot.plot(timeseries.index.values, timeseries[column],
                     **plot_kw)[0]
            handle_kw = multi_get_copy("handle_kw", kwargs, {})
            handle_kw["mfc"] = plot.get_color()
            handle = subplot.plot([-10, -10], [-10, -10], **handle_kw)[0]
            if handles is not None and label is not None:
                handles[label] = handle

        # Draw pdist
        if draw_pdist:

            if not hasattr(dataset, "pdist_df"):
                warn("'draw_pdist' is enabled but dataset does not have the "
                     "necessary attribute 'pdist_df', skipping.")
            else:

                # Add subplot if not already present
                if not hasattr(subplot, "_mps_partner_subplot"):
                    from .myplotspec.axes import add_partner_subplot

                    add_partner_subplot(subplot, **kwargs)

                pdist = dataset.pdist_df[column]
                pdist_kw = plot_kw.copy()
                pdist_kw.update(kwargs.get("pdist_kw", {}))

                pd_x = pdist.index.values
                pd_y = np.squeeze(pdist.values)

                subplot._mps_partner_subplot.plot(pd_y, pd_x, **pdist_kw)
                pdist_rescale = True
                if pdist_rescale:
                    pdist_max = pd_y.max()
                    x_max = subplot._mps_partner_subplot.get_xbound()[1]
                    if (pdist_max > x_max / 1.25
                    or not hasattr(subplot, "_mps_rescaled")):
                        subplot._mps_partner_subplot.set_xbound(0,
                        pdist_max*1.25)
                        xticks = [0, pdist_max*0.25, pdist_max*0.50,
                          pdist_max*0.75, pdist_max, pdist_max*1.25]
                        subplot._mps_partner_subplot.set_xticks(xticks)
                        subplot._mps_rescaled = True
                    if draw_mean:
                        mean_kw = plot_kw.copy()
                        mean_kw.update(kwargs.get("mean_kw", {}))
                        mean = np.sum(np.array(pd_x, np.float64)
                                     *np.array(pd_y, np.float64))
                        subplot._mps_partner_subplot.plot(
                          pd_y[np.abs(pd_x - mean).argmin()], mean, **mean_kw)

            if draw_fill_between:
                subplot._mps_partner_subplot.fill_between(fb_x, fb_ylb,
                  fb_yub, **fill_between_kw)

#################################### MAIN #####################################
if __name__ == "__main__":
    TimeSeriesFigureManager().main()
