#!/usr/bin/python
# -*- coding: utf-8 -*-
#   moldynplot.TimeSeriesFigureManager.py
#
#   Copyright (C) 2015 Karl T Debiec
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
################################### CLASSES ###################################
class TimeSeriesFigureManager(FigureManager):
    """
    Manages the generation of time series figures
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
            legend_kw:
              frameon: False
              loc: 9
              numpoints: 1
              handletextpad: 0.0
        draw_subplot:
          title_kw:
            verticalalignment: bottom
          xlabel: Time
          tick_params:
            left: on
            right: off
            bottom: on
            top: off
            direction: out
          grid: True
          grid_kw:
            b: True
            color: [0.8,0.8,0.8]
            linestyle: '-'
        draw_dataset:
          partner_kw:
            position: right
            tick_params:
              bottom: on
              top: off
              right: off
              left: off
              direction: out
            grid: True
            grid_kw:
              b: True
              color: [0.8,0.8,0.8]
              linestyle: '-'
          plot_kw:
            zorder: 10
          handle_kw:
            marker: s
            ls: none
            mec: black
    """

    available_presets = """
      rmsd:
        class: content
        help: Root Mean Standard Deviation (RMSD) vs. time
        draw_subplot:
          ylabel: RMSD (Å)
        draw_dataset:
          dataset_kw:
            read_csv_kw:
              delim_whitespace: True
              index_col: False
              names: [time, rmsd]
              header: 0
          ykey: rmsd
      rg:
        class: content
        help: Radius of Gyration (Rg) vs. time
        draw_subplot:
          ylabel: $R_g$ (Å)
        draw_dataset:
          dataset_kw:
            read_csv_kw:
              delim_whitespace: True
              index_col: False
              names: [time, rg, rgmax]
              header: 0
          ykey: rg
      presentation:
        class: target
        inherits: presentation
        draw_figure:
          left:       1.20
          sub_width:  7.00
          bottom:     2.90
          sub_height: 3.00
          shared_legend: True
          shared_legend_kw:
            left:       1.20
            sub_width:  7.00
            bottom:     1.70
            sub_height: 0.50
            legend_kw:
              labelspacing: 0.5
              ncol: 2
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
      manuscript:
        class: target
        inherits: manuscript
        draw_figure:
          left:       0.50
          sub_width:  4.40
          right:      0.20
          bottom:     0.70
          sub_height: 1.80
          top:        0.20
          shared_legend: True
          shared_legend_kw:
            left:       0.50
            sub_width:  4.40
            sub_height: 0.30
            bottom:     0.00
            legend_kw:
              labelspacing: 0.5
              legend_fp: 7r
              ncol: 5
        draw_dataset:
          partner_kw:
            sub_width: 0.8
            title_fp: 8b
            xlabel_kw:
              labelpad: 12.5
            label_fp: 8b
            tick_fp: 6r
            xticks:
            tick_params:
              direction: out
              length: 2
              pad: 3
              width: 1
          plot_kw:
            lw: 1.0
          handle_kw:
            ms: 6
            mew: 1
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
          pdist: True
          partner_kw:
            xlabel:      Distribution
            xticks:      [0,0.000001]
            xticklabels: []
            yticks:      [0,1,2,3,4,5,6]
            yticklabels: []
    """

    @manage_defaults_presets()
    @manage_kwargs()
    def draw_dataset(self, subplot, dt=None, downsample=None, label=None,
        xoffset=0, handles=None, ykey=None, pdist=False, fill_between=False,
        verbose=1, debug=0, **kwargs):
        from .myplotspec import get_color, multi_get_copy
        from .myplotspec.Dataset import Dataset
        import pandas as pd
        import numpy as np
        from os.path import expandvars

        # Load data
        dataset_kw = multi_get_copy("dataset_kw", kwargs, {})
        if "infile" in kwargs:
            dataset_kw["infile"] = kwargs["infile"]
        dataframe= self.load_dataset(Dataset, verbose=verbose, debug=debug,
          **dataset_kw).data

        # Scale:
        if dt is not None:
            dataframe["time"] *= dt

        # Configure plot settings
        plot_kw = multi_get_copy("plot_kw", kwargs, {})
        if "color" in plot_kw:
            plot_kw["color"] = get_color(plot_kw["color"])
        elif "color" in kwargs:
            plot_kw["color"] = get_color(kwargs.pop("color"))
        fill_between_kw = multi_get_copy("fill_between_kw", kwargs, {})
        if "color" in fill_between_kw:
            fill_between_kw["color"] = get_color(fill_between_kw["color"])
        elif "color" in plot_kw:
            fill_between_kw["color"] = get_color(plot_kw["color"])

        # Downsample
        if downsample is not None:
            full_size = dataframe.shape[0]
            reduced_size = int(full_size / downsample)
            reduced = pd.DataFrame(0.0, index=range(0, reduced_size),
              columns=dataframe.columns, dtype=np.int64)
            for i in range(0, reduced_size):
                reduced.loc[i] = dataframe[
                  i*downsample:(i+1)*downsample+1].mean()
            dataframe = reduced

        # Offset
        dataframe["time"] += xoffset

        # Plot pdist
        if pdist:
            from sklearn.neighbors import KernelDensity

            if not hasattr(subplot, "_mps_partner_subplot"):
                from .myplotspec.axes import add_partner_subplot
                add_partner_subplot(subplot, **kwargs)
            kde_kw = multi_get_copy("kde_kw", kwargs, {"bandwidth": 0.1})
            grid = kwargs.get("grid", np.linspace(0,6,100))
            kde = KernelDensity(**kde_kw)
            kde.fit(dataframe[ykey][:, np.newaxis])
            pdf = np.exp(kde.score_samples(grid[:, np.newaxis]))
            pdf /= pdf.sum()
            pdist_kw = plot_kw.copy()
            pdist_kw.update(kwargs.get("pdist_kw", {}))
            subplot._mps_partner_subplot.plot(pdf, grid, **pdist_kw)
            pdf_max = pdf.max()
            x_max = subplot._mps_partner_subplot.get_xbound()[1]
            if pdf_max > x_max / 1.25:
                subplot._mps_partner_subplot.set_xbound(0, pdf_max * 1.25)
                xticks = [0, pdf_max*0.3125, pdf_max*0.625, pdf_max*0.9375,
                          pdf_max*1.25]
                subplot._mps_partner_subplot.set_xticks(xticks)

        # Plot
        if fill_between:
            fill_between_lb_key = kwargs.get("fill_between_lb_key")
            fill_between_ub_key = kwargs.get("fill_between_ub_key")
            subplot.fill_between(dataframe["time"], 
              dataframe[fill_between_lb_key],
              dataframe[fill_between_ub_key], **fill_between_kw)
        plot = subplot.plot(dataframe["time"], dataframe[ykey], **plot_kw)[0]
        handle_kw = multi_get_copy("handle_kw", kwargs, {})
        handle_kw["mfc"] = plot.get_color()
        handle = subplot.plot([-10, -10], [-10, -10], **handle_kw)[0]
        if handles is not None and label is not None:
            handles[label] = handle

        # Plot experiment
#        if (kwargs.get("experiment", False)
#        and not hasattr(subplot, "_mps_experiment")):
#            subplot._mps_experiment = subplot.axhspan(20.04, 21.54, lw=0,
#            color=[0.7,0.7,0.7])
#            subplot._mps_partner_subplot.axhspan(20.04, 21.54, lw=0,
#            color=[0.7,0.7,0.7])
#            handles["Experiment"] = subplot.plot([-10, -10], [-10, -10],
#              color=[0.7, 0.7, 0.7], lw=2)[0]

#################################### MAIN #####################################
if __name__ == "__main__":
    TimeSeriesFigureManager().main()
