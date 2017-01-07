#!/usr/bin/python
# -*- coding: utf-8 -*-
#   moldynplot.SequenceFigureManager.py
#
#   Copyright (C) 2015-2017 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Generates one or more sequence figures to specifications in a YAML file.
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
class SequenceFigureManager(FigureManager):
    """
    Manages the generation of sequence figures.

    **Supported Presets:**

    Relaxation (``relax_4_s2``):

    .. image:: _static/gb3/relax.png
        :scale: 35
        :align: center
    """

    defaults = """
        draw_figure:
          subplot_kw:
            autoscale_on: True
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
          xlabel: Residue
          ylabel_kw:
            va: center
          tick_params:
            bottom: on
            top: off
            right: off
            left: on
            direction: out
          grid: True
          grid_kw:
            b: True
            linestyle: '-'
            color: [0.8,0.8,0.8]
          label_kw:
            zorder: 10
            horizontalalignment: left
            verticalalignment: top
        draw_dataset:
          dataset_kw:
            cls: moldynplot.dataset.SequenceDataset.SequenceDataset
            read_csv_kw:
              delimiter: "\\\\s\\\\s+"
              index_col: 0
          fill_between_kw:
            zorder: 9
          errorbar_kw:
            ls: None
            zorder: 10
          handle_kw:
            ls: none
            marker: s
            mec: black
    """
    available_presets = """
      secstruct:
        class: content
        help: Plot secondary structure as a function of sequence
        draw_subplot:
          ylabel: "% α Helix"
          yticks: [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
        draw_dataset:
          kind: ss
          column: Alpha
      r1:
        class: content
        help: Format subplot for R1 relaxation
        draw_subplot:
          ylabel: "$R_1$"
          yticks: [0,1,2,3,4,5]
        draw_dataset:
          dataset_kw:
            cls: moldynplot.dataset.SequenceDataset.RelaxDataset
          column:    "r1"
          column_se: "r1 se"
      r2:
        class: content
        help: Format subplot for R2 relaxation
        draw_subplot:
          ylabel: "$R_2$"
          yticks: [0,2,4,6,8,10,12,14,16,18,20]
        draw_dataset:
          dataset_kw:
            cls: moldynplot.dataset.SequenceDataset.RelaxDataset
          column:    "r2"
          column_se: "r2 se"
      r2/r1:
        class: content
        help: Format subplot for R2/R1 relaxation
        draw_subplot:
          ylabel: "$\\\\frac{R_2}{R_1}$"
          ylabel_fp: 12b
          yticks: [3,4,5,6,7,8,9,10,11]
        draw_dataset:
          dataset_kw:
            cls: moldynplot.dataset.SequenceDataset.RelaxDataset
          column:    "r2/r1"
          column_se: "r2/r1 se"
      hetnoe:
        class: content
        help: Format subplot for Heteronuclear NOE relaxation
        draw_subplot:
          ylabel: "Het\n\nNOE"
          yticks: [0.0,0.2,0.4,0.6,0.8,1.0]
        draw_dataset:
          dataset_kw:
            cls: moldynplot.dataset.SequenceDataset.RelaxDataset
          column:    "noe"
          column_se: "noe se"
      s2:
        class: content
        help: Format subplot for S2 order parameter
        draw_subplot:
          ylabel: "$S^2$"
          yticks: [0.0,0.2,0.4,0.6,0.8,1.0]
        draw_dataset:
          dataset_kw:
            cls: moldynplot.dataset.SequenceDataset.RelaxDataset
          column:    "s2"
          column_se: "s2 se"
      error:
        class: content
        help: Format subplot for normalized error
        draw_subplot:
          ylabel: "Normalized\\nerror"
          yticks: [0.0,0.2,0.4,0.6,0.8,1.0]
          ylabel_kw:
            rotation: vertical
        draw_dataset:
          draw_fill_between: False
          draw_errorbar:     False
          draw_plot:         True
      pre:
        class: content
        help: Format subplot for Paramagnetic Relaxation Enhancement
        draw_subplot:
          ylabel: "$\\\\frac{diamagnetic}{paramagnetic}$"
          ylabel_fp: 12b
          ylabel_kw:
            rotation: vertical
            labelpad: 12
          yticks: [0.0,0.2,0.4,0.6,0.8,1.0]
        draw_dataset:
          column:    "pre"
          column_se: "pre se"
      deltacs:
        class: content
        help: Format subplot for chemical shift change
        draw_subplot:
          ylabel: "$\\\\Delta \\\\delta$"
          ylabel_kw:
            rotation: vertical
          yticks: [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
        draw_dataset:
          column: "delta cs"
      relax_3:
        class: content
        help: Three stacked plots including R1, R2, and HetNOE
        draw_figure:
          multiplot: True
          nrows: 3
          subplots:
            0:
              preset: "r1"
            1:
              preset: "r2"
            2:
              preset: "hetnoe"
      relax_4:
        class: content
        help: Four stacked plots including R1, R2, R2/R1, and HetNOE
        draw_figure:
          multiplot: True
          nrows: 4
          subplots:
            0:
              preset: "r1"
            1:
              preset: "r2"
            2:
              preset: "r2/r1"
            3:
             preset: "hetnoe"
      relax_4_s2:
        class: content
        help: Four stacked plots including R1, R2, HetNOE, and S2
        draw_figure:
          multiplot: True
          nrows: 4
          subplots:
            0:
              preset: "r1"
            1:
              preset: "r2"
            2:
             preset: "hetnoe"
            3:
             preset: "s2"
      manuscript:
        class: target
        inherits: manuscript
        draw_figure:
          left:       0.60
          sub_width:  6.00
          wspace:     0.05
          right:      0.12
          bottom:     0.70
          sub_height: 1.00
          hspace:     0.10
          top:        0.25
          title_kw:
            top:     -0.125
          shared_legend_kw:
            left:       0.60
            sub_width:  6.00
            bottom:     0.00
            sub_height: 0.30
            legend_kw:
              ncol: 5
        draw_subplot:
          xlabel_kw:
            labelpad: 3
          ylabel_kw:
            rotation: horizontal
            labelpad: 3
          grid_kw:
            alpha: 0.3
          draw_label: True
          label_kw:
            border_lw: 1
            xabs:  0.02
            yabs: -0.03
        draw_dataset:
          fill_between_kw:
            lw: 1
          errorbar_kw:
            capsize: 0
            elinewidth: 1.5
            marker: 'o'
            mew: 0
            ms: 4
      manuscript_tight:
        class: target
        extends: manuscript
        help: Tighter formatting for two columns
        draw_figure:
          left:       0.55
          sub_width:  2.85
          right:      0.10
          top:        0.25
          title_kw:
            top:     -0.125
          shared_legend_kw:
            left:       0.55
            sub_width:  2.85
            bottom:     0.05
            sub_height: 0.30
            legend_kw:
              loc: 10
        draw_dataset:
          fill_between_kw:
            lw: 0.5
          errorbar_kw:
            capsize: 0
            elinewidth: 0.75
            marker: 'o'
            mew: 0
            ms: 2
      notebook:
        class: target
        inherits: notebook
        draw_figure:
          left:       0.70
          sub_width:  5.60
          right:      0.20
          bottom:     0.40
          sub_height: 1.00
          hspace:     0.00
          top:        0.30
        draw_subplot:
          ylabel_kw:
            labelpad: 15
        draw_dataset:
          bar_kw:
            lw: 0.5
      presentation:
        class: target
        inherits: presentation
        draw_figure:
          left:       1.50
          sub_width:  6.00
          bottom:     2.50
          sub_height: 1.20
          hspace:     0.10
          shared_legend_kw:
            left: 7.6
            sub_width: 2.00
            bottom: 2.30
            sub_height:  2.5
            legend_kw:
                loc: 2
        draw_subplot:
          xlabel_kw:
            labelpad: 3
          ylabel_kw:
            rotation: horizontal
            labelpad: 3
        draw_dataset:
          fill_between_kw:
            lw: 1
            zorder: 9
          errorbar_kw:
            capsize: 0
            elinewidth: 2
            marker: 'o'
            mew: 0
            ms: 5
          handle_kw:
            ms: 12
            mew: 2
    """

    @manage_defaults_presets()
    @manage_kwargs()
    def draw_dataset(self, subplot,
        column=None, column_se=None, label=None,
        handles=None,
        draw_fill_between=False, draw_errorbar=True, draw_plot=False,
        draw_handle=True, draw_label=False,
        verbose=1, debug=0, **kwargs):
        import numpy as np
        from .myplotspec import get_colors, multi_get_copy

        # Load data
        dataset_kw = multi_get_copy("dataset_kw", kwargs, {})
        if "infile" in kwargs:
            dataset_kw["infile"] = kwargs["infile"]
        dataframe = self.load_dataset(verbose=verbose, debug=debug,
          **dataset_kw).sequence_df
        x = np.array([filter(lambda x: x in '0123456789.', s)
              for s in dataframe.index.values], np.int)

        # Configure plot settings
        plot_kw = multi_get_copy("plot_kw", kwargs, {})
        get_colors(plot_kw, kwargs)

        # Plot fill_between
        if draw_fill_between:
            fill_between_kw = multi_get_copy("fill_between_kw", kwargs, {})
            get_colors(fill_between_kw, plot_kw)
            fb_x, fb_lb, fb_ub = [], [], []
            yse_min = kwargs.get("yse_min")
            for residue in range(x.min(), x.max()+1):
                if residue in x:
                    fb_x.extend([residue-0.5, residue+0.5])
                    index = np.argmax(x==residue)
                    if column_se is not None:
                        yse = dataframe[column_se][index]
                        if yse_min is not None and yse < yse_min:
                            yse = yse_min
                    elif yse_min is not None:
                        yse = yse_min
                    else:
                        yse = 0
                    fb_lb.extend([
                      dataframe[column][index] - yse * 1.96,
                      dataframe[column][index] - yse * 1.96])
                    fb_ub.extend([
                      dataframe[column][index] + yse * 1.96,
                      dataframe[column][index] + yse * 1.96])
                else:
                    fb_x.append(None)
                    fb_lb.append(None)
                    fb_ub.append(None)
            fb_x  = np.array(fb_x, np.float)
            fb_lb = np.array(fb_lb, np.float)
            fb_ub = np.array(fb_ub, np.float)
                
            subplot.fill_between(fb_x, fb_lb, fb_ub, **fill_between_kw)

        # Plot error bar
        if draw_errorbar:
            errorbar_kw = multi_get_copy("errorbar_kw", kwargs, {})
            get_colors(errorbar_kw, plot_kw)
            subplot.errorbar(x, y=dataframe[column],
              yerr=dataframe[column_se]*1.96,
              **errorbar_kw)

        # Plot series
        if draw_plot:
            p_x, p_y = [], []
            for residue in range(x.min(), x.max()+1):
                if residue in x:
                    p_x.extend([residue-0.5, residue+0.5])
                    index = np.argmax(x==residue)
                    p_y.extend([dataframe[column][index],
                                dataframe[column][index]])
                else:
                    p_x.append(None)
                    p_y.append(None)
            p_x = np.array(p_x, np.float)
            p_y = np.array(p_y, np.float)
            plot = subplot.plot(p_x, p_y, **plot_kw)[0]
            if verbose >= 2:
                print(column, np.nanmean(p_y))
                print(p_x)
                print(p_y)

        # Plot handle
        if draw_handle:
            handle_kw = multi_get_copy("handle_kw", kwargs, {})
            handle_kw["mfc"] = plot_kw["color"]
            handle = subplot.plot([-10, -10], [-10, -10], **handle_kw)[0]
            if handles is not None and label is not None:
                handles[label] = handle

#################################### MAIN #####################################
if __name__ == "__main__":
    SequenceFigureManager().main()
