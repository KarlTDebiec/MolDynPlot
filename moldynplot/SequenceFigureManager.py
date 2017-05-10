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
from __future__ import (absolute_import, division, print_function,
    unicode_literals)

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
            lw: 0
            zorder: 1
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
          y_key: Alpha
      r1:
        class: content
        help: Format subplot for R1 relaxation
        draw_subplot:
          ylabel: "$R_1$"
          yticks: [0,1,2,3,4,5]
        draw_dataset:
          dataset_kw:
            cls: moldynplot.dataset.RelaxDataset.RelaxDataset
          y_key:    "r1"
          yse_key: "r1 se"
      r2:
        class: content
        help: Format subplot for R2 relaxation
        draw_subplot:
          ylabel: "$R_2$"
          yticks: [0,2,4,6,8,10,12,14,16,18,20]
        draw_dataset:
          dataset_kw:
            cls: moldynplot.dataset.RelaxDataset.RelaxDataset
          y_key:    "r2"
          yse_key: "r2 se"
      r2/r1:
        class: content
        help: Format subplot for R2/R1 relaxation
        draw_subplot:
          ylabel: "$\\\\frac{R_2}{R_1}$"
          ylabel_fp: 12b
          yticks: [3,4,5,6,7,8,9,10,11]
        draw_dataset:
          dataset_kw:
            cls: moldynplot.dataset.RelaxDataset.RelaxDataset
          y_key:    "r2/r1"
          yse_key: "r2/r1 se"
      hetnoe:
        class: content
        help: Format subplot for Heteronuclear NOE relaxation
        draw_subplot:
          ylabel: "Het\n\nNOE"
          yticks: [0.0,0.2,0.4,0.6,0.8,1.0]
        draw_dataset:
          dataset_kw:
            cls: moldynplot.dataset.RelaxDataset.RelaxDataset
          y_key:    "noe"
          yse_key: "noe se"
      s2:
        class: content
        help: Format subplot for S2 order parameter
        draw_subplot:
          ylabel: "$S^2$"
          yticks: [0.0,0.2,0.4,0.6,0.8,1.0]
        draw_dataset:
          dataset_kw:
            cls: moldynplot.dataset.RelaxDataset.RelaxDataset
          y_key:    "s2"
          yse_key: "s2 se"
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
      pre_r20_r2:
        class: content
        help: Format subplot for paramagnetic relaxation enhancement r20 / r2
        draw_subplot:
          ylabel: "$\\\\frac{r_{2,0}}{r_2}$"
          yticks: [0.0,0.2,0.4,0.6,0.8,1.0]
        draw_dataset:
          y_key:    "r20/r2"
          yse_key: "r20/r2 se"
      pre_I_I0:
        class: content
        help: Format subplot for paramagnetic relaxation enhancement I/I0
        draw_subplot:
          ylabel: "$\\\\frac{I}{I_0}$"
          yticks: [0.0,0.2,0.4,0.6,0.8,1.0]
        draw_dataset:
          y_key:    "I/I0"
          yse_key: "I/I0 se"
      pre_gamma2:
        class: content
        help: Format subplot for paramagnetic relaxation enhancement Γ2
        draw_subplot:
          ylabel: "$Γ_2$"
          yticks: [0,20,40,60,80,100]
        draw_dataset:
          y_key:    "rho2"
          yse_key: "rho2 se"
      deltacs:
        class: content
        help: Format subplot for chemical shift difference
        draw_subplot:
          ylabel: "$\\\\Delta \\\\delta$"
          ylabel_kw:
            rotation: vertical
          yticks: [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
        draw_dataset:
          y_key: "delta cs"
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
          draw_label: True
          label_kw:
            border_lw: 1
            xabs:  0.02
            yabs: -0.03
        draw_dataset:
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
    def draw_dataset(self, subplot, y_key=None, yse_key=None, yse_min=None,
      ylb_key=None, yub_key=None, label=None, handles=None,
      draw_fill_between=False, draw_errorbar=False, draw_plot=False,
      draw_handle=True, **kwargs):
        """
        Draws a dataset on a subplot

        Loaded datasets should have attribue `sequence_df`

        Args:
            subplot (Axes): :class:`Axes<matplotlib.axes.Axes>` on
              which to draw
            dataset_kw (dict): Keyword arguments passed to
              :meth:`load_dataset
              <myplotspec.FigureManager.FigureManager.load_dataset>`
            plot_kw (dict): Keyword arguments passed to methods of
              :class:`Axes<matplotlib.axes.Axes>`
            y_key (string): Name of column containing y
            yse_key (string): Name of column containing y standard error
            yse_min (float): Minimum y standard error (fill_between only);
              useful for making sure filled area does not become too narrow
              to see clearly
            label (string): Label of dataset
            handles (OrderedDcit): Nascent OrderedDict of handles
            draw_fill_between (bool): Draw fill_between element
            draw_errorbar (bool): Draw errorbar element
            draw_plot (bool): Draw plot element
            draw_handle (bool): Draw additional handle on plot for legend
            verbose (int): Level of verbose output
          kwargs (dict): Additional keyword arguments
        """
        from warnings import warn
        import numpy as np
        from .myplotspec import get_colors, multi_get_copy

        # Process arguments and load data
        verbose = kwargs.get("verbose", 1)
        dataset_kw = multi_get_copy("dataset_kw", kwargs, {})
        if "infile" in kwargs:
            dataset_kw["infile"] = kwargs["infile"]
        dataset = self.load_dataset(verbose=verbose, **dataset_kw)
        if not hasattr(dataset, "sequence_df"):
            warn("""Dataset does not have necessary attribute 'sequence_df',
            skipping.""")
            return
        df = dataset.sequence_df
        x = np.array(
          [filter(lambda x: x in "0123456789.", s) for s in df.index.values],
          np.int)

        # Verbose output
        if verbose >= 2:
            print(df)
            print(x)

        # Configure plot settings
        plot_kw = multi_get_copy("plot_kw", kwargs, {})
        get_colors(plot_kw, kwargs)

        # Plot fill_between
        if draw_fill_between:
            if y_key is None:
                warn("""'draw_fill_between' is enabled but the necessary
                parameter 'y_key' has not been provided, skipping.""")
            # //@formatter:off
            elif ((ylb_key is None or yub_key is None)
            and yse_key is None
            and yse_min is None):
                warn("""'draw_fill_between' is enabled but the necessary
                parameters ('ylb_key' and 'yub_key'), 'yse_key', or 'ymin_key'
                has not been provided, skipping.""")
            # //@formatter:on
            else:
                fill_between_kw = multi_get_copy("fill_between_kw", kwargs, {})
                get_colors(fill_between_kw, plot_kw)
                fb_x, fb_lb, fb_ub = [], [], []

                for residue in range(x.min(), x.max() + 1):
                    if residue in x:
                        index = np.argmax(x == residue)
                        fb_x.extend([residue - 0.5, residue + 0.5])
                        # Priority given to explicit upper and lower bounds
                        if ylb_key is not None and yub_key is not None:
                            fb_lb.extend(
                              [df[ylb_key][index], df[ylb_key][index]])
                            fb_ub.extend(
                              [df[yub_key][index], df[yub_key][index]])
                        # Otherwise try for 95% CI using standard error
                        else:
                            if yse_key is not None:
                                yse = df[yse_key][index]
                                if yse_min is not None and yse < yse_min:
                                    yse = yse_min
                            elif yse_min is not None:
                                yse = yse_min
                            else:
                                yse = 0
                            fb_lb.extend([df[y_key][index] - yse * 1.96,
                                df[y_key][index] - yse * 1.96])
                            fb_ub.extend([df[y_key][index] + yse * 1.96,
                                df[y_key][index] + yse * 1.96])
                    # Add explicit gaps in sequence if no data is available
                    else:
                        fb_x.append(None)
                        fb_lb.append(None)
                        fb_ub.append(None)
            fb_x = np.array(fb_x, np.float)
            fb_lb = np.array(fb_lb, np.float)
            fb_ub = np.array(fb_ub, np.float)

            subplot.fill_between(fb_x, fb_lb, fb_ub, **fill_between_kw)

        # Plot error bar
        if draw_errorbar:
            if y_key is None:
                warn("""'draw_errorbar' is enabled but the necessary
                parameter 'y_key' has not been provided, skipping.""")
            elif yse_key is None:
                warn("""'draw_errorbar' is enabled but the necessary
                parameter 'yse_key' has not been provided,
                skipping.""")
            else:
                errorbar_kw = multi_get_copy("errorbar_kw", kwargs, {})
                get_colors(errorbar_kw, plot_kw)
                subplot.errorbar(x, y=df[y_key], yerr=df[yse_key] * 1.96,
                  **errorbar_kw)

        # Plot series
        if draw_plot:
            if y_key is None:
                warn("""'draw_plot' is enabled but the necessary
                parameter 'y_key' has not been provided, skipping.""")
            else:
                p_x, p_y = [], []
                stepline = kwargs.get("stepline", True)
                gapline = kwargs.get("gapline", False)
                for residue in range(x.min(), x.max() + 1):
                    if residue in x:
                        index = np.argmax(x == residue)
                        if stepline:
                            p_x.extend([residue - 0.5, residue + 0.5])
                            p_y.extend([df[y_key][index], df[y_key][index]])
                        else:
                            p_x.append(residue)
                            p_y.append(df[y_key][index])
                        if gapline:
                            p_x.append(None)
                            p_y.append(None)
                    else:
                        p_x.append(None)
                        p_y.append(None)
                p_x = np.array(p_x, np.float)
                p_y = np.array(p_y, np.float)
                plot = subplot.plot(p_x, p_y, **plot_kw)[0]
                if verbose >= 2:
                    print(y_key, np.nanmean(p_y))
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
