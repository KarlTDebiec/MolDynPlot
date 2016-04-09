#!/usr/bin/python
# -*- coding: utf-8 -*-
#   moldynplot.MDGXFigureManager.py
#
#   Copyright (C) 2015-2016 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Generates one or more MDGX figures to specifications provided in a YAML
file.
"""
################################### MODULES ###################################
from __future__ import absolute_import,division,print_function,unicode_literals
import matplotlib
if __name__ == "__main__":
    __package__ = str("moldynplot")
    import moldynplot
from myplotspec.FigureManager import FigureManager
################################### CLASSES ###################################
class MDGXFigureManager(FigureManager):
    """
    Manages the generation of MDGX figures.
    """
    from .myplotspec.manage_defaults_presets import manage_defaults_presets
    from .myplotspec.manage_kwargs import manage_kwargs
    from myplotspec.debug import debug_arguments

    defaults = """
        draw_figure:
          subplot_kw:
            autoscale_on: False
          multi_tick_params:
            left: on
            right: off
            bottom: off
            top: off
        draw_subplot:
          title_kw:
            verticalalignment: bottom
          yticks: [0,1,2,3,4]
          ylabel: '$\\left|{U_{QM}-U_{MM}}\\right|$\n\n(kcal/mol)'
          tick_params:
            direction: out
            left: on
            right: off
            bottom: off
            top: off
          grid: True
          grid_kw:
            axis: y
            b: True
            color: [0.8,0.8,0.8]
            linestyle: '-'
        draw_dataset:
          dataset_kw:
            cls: moldynplot.MDGXDataset.MDGXDataset
            read_csv_kw:
              delim_whitespace: True
              names: [topology, restart, qm_energy, mm_energy]
              header: 0
    """

    available_presets = """
      manuscript:
        class: target
        inherits: manuscript
        draw_figure:
          left:       0.59
          sub_width:  6.29
          right:      0.12
          bottom:     0.30
          sub_height: 1.00
          top:        0.20
        draw_subplot:
          xlabel_kw:
            labelpad: 3
          ylabel_kw:
            labelpad: 8
          grid_kw:
            alpha: 0.3
        draw_dataset:
          violin_kw:
            points: 1000
            widths: 0.8
            showmeans: False
            showextrema: False
          body_kw:
            alpha: 1
            zorder: 10
          median_kw:
            color: black
            lw: 1
          percentile_kw:
            color: black
            lw: 0.5
            zorder: 11
          mean_kw:
            linestyle: none
            marker: o
            ms: 2
            mfc: white
            mew: 0
            zorder: 12
          edge_kw:
            alpha: 1
            zorder: 13
            lw: 1
        """

    @manage_defaults_presets()
    @manage_kwargs()
    def draw_dataset(self, subplot,
        label=None,
        handles=None,
        draw_body=True, draw_percentile=True, draw_mean=True,
        draw_edge=True,
        verbose=1, debug=0, **kwargs):
        """
        .. todo:
            - Use violinplot function only to get vertices, and then
              draw manually, rather than calculating KDE twice. May
              alternatively calculate KDE myself in MDGXDataset.
        """
        from .myplotspec import get_color, get_colors, multi_get_copy
        import numpy as np

        # Load data
        dataset_kw = multi_get_copy("dataset_kw", kwargs, {})
        if "infile" in kwargs:
            dataset_kw["infile"] = kwargs["infile"]
        dataset = self.load_dataset(verbose=verbose, debug=debug, **dataset_kw)
        dataframe = dataset.dataframe

        # Configure plot settings
        plot_kw = multi_get_copy("plot_kw", kwargs, {})
        get_colors(plot_kw, kwargs)
        violin_kw = multi_get_copy("violin_kw", kwargs, {})
        get_colors(violin_kw, plot_kw)
        y = [d["error"] for d in dataset.selections]
        x = kwargs.get("x", np.arange(0.5,len(dataset.selections)+0.5,1))
        if draw_percentile:
            percentile_kw = plot_kw.copy()
            percentile_kw.update(kwargs.get("percentile_kw", {}))
            median_kw = percentile_kw.copy()
            median_kw.update(kwargs.get("median_kw", {}))
        if draw_mean:
            mean_kw = plot_kw.copy()
            mean_kw.update(kwargs.get("mean_kw", {}))

        # Plot bodies
        if draw_body:
            body_kw = plot_kw.copy()
            body_kw.update(kwargs.get("body_kw", {}))
            colors = multi_get_copy(["color", "colors", "facecolor",
              "facecolors"], body_kw)
            if colors is not None:
                colors = map(get_color, colors)

            bodies = subplot.violinplot(y, x, **violin_kw)["bodies"]
            for i, body in enumerate(bodies):
                body.set_edgecolors("none")
                if isinstance(colors, list):
                    body.set_facecolors(colors[i])
                elif colors is not None:
                    body.set_facecolor(colors)
                if "alpha" in body_kw:
                    body.set_alpha(body_kw["alpha"])
                if "zorder" in body_kw:
                    body.set_zorder(body_kw["zorder"])
                vertices = body.get_paths()[0].vertices

                # Plot percentiles
                if draw_percentile:
                    prc_y = np.percentile(y[i], 25)
                    prc_xmin = vertices[np.abs(
                      vertices[:vertices.shape[0],1]-prc_y).argmin(),0]
                    subplot.plot([prc_xmin, x[i]+(x[i]-prc_xmin)],
                      [prc_y, prc_y], **percentile_kw)
                    prc_y = np.percentile(y[i], 50)
                    prc_xmin = vertices[np.abs(
                      vertices[:vertices.shape[0],1]-prc_y).argmin(),0]
                    subplot.plot([prc_xmin, x[i]+(x[i]-prc_xmin)],
                      [prc_y, prc_y], **median_kw)
                    prc_y = np.percentile(y[i], 75)
                    prc_xmin = vertices[np.abs(
                      vertices[:vertices.shape[0],1]-prc_y).argmin(),0]
                    subplot.plot([prc_xmin, x[i]+(x[i]-prc_xmin)],
                      [prc_y, prc_y], **percentile_kw)

                if draw_mean:
                    subplot.plot([x[i]], [y[i].mean()], **mean_kw)

        # Plot edges
        if draw_edge:
            edge_kw = plot_kw.copy()
            edge_kw.update(kwargs.get("edge_kw", {}))
            colors = multi_get_copy(["color", "colors", "edgecolor",
              "edgecolors"], edge_kw)
            if colors is not None:
                colors = map(get_color, colors)
            bodies = subplot.violinplot(y, x, **violin_kw)["bodies"]
            for i, body in enumerate(bodies):
                body.set_facecolors("none")
                if isinstance(colors, list):
                    body.set_edgecolors(colors[i])
                elif colors is not None:
                    body.set_edgecolor(colors)
                if "alpha" in edge_kw:
                    body.set_alpha(edge_kw["alpha"])
                if "zorder" in edge_kw:
                    body.set_zorder(edge_kw["zorder"])
                if "linewidth" in edge_kw:
                    body.set_linewidth(edge_kw["linewidth"])
                elif "lw" in edge_kw:
                    body.set_linewidth(edge_kw["lw"])
                body.get_paths()[0].codes[-1] = 1 # Open closed polygon

#################################### MAIN #####################################
if __name__ == "__main__":
    MDGXFigureManager().main()
