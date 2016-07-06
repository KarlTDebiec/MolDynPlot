#!/usr/bin/python
# -*- coding: utf-8 -*-
#   myplotspec_nmr.HSQCFigureManager.py
#
#   Copyright (C) 2015 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Generates one or more HSQC figures to specifications provided in a YAML
file.
"""
################################### MODULES ###################################
from __future__ import absolute_import,division,print_function,unicode_literals
if __name__ == "__main__":
    __package__ = str("myplotspec_nmr")
    import myplotspec_nmr
from .myplotspec.FigureManager import FigureManager
################################### CLASSES ###################################
class HSQCFigureManager(FigureManager):
    """
    Manages the generation of HSQC figures.
    """
    from .myplotspec.manage_defaults_presets import manage_defaults_presets
    from .myplotspec.manage_kwargs import manage_kwargs

    defaults = """
        draw_figure:
          subplot_kw:
            autoscale_on: True
            axisbg: none
        draw_subplot:
          xlabel: $^{1}H$ (ppm)
          ylabel: $^{15}N$ (ppm)
          xticks: [5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0,10.5,11.0,11.5]
          yticks: [100,105,110,115,120,125,130,135]
          grid: True
          grid_kw:
            b: True
            linestyle: '-'
            color: [0.8,0.8,0.8]
          legend_kw:
            loc: 2
            handletextpad: 0.0
            numpoints: 1
            frameon: True
          tick_params:
            bottom: on
            top: off
            right: off
            left: on
            direction: out
          title_kw:
            verticalalignment: bottom
        draw_dataset:
          cutoff: 0.970
          plot_kw:
            zorder: 10
          handle_kw:
            marker: o
            ls: none
            mew: 0
    """
    available_presets = """
      letter:
        class: target
        inherits: letter
        draw_figure:
          left:       0.9
          sub_width:  7.8
          bottom:     0.7
          sub_height: 5.4
        draw_subplot:
          xticklabels: [5.0,'',6.0,'',7.0,'',8.0,'',9.0,'',10.0,'',11.0]
          legend: True
        draw_dataset:
          plot_kw:
            linewidths: 0.75
          handle_kw:
            ms: 8
      notebook:
        class: target
        inherits: notebook
        draw_figure:
          left:       0.7
          sub_width:  5.5
          right:      0.3
          bottom:     0.6
          sub_height: 3.8
          top:        0.3
        draw_subplot:
          xticklabels: [5.0,'',6.0,'',7.0,'',8.0,'',9.0,'',10.0,'',11.0]
          legend: True
        draw_dataset:
          plot_kw:
            linewidths: 0.5
          handle_kw:
            ms: 5
      presentation:
        class: target
        inherits: presentation
        draw_figure:
          left:       2.12
          sub_width:  6.00
          sub_height: 4.00
          bottom:     2.10
        draw_subplot:
          xticks: [6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0]
          xticklabels: [6,'',7,'',8,'',9,'',10]
          yticks: [105,110,115,120,125,130]
          yticklabels: [105,110,115,120,125,130]
          legend: False
        draw_dataset:
          plot_kw:
            linewidths: 1
    """

    @manage_defaults_presets()
    @manage_kwargs()
    def draw_dataset(self, subplot, infile, peaklist=None, sequence=None,
        label=None, handles=None, xoffset=0, yoffset=0, **kwargs):
        """
        Draws a dataset.

        Arguments:
          subplot (Axes): Axes on which to draw
          infile (str): Path to input file; NMRPipe format
          label (str, optional): Dataset label
            handles: Nascent list of dataset handles on subplot
          color (str, list, ndarray, float, optional): Dataset color
          xoffset (float, optional): Offset added to x coordinates
          yoffset (float, optional): Offset added to y coordinates
          plot_kw (dict, optional): Additional keyword arguments passed
            to subplot.plot()
          handles (OrderedDict, optional): Nascent OrderedDict of
            [labels]: handles on subplot
          kwargs (dict): Additional keyword arguments
        """
        import numpy as np
        from . import get_cmap, get_contour_levels
        from .myplotspec import get_color, multi_get_copy
        import nmrglue

        # Cheap way to invert axes without overriding draw_subplot
        if not subplot.xaxis_inverted():
            subplot.invert_xaxis()
        if not subplot.yaxis_inverted():
            subplot.invert_yaxis()

        def get_slice(collection, lower_bound, upper_bound):
            try:
                min_index = np.where(collection < lower_bound)[0][0] + 1
            except IndexError:
                min_index = collection.size
            try:
                max_index = max(np.where(collection > upper_bound)[0][-1],0)
            except IndexError:
                max_index = 0
            return slice(max_index, min_index, 1)

        # Load HSQC data
        parameters, intensity = nmrglue.pipe.read(infile)
        hydrogen  = nmrglue.pipe.make_uc(parameters, intensity,
                       dim=1).ppm_scale()
        nitrogen  = nmrglue.pipe.make_uc(parameters, intensity,
                        dim=0).ppm_scale()

        hydrogen += xoffset
        nitrogen += yoffset


#        # Load peak list
#        if peaklist is not None:
#            from . import parse_ccpnmr_peaks
#            peaklist = parse_ccpnmr_peaks(peaklist)

#        # Load sequence
#        if sequence is not None:
#            import re
#
#            sequence_dict = {}
#            with open(sequence, "r") as sequence_infile:
#                lines = sequence_infile.readlines()
#            index = 1
#            for line in lines:
#                if line.startswith(">"):
#                    index = int(line.split()[2])
#                else:
#                    for char in re.sub("[^A-Z]","", line):
#                        sequence_dict[index] = char
#                        index += 1

        # Configure plot settings
        plot_kw = multi_get_copy("plot_kw", kwargs, {})
        if "levels" not in plot_kw:
            plot_kw["levels"] = get_contour_levels(intensity, **kwargs)
        if "cmap" not in plot_kw:
            if "color" in plot_kw:
                plot_kw["color"] = get_color(plot_kw.pop("color"))
            elif "color" in kwargs:
                plot_kw["color"] = get_color(kwargs.pop("color"))
            plot_kw["cmap"] = get_cmap(plot_kw["color"])
        if label is not None:
            plot_kw["label"] = label

        # Cut data that lies outside boundaries
#        H_slice = get_slice(hydrogen, *subplot.get_xbound())
#        N_slice = get_slice(nitrogen, *subplot.get_ybound())

        # Plot
        subplot.contour(hydrogen, nitrogen, intensity, **plot_kw)
#        subplot.contour(hydrogen[H_slice], nitrogen[N_slice],
#          intensity[N_slice, H_slice], **plot_kw)
        if handles is not None and label is not None:
            handle_kw = multi_get_copy("handle_kw", kwargs, {})
            handle_kw["mfc"] = plot_kw.get("color")
            handles[label] = subplot.plot([0,0], [0,0], **handle_kw)[0]
#
#        # Draw peak labels
#        if peaklist is not None:
#            from .myplotspec.text import set_text
#            from . import three_one
#            for peak in peaklist.iterrows():
#                peak_x = peak[1].loc[2]
#                peak_y = peak[1].loc[3]
#                assign = peak[1].loc[5].strip()
#                if not assign.startswith("B"):
#                    continue
#                try:
#                    float(peak_x)
#                except:
#                    continue
#                label = "{0}{1}".format(three_one(assign[-4:-1]), assign[1:-4])
#                print(peak_x, peak_y, assign, label)
##                if sequence is not None and peak_index in sequence_dict:
##                    peak_label = "{0}{1}".format(sequence_dict[peak_index],
##                      peak_index)
#                set_text(subplot, s=label, x=peak_x, y=peak_y, fp="6r",
#                  text_kw=dict(ha="center", va="center"))

#################################### MAIN #####################################
if __name__ == "__main__":
    HSQCFigureManager().main()
