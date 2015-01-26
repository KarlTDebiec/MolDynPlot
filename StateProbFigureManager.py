#!/usr/bin/python
# -*- coding: utf-8 -*-
#   myplotspec_sim.StateProbFigureManager.py
#
#   Copyright (C) 2015 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Generates one or more state probability figures to specifications in a
YAML file.
"""
################################### MODULES ###################################
from __future__ import absolute_import,division,print_function,unicode_literals
if __name__ == "__main__":
    __package__ = str("myplotspec_sim")
    import myplotspec_sim
from .myplotspec.FigureManager import FigureManager
################################### CLASSES ###################################
class StateProbFigureManager(FigureManager):
    """
    Class to manage the generation of probability distribution figures

    Attributes:
      defaults (str, dict): Default arguments to :func:`draw_report`,
        :func:`draw_figure`, :func:`draw_subplot`, and
        :func:`draw_dataset` functions, in yaml format. Outer level (of
        indentation or keys) provides function names, and inner level
        provides default arguments to each function::

          defaults = \"\"\"
              method_1:
                method_1_arg_1: 1000
                method_1_arg_2: abcd
              method_2
                method_2_arg_1: 2000
                method_2_arg_2: efgh
              ...
          \"\"\"

      presets (str, dict): Available sets of preset arguments to
        :func:`draw_report`, :func:`draw_figure`, :func:`draw_subplot`,
        and :func:`draw_dataset` functions, in yaml format. Outer level
        (of indentation or keys) provides preset names, middle level
        provides function names, and inner level provides arguments to
        pass to each function when preset is active::

          presets = \"\"\"
            preset_1:
              method_1:
                method_1_arg_1: 1001
                method_1_arg_2: abcde
              method_2
                method_2_arg_1: 2001
                method_2_arg_2: efghi
            preset_2:
              method_1:
                method_1_arg_1: 1002
                method_1_arg_2: abcdef
              method_2
                method_2_arg_1: 2002
                method_2_arg_2: efghij
          \"\"\"
    """

    from .myplotspec.manage_defaults_presets import manage_defaults_presets
    from .myplotspec.manage_kwargs import manage_kwargs

    defaults = """
        draw_figure:
          fig_width:     5.00
          left:          1.20
          sub_width:     3.25
          right:         0.55
          fig_height:    7.50
          top:           1.35
          sub_height:    3.25
          bottom:        2.90
          subplot_kw:
            autoscale_on: False
          shared_legend:
            left:       1.20
            sub_width:  3.25
            sub_height: 2.30
            bottom:     0.50
            legend_kw:
              frameon:      False
              labelspacing: 0.5
              loc:          9
            legend_lw:  5
        draw_subplot:
          ylabel:       $P_{bound}$
          xticklabels:  []
          yticks:       [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
          yticklabels:  [0.0,"",0.2,"",0.4,"",0.6,"",0.8,"",1.0]
          tick_params:
            bottom:     off
            top:        off
          lw:           2
        draw_dataset:
          plot_kw:
            lw:         2.0
            align:      center
            width:      0.6
            ecolor:     black
            capsize:    4
            error_kw:
              elinewidth: 2.0
              capthick:   2.0
    """
    presets  = """
      presentation:
        draw_figure:
          fig_width:     5.00
          left:          1.20
          sub_width:     3.25
          right:         0.55
          fig_height:    7.50
          top:           1.35
          sub_height:    3.25
          bottom:        2.90
          title_fp:     24r
          label_fp:     24r
          shared_legend:
            left:       1.20
            sub_width:  3.24
            sub_height: 2.30
            bottom:     0.50
            legend_fp:  14r
            legend_kw:
              frameon:      False
              labelspacing: 0.5
              loc:          9
        draw_subplot:
          title_fp:     24r
          label_fp:     24r
          tick_fp:      18r
          legend_fp:    18r
          lw:           2
        draw_dataset:
          plot_kw:
              lw:       2
      presentation_three:
        draw_figure:
          ncols:      3
          fig_width:  10.00
          left:        1.00
          sub_width:   2.50
          wspace:      0.50
          right:       0.50
          fig_height:  7.50
          bottom:      3.00
          sub_height:  2.50
          top:         2.00
          subplots:
            1:
              ylabel: ""
            2:
              ylabel: ""
          shared_legend:
            left:       6.80
            sub_width:  2.50
            sub_height: 2.50
            bottom:     0.50
            legend_fp:  14r
            legend_kw:
              loc:          2
              labelspacing: 0.5
        draw_subplot:
          title_fp:     16r
          label_fp:     24r
          tick_fp:      16r
          legend_fp:    18r
          lw:           2
        draw_dataset:
          plot_kw:
            lw:         2
    """

    @manage_defaults_presets()
    @manage_kwargs()
    def draw_subplot(self, subplot, experiment=None, title=None, legend=None,
        shared_handles=None, **in_kwargs):
        """
        Draws a subplot.

        Subplot will typically plot one or more datasets, whose
        specifications are provided in a dict structured as follows::

            datasets = {
                'all': {
                    'lw': 2,
                    ...
                },
                '0': {
                    'label':  'Dataset 1',
                    'infile': '/path/to/dataset_1.txt',
                    'color':  'red',
                    ...
                },
                '1': {
                    'label':  'Dataset 2',
                    'infile': '/path/to/dataset_2.txt',
                    'color':  'green',
                    ...
                },
                '2': {
                    'label':  'Dataset 3',
                    'infile': '/path/to/dataset_3.txt',
                    'color':  'blue',
                    ...
                },
                ...
            }

        The values stored at each integer key (0-indexed) provide the
        arguments to be passed to :func:`draw_dataset` for each of a
        series of datasets. Values stored at 'all' are passed to each
        dataset, but overridden by values specific to that dataset.

        Subplot may be formatted by adjusting or labeling the x and y
        axes, or drawing a title or a legend.

        Arguments:
          subplot (Axes): Axes on which to act
          experiment (list): Confidence interval of experimental y
            value; region shaded in gray
          datasets (dict): Dataset specifications
          title (str, optional): Subplot title
          legend (bool, optional): Draw legend on subplot
          shared_handles (OrderedDict, optional): Nascent OrderedDict of
            [labels]:handles shared among subplots of host figure; used
            to draw shared legend
          in_kwargs (dict): Additional keyword arguments
        """
        from collections import OrderedDict
        from .myplotspec.axes import set_xaxis, set_yaxis
        from .myplotspec.legend import set_legend
        from .myplotspec.text import set_title

        # Format
        set_xaxis(subplot, **in_kwargs)
        set_yaxis(subplot, **in_kwargs)
        if str(title) != "None":
            set_title(subplot, title = title, **in_kwargs)

        # Configure and plot datasets
        handles         = OrderedDict()
        dataset_specs   = in_kwargs.pop("datasets", {})
        dataset_indexes = sorted([i for i in dataset_specs.keys()
                            if str(i).isdigit()])

        # Plot experimental data
        if experiment is not None:
            subplot.axhspan(experiment[0], experiment[1], lw = 0,
              color = [0.7, 0.7, 0.7])
            handles["Experiment"] = subplot.plot([-10, -10], [-10, -10],
              color = [0.7, 0.7, 0.7], lw = 5)[0]

        for i in dataset_indexes:
            out_kwargs              = {}
            out_kwargs["debug"]     = in_kwargs.get("debug",    False)
            out_kwargs["preset"]    = in_kwargs.get("preset",    [])
            out_kwargs["yaml_dict"] = in_kwargs.get("yaml_dict", {})
            out_kwargs["yaml_keys"] = [key
              for key2 in [[key3 + ["datasets", "all"],
                            key3 + ["datasets", i]]
              for key3 in in_kwargs.get("yaml_keys")]
              for key  in key2]
            out_kwargs["x"]         = i + 0.75
            self.draw_dataset(subplot = subplot, handles = handles,
              **out_kwargs)

        # Draw legend
        if str(legend) != "None":
            set_legend(subplot, handles = handles, **in_kwargs)
        if str(shared_handles) != "None":
            for label, handle in handles.items():
                if label not in shared_handles:
                    shared_handles[label] = handle

    @manage_defaults_presets()
    @manage_kwargs()
    def draw_dataset(self, subplot, x=None, label="", handles=None,
        debug=False, **kwargs):
        """
        Draws a dataset.

        Arguments:
          subplot (Axes): Axes on which to draw
          x (float): X coordinate of bar
          label (str, optional): Dataset label
          color (str, list, ndarray, float, optional): Dataset color
          plot_kw (dict, optional): Additional keyword arguments passed
            to subplot.plot()
          handles (OrderedDict, optional): Nascent OrderedDict of
            [labels]: handles on subplot
          kwargs (dict): Additional keyword arguments
        """
        from .myplotspec import get_color
        from .H5Dataset import H5Dataset

        plot_kw = kwargs.get("plot_kw", {})
        if "color" in plot_kw:
            plot_kw["color"] = get_color(plot_kw["color"])
        elif "color" in kwargs:
            plot_kw["color"] = get_color(kwargs.pop("color"))

        # Load dataset, x, and y
        if x is None:
            if handles is not None:
                x = 1 + len(handles)
            else:
                x = 1
        if not ("P unbound" in kwargs and "P unbound se" in kwargs):
            dataset = H5Dataset(
                        default_address = "assign/stateprobs",
                        default_key = "pbound", **kwargs)
            y    = 1.0 - dataset.datasets["pbound"]["P unbound"][0]
            yerr = dataset.datasets["pbound"]["P unbound se"][0] * 1.96
        else:
            y    = 1.0 - kwargs.pop("P unbound")
            yerr = kwargs.pop("P unbound se") * 1.96

        # Plot
        bar    = subplot.bar(x, y, yerr = yerr, **plot_kw)
        color  = bar.patches[0].get_facecolor()
        handle = subplot.plot([-10, -10], [-10, -10], color = color)[0]
        if handles is not None:
            handles[label] = handle

#################################### MAIN #####################################
if __name__ == "__main__":
    StateProbFigureManager().main()
