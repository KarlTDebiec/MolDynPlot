#!/usr/bin/python
# -*- coding: utf-8 -*-
#   moldynplot.dataset.SAXSDataset.py
#
#   Copyright (C) 2015-2016 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Moldynplot includes several dataset classes that build on
:class:`Dataset<myplotspec.Dataset.Dataset>` with additions
specific for molecular dynamics simulation data.

.. todo:
  - Fix separation and ordering of argument groups: input, action, output
  - Move SAXS parsing here
"""
################################### MODULES ###################################
from __future__ import absolute_import,division,print_function,unicode_literals
if __name__ == "__main__":
    __package__ = str("moldynplot")
    import moldynplot
from IPython import embed
import h5py
import numpy as np
import pandas as pd
from ..myplotspec.Dataset import Dataset
from ..myplotspec import sformat, wiprint
################################### CLASSES ###################################
class SAXSDataset(Dataset):
    """
    Represents Small-Angle X-ray Scattering Data.
    """

    def scale(self, scale, **kwargs):
        """
        Scales SAXS intensity, either by a constant or to match the
        intensity of a target dataset.

        Arguments:
          scale (float, str): If float, proportion by which to scale
            intensity; if str, path to input file to which intensity
            will be scaled, may contain environment variables
          curve_fit_kw (dict): Keyword arguments passed to
            scipy.optimize.curve_fit (scale to match target dataset
            only)
          verbose (int): Level of verbose output
          kwargs (dict): Additional keyword arguments
        """
        from os.path import expandvars, isfile
        from scipy.interpolate import interp1d
        from scipy.optimize import curve_fit
        import six

        verbose = kwargs.get("verbose", 1)

        # Scale by constant
        if (isinstance(scale, float)
        or (isinstance(scale, int) and not isinstance(scale, bool))):
            scale = float(scale)
        # Scale to match target
        elif isinstance(scale, six.string_types):
            if not isfile(expandvars(scale)):
                if verbose >= 1:
                   wiprint("scale target '{0}' ".format(scale) +
                          "not found, not scaling.")
                return

            # Prepare target
            target = self.load_dataset(infile=expandvars(scale),
                       loose=True).dataframe
            if "intensity_se" in target.columns.values:
                scale_se = True
            else:
                scale_se = False
            target_q = np.array(target.index.values, np.float64)
            target_I = np.array(target["intensity"], np.float64)
            if scale_se:
                target_Ise = np.array(target["intensity_se"], np.float64)

            # Prepare own values over x range of target
            template = self.dataframe
            template_q  = np.array(template.index.values, np.float64)
            template_I  = np.array(template["intensity"].values, np.float64)
            indexes = np.logical_and(template_q > target_q.min(),
                                     template_q < target_q.max())
            template_q = template_q[indexes]
            template_I = template_I[indexes]

            # Update target
            target_I = interp1d(target_q, target_I, kind="cubic")(template_q)
            if scale_se:
                target_Ise = interp1d(target_q, target_Ise,
                               kind="cubic")(template_q)

            def scale_I(_, a):
                return a * template_I
#            curve_fit_kw = dict(p0=(1), bounds=(0.0,0.35))
            curve_fit_kw = dict(p0=(1))     # Not clear why bounds broke
            curve_fit_kw.update(kwargs.get("curve_fit_kw", {}))
            if scale_se:
                curve_fit_kw["sigma"] = target_Ise
            scale = curve_fit(scale_I, template_q, target_I,
                      **curve_fit_kw)[0][0]
        # 'scale' argument not understood
        else:
            if verbose >= 1:
                wiprint("scale target '{0}' ".format(scale) +
                      "not understood, not scaling.")
            return

        if verbose >= 1:
            wiprint("scaling by factor of {0}".format(scale))
        self.dataframe["intensity"] *= scale
        if "intensity_se" in self.dataframe.columns.values:
            self.dataframe["intensity_se"] *= scale

        return scale


class SAXSExperimentDataset(SAXSDataset):
    """
    Represents Small Angle X-ray Scattering experimental data.
    """

    def __init__(self, scale=False, **kwargs):
        """
        Arguments:
          infile (str): Path to input file, may contain environment
            variables
          verbose (int): Level of verbose output
          kwargs (dict): Additional keyword arguments
        """
        from os.path import expandvars

        # Load
        super(SAXSExperimentDataset, self).__init__(**kwargs)
        dataframe = self.dataframe

        # Scale
        if scale:
            self.scale(scale, **kwargs)


class SAXSDiffDataset(SAXSDataset):
    """
    Represents Small Angle X-ray Scattering difference data.
    """

#    @classmethod
#    def get_cache_key(cls, dataset_classes=None, *args, **kwargs):
#        """
#        Generates tuple of arguments to be used as key for dataset
#        cache.
#
#        Arguments documented under :func:`__init__`.
#        """
#        from ..myplotspec import multi_get_copy
#
#        minuend_kw = multi_get_copy(["minuend", "minuend_kw"], kwargs, {})
#        minuend_class = dataset_classes[minuend_kw["kind"].lower()]
#        key = [cls, mask_cutoff, minuend_class.get_cache_key(**minuend_kw)]
#
#        subtrahend_kw = multi_get_copy(["subtrahend", "subtrahend_kw"],
#          kwargs, {})
#        if isinstance(subtrahend_kw, dict):
#            subtrahend_kw = [subtrahend_kw]
#        for sh_kw in subtrahend_kw:
#            sh_class = dataset_classes[sh_kw.pop("kind").lower()]
#            key.append(sh_class.get_cache_key(**sh_kw))
#
#        return tuple(key)

    def __init__(self, dataset_cache=None, **kwargs):
        """
        """
        from sys import exit
        from ..myplotspec import multi_get_copy

        self.dataset_cache = dataset_cache

        minuend_kw    = multi_get_copy(["minuend", "minuend_kw"],
                          kwargs, {})
        subtrahend_kw = multi_get_copy(["subtrahend", "subtrahend_kw"],
                          kwargs, {})
        minuend    = self.load_dataset(loose=True, **minuend_kw)
        subtrahend = self.load_dataset(loose=True, **subtrahend_kw)
        m_I    = minuend.dataframe["intensity"]
        m_I_se = minuend.dataframe["intensity_se"]
        s_I    = subtrahend.dataframe["intensity"]
        s_I_se = subtrahend.dataframe["intensity_se"]
        diff_I    = (m_I - s_I)
        diff_I_se = np.sqrt(m_I_se**2 +s_I_se**2)
        diff_I.name = "intensity"
        diff_I_se.name = "intensity_se"
        self.dataframe = pd.concat([diff_I, diff_I_se], axis=1)

#################################### MAIN #####################################
if __name__ == "__main__":
    import argparse

    # Prepare argument parser
    parser = argparse.ArgumentParser(
      description = """Processes datasets""")
    subparsers = parser.add_subparsers(
      dest        = "mode",
      description = "")

    SAXSDataset.construct_argparser(subparsers)
    SAXSExperimentDataset.construct_argparser(subparsers)
    SAXSDiffDataset.construct_argparser(subparsers)

    kwargs  = vars(parser.parse_args())
    kwargs.pop("cls")(**kwargs)
