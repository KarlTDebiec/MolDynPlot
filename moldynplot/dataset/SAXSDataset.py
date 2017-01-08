#!/usr/bin/python
# -*- coding: utf-8 -*-
#   moldynplot.dataset.SAXSDataset.py
#
#   Copyright (C) 2015-2017 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Represents Small-Angle X-ray Scattering (SAXS) data
"""
################################### MODULES ###################################
from __future__ import (absolute_import, division, print_function,
    unicode_literals)

if __name__ == "__main__":
    __package__ = str("moldynplot.dataset")
    import moldynplot.dataset
import numpy as np
import pandas as pd
from ..myplotspec.Dataset import Dataset
from ..myplotspec import wiprint


################################### CLASSES ###################################
class SAXSDataset(Dataset):
    """
    Represents Small-Angle X-ray Scattering (SAXS) data
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

        # Processs arguments
        verbose = kwargs.get("verbose", 1)

        # Scale by constant
        if (isinstance(scale, float) or (
              isinstance(scale, int) and not isinstance(scale, bool))):
            scale = float(scale)
        # Scale to match target
        elif isinstance(scale, six.string_types):
            if not isfile(expandvars(scale)):
                if verbose >= 1:
                    wiprint("scale target '{0}' ".format(
                      scale) + "not found, not scaling.")
                return

            # Prepare target
            target = self.load_dataset(infile=expandvars(scale), loose=True).df
            if "intensity se" in target.columns.values:
                scale_se = True
            else:
                scale_se = False
            target_q = np.array(target.index.values, np.float64)
            target_I = np.array(target["intensity"], np.float64)
            if scale_se:
                target_Ise = np.array(target["intensity se"], np.float64)

            # Prepare own values over x range of target
            template = self.df
            template_q = np.array(template.index.values, np.float64)
            template_I = np.array(template["intensity"].values, np.float64)
            indexes = np.logical_and(template_q > target_q.min(),
              template_q < target_q.max())
            template_q = template_q[indexes]
            template_I = template_I[indexes]

            # Update target
            target_I = interp1d(target_q, target_I, kind="cubic")(template_q)
            if scale_se:
                target_Ise = interp1d(target_q, target_Ise, kind="cubic")(
                  template_q)

            def scale_I(_, a):
                return a * template_I

            #            curve_fit_kw = dict(p0=(1), bounds=(0.0,0.35))
            curve_fit_kw = dict(p0=(1))  # Not clear why bounds broke
            curve_fit_kw.update(kwargs.get("curve_fit_kw", {}))
            if scale_se:
                curve_fit_kw["sigma"] = target_Ise
            scale = \
                curve_fit(scale_I, template_q, target_I, **curve_fit_kw)[0][0]
        # 'scale' argument not understood
        else:
            if verbose >= 1:
                wiprint("scale target '{0}' ".format(
                  scale) + "not understood, not scaling.")
            return

        if verbose >= 1:
            wiprint("scaling by factor of {0}".format(scale))
        self.df["intensity"] *= scale
        if "intensity se" in self.df.columns.values:
            self.df["intensity se"] *= scale

        return scale


#################################### MAIN #####################################
if __name__ == "__main__":
    SAXSDataset.main()