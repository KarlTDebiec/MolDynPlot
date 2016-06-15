# -*- coding: utf-8 -*-
#   moldynplot.Dataset.py
#
#   Copyright (C) 2015-2016 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Manages Moldynplot datasets.
"""
################################### MODULES ###################################
from __future__ import absolute_import,division,print_function,unicode_literals
import h5py
import numpy as np
import pandas as pd
from .myplotspec.Dataset import Dataset
################################### CLASSES ###################################
class TimeSeriesDataset(Dataset):
    """
    Manages time series datasets.
    """

    def __init__(self, downsample=None, calc_pdist=False, **kwargs):
        """
        Initializes dataset.

        Arguments:
          infile (str): Path to input file, may contain environment
            variables
          usecols (list): Columns to select from DataFrame, once
            dataframe has already been loaded
          dt (float): Time interval between points; units unspecified
          toffset (float): Time offset to be added to all points (i.e.
            time of first point)
          downsample (int): Interval by which to downsample points
          downsample_mode (str): Method of downsampling; may be 'mean'
            or 'mode'
          pdist (bool): Calculate probability distribution
          pdist_key (str): Column of which to calculate probability
            distribution
          kde_kw (dict): Keyword arguments passed to
            sklearn.neighbors.KernelDensity; key argument is 'bandwidth'
          grid (ndarray): Grid on which to calculate probability
            distribution
          verbose (int): Level of verbose output
          kwargs (dict): Additional keyword arguments

          .. todo:
            - 'y' argument does not belong here; make sure it is
              removable
            - move downsampling to function
            - Calculate pdist for multiple columns
            - Make pdist a dataframe rather than explicit x and y
            - Calculate pdist using histogram
            - Verbose pdist
        """

        # Arguments
        verbose = kwargs.get("verbose", 1)

        # Load
        super(TimeSeriesDataset, self).__init__( **kwargs)
        timeseries = self.timeseries = self.dataframe
        timeseries.index.name = "time"

        if "usecols" in kwargs:
            timeseries = timeseries[timeseries.columns[kwargs.pop("usecols")]]

        # Convert from frame index to time
        if "dt" in kwargs:
            timeseries.index *= kwargs.pop("dt")

        # Offset time
        if "toffset" in kwargs:
            timeseries.index += kwargs.pop("toffset")
#
#        # Store y, if applicable
#        if "y" in kwargs:
#            self.y = kwargs.pop("y")

        # Downsample
        if downsample:
            self.downsample(downsample, **kwargs)

        if calc_pdist:
            from sklearn.neighbors import KernelDensity

            pdist_key = kwargs.get("pdist_cols", dataframe.columns.values[0])

            kde_kw = kwargs.get("pdist_kde_kw",  {"bandwidth": 0.1})
            grid = kwargs.get("pdist_grid", np.linspace(
                dataframe.values.min(), dataframe.values.max(), 100))
            kde = KernelDensity(**kde_kw)
            kde.fit(dataframe[pdist_key][:, np.newaxis])
            pdf = np.exp(kde.score_samples(grid[:, np.newaxis]))
            pdf /= pdf.sum()
            self.pdist_x = grid
            self.pdist_y = pdf

    def downsample(self, downsample, downsample_mode="mean", **kwargs):
        """
        Downsamples time series.

        Arguments:
          downsample (int): Interval by which to downsample points
          downsample_mode (str): Method of downsampling; may be 'mean'
            or 'mode'
          verbose (int): Level of verbose output
          kwargs (dict): Additional keyword arguments
        """
        from scipy.stats.mstats import mode

        # Arguments
        verbose = kwargs.get("verbose", 1)
        timeseries=self.timeseries

        # Truncate dataset
        reduced = timeseries.values[
          :timeseries.shape[0] - (timeseries.shape[0] % downsample),:]
        new_shape = (int(reduced.shape[0]/downsample), downsample,
          reduced.shape[1])
        index = np.reshape(timeseries.index.values[
          :timeseries.shape[0]-(timeseries.shape[0] % downsample)],
          new_shape[:-1]).mean(axis=1)
        reduced = np.reshape(reduced, new_shape)

        # Downsample
        if downsample_mode == "mean":
            if verbose >= 1:
                print("downsampling by factor of {0} using mean".format(
                  downsample))
            reduced = np.squeeze(reduced.mean(axis=1))
        elif downsample_mode == "mode":
            if verbose >= 1:
                print("downsampling by factor of {0} using mode".format(
                  downsample))
            reduced = np.squeeze(mode(reduced, axis=1)[0])

        # Store downsampled time series
        reduced = pd.DataFrame(data=reduced, index=index,
          columns=timeseries.columns.values)
        reduced.index.name = "time"
        timeseries = self.timeseries = reduced

        return timeseries

    def calc_error(self, error_mode="std", **kwargs):
        """
        Calculates standard error using time series data.

        .. todo:
            - Support breaking into sections (essentially downsampling,
              then calculating standard error
        """
        from .fpblockaverager.FPBlockAverager import FPBlockAverager

        # Arguments
        verbose = kwargs.get("verbose", 1)
        timeseries=self.timeseries

        # Calculate standard error
        if error_mode == "std":
            if verbose >= 1:
                print("calculating standard error using standard deviation")
            se = timeseries.std()
        elif error_mode == "block":
            if verbose >= 1:
                print("calculating standard error using block averaging")
            ba = FPBlockAverager(timeseries, **kwargs)
            se = ba.parameters.loc["exp", "a (se)"]
        else:
            if verbose >= 1:
                print("error_mode '{0}' not understood, ".format(scale) +
                      "must be one of 'std', 'block'; not calculating error.")
            return

        return se

class SAXSDataset(Dataset):
    """
    Manages Small Angle X-ray Scattering datasets.
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
                    print("scale target '{0}' ".format(scale) +
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
                print("scale target '{0}' ".format(scale) +
                      "not understood, not scaling.")
            return

        if verbose >= 1:
            print("scaling by factor of {0}".format(scale))
        self.dataframe["intensity"] *= scale
        if "intensity_se" in self.dataframe.columns.values:
            self.dataframe["intensity_se"] *= scale

        return scale

class NatConDataset(TimeSeriesDataset):
    """
    Manages native contact datasets.
    """

    def __init__(self, downsample=None, calc_pdist=True, **kwargs):
        """
        Initializes dataset.

        Arguments:
          infile (str): Path to input file, may contain environment
            variables
          usecols (list): Columns to select from DataFrame, once
            dataframe has already been loaded
          dt (float): Time interval between points; units unspecified
          toffset (float): Time offset to be added to all points (i.e.
            time of first point)
          cutoff (float): Minimum distance within which a contact is
            considered to be formed
          downsample (int): Interval by which to downsample points using
            mode
          pdist (bool): Calculate probability distribution
          verbose (int): Level of verbose output
          kwargs (dict): Additional keyword arguments
        """
        verbose = kwargs.get("verbose", 1)

        # Load
        super(NatConDataset, self).__init__(**kwargs)
        dataframe = self.dataframe
        n_contacts = self.dataframe.shape[1]

        # Convert minimum distances to percent native contacts
        cutoff = kwargs.get("cutoff", 5.5)
        percent = pd.DataFrame(data=(dataframe.values <= cutoff).sum(axis=1)
          / dataframe.shape[1], index=dataframe.index,
          columns=["percent_native_contacts"])
        dataframe = self.dataframe = percent

        # Downsample; flag included in function definition to prevent
        #   superclass from downsampling before applying cutoff
        if downsample is not None:
            from scipy.stats.mstats import mode

            if verbose >= 1:
                print("downsampling by factor of {0} using mode".format(
                  downsample))

            reduced = dataframe.values[
              :dataframe.shape[0]-(dataframe.shape[0] % downsample),:]
            new_shape=(int(reduced.shape[0]/ downsample),
                downsample, reduced.shape[1])
            index = np.reshape(dataframe.index.values[
              :dataframe.shape[0]-(dataframe.shape[0] % downsample)],
              new_shape[:-1]).mean(axis=1)
            reduced = np.reshape(reduced, new_shape)
            reduced = np.squeeze(mode(reduced, axis=1)[0])
            reduced = pd.DataFrame(data=reduced, index=index,
              columns=dataframe.columns.values)
            reduced.index.name = "time"
            dataframe = self.dataframe = reduced

        # Calculate probability distribution
        if calc_pdist:
            if verbose >= 1:
                print("calculating probability distribution using histogram")
            bins = np.linspace(0-((1/n_contacts)/2), 1+((1/n_contacts)/2),
              n_contacts+1)
            pdist, _ = np.histogram(self.dataframe.values, bins)
            pdist    =  np.array(pdist, np.float) / pdist.sum()
            pdist_x = np.zeros(bins.size*2)
            pdist_y = np.zeros(bins.size*2)
            pdist_x[::2]    = pdist_x[1::2]   = bins
            pdist_y[1:-1:2] = pdist_y[2:-1:2] = pdist
            self.pdist_x = pdist_x
            self.pdist_y = pdist_y

class SAXSTimeSeriesDataset(TimeSeriesDataset, SAXSDataset):
    """
    Manages Small Angle X-ray Scattering time series datasets.
    """

    def __init__(self, infile, address="saxs", downsample=None,
        calc_mean=False, calc_error=True, error_mode="std", scale=False,
        **kwargs):
        """
        Initializes dataset.

        Arguments:
          infile (str): Path to input file, may contain environment
            variables
          usecols (list): Columns to select from DataFrame, once
            dataframe has already been loaded
          dt (float): Time interval between points; units unspecified
          toffset (float): Time offset to be added to all points (i.e.
            time of first point)
          downsample (int): Interval by which to downsample points
          verbose (int): Level of verbose output
          kwargs (dict): Additional keyword arguments

        .. todo:
          - Calculate error
          - Shift downsampling to superclass
        """
        from os.path import expandvars

        # Arguments
        verbose = kwargs.get("verbose", 1)

        # Load
        with h5py.File(expandvars(infile)) as h5_file:
            q = ["{0:5.3f}".format(a) for a in np.array(h5_file[address+"/q"])]
        super(SAXSTimeSeriesDataset, self).__init__(infile=infile,
          address=address+"/intensity", dataframe_kw=dict(columns=q), **kwargs)
        timeseries = self.timeseries = self.dataframe

        # Downsample
        if downsample:
            self.downsample(downsample, downsample_mode="mean", **kwargs)

        # Average over time series
        if calc_mean:
            self.dataframe = dataframe = pd.DataFrame(
              data=timeseries.mean(axis=0), columns=["intensity"])
            dataframe.index.name = "q"
            dataframe.index = np.array(timeseries.columns.values, np.float)

            # Scale
            if scale:
#                curve_fit_kw = dict(p0=(2e-9), bounds=(0.0,0.35))
                curve_fit_kw = dict(p0=(2e-9))  # Not clear why bounds broke
                curve_fit_kw.update(kwargs.get("curve_fit_kw", {}))
                scale = self.scale(scale, curve_fit_kw=curve_fit_kw, **kwargs)
                self.timeseries *= scale
        elif scale:
            self.timeseries *= scale
        if calc_error:
            se = self.calc_error(error_mode="block", **kwargs)
            se.name = "intensity_se"
            dataframe = self.dataframe = pd.concat([dataframe, se], axis=1)

class SAXSExperimentDataset(SAXSDataset):
    """
    Manages Small Angle X-ray Scattering experimental datasets.
    """

    def __init__(self, scale=False, **kwargs):
        """
        Initializes dataset.

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

class H5Dataset(object):
    """
    Class for managing hdf5 datasets

    .. todo:
      - Reimplement or depreciate
    """

    def __init__(self, **kwargs):
        """
        Initializes dataset.

        Arguments:
          infiles (list): List of infiles
          infile (str): Alternatively, single infile
        """
        self.default_address = kwargs.get("default_address", "")
        self.default_key     = kwargs.get("default_key",     "key")
        self.datasets        = {}
        self.attrs           = {}

        if   "infiles" in kwargs:
            self.load(infiles = kwargs.pop("infiles"))
        elif "infile"  in kwargs:
            self.load(infiles = [kwargs.pop("infile")])

    def load(self, infiles, **kwargs):
        """
        Loads data from h5 files.

        Arguments:
          infiles (list): infiles
        """
        from os.path import expandvars, isfile
        from h5py import File as h5
        import numpy as np
        import six

        for infile in infiles:
            if isinstance(infile, six.string_types):
                path    = expandvars(infile)
                address = self.default_address
                key     = self.default_key
            elif isinstance(infile, dict):
                path    = expandvars(infile.pop("path"))
                address = infile.pop("address", self.default_address)
                key     = infile.pop("key",     self.default_key)
            elif isinstance(infile, list):
                if len(infile) >= 1:
                    path = expandvars(infile[0])
                else:
                    raise OSError("Path to infile not provided")
                if len(infile) >= 2:
                    address = infile[1]
                else:
                    address = self.default_address
                if len(infile) >= 3:
                    key     = infile[2]
                else:
                    key     = self.default_key

            if not isfile(path):
                raise OSError("h5 file '{0}' does not exist".format(path))

            with h5(path) as in_h5:
                if address not in in_h5:
                    raise KeyError("Dataset {0}[{1}] not found".format(path,
                      address))
                dataset            = in_h5[address]
                self.datasets[key] = np.array(dataset)
                self.attrs[key]    = dict(dataset.attrs)
            print("Loaded Dataset {0}[{1}]; Stored at {2}".format(
              path, address, key))
