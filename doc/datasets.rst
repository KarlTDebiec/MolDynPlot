Datasets
========
Moldynplot includes several dataset classes that build on
:class:`Dataset<myplotspec:myplotspec.Dataset.Dataset>` will additions
specific for molecular dynamics simulation data.

SequenceDataset
---------------
.. autoclass::  moldynplot.Dataset.SequenceDataset

IREDSequenceDataset
___________________
.. autoclass::  moldynplot.Dataset.IREDSequenceDataset
.. automethod:: moldynplot.Dataset.IREDSequenceDataset.construct_argparser
.. automethod:: moldynplot.Dataset.IREDSequenceDataset.identify_infile
.. automethod:: moldynplot.Dataset.IREDSequenceDataset.parse_ired
.. automethod:: moldynplot.Dataset.IREDSequenceDataset.average_independent

TimeSeriesDataset
-----------------
.. autoclass::  moldynplot.Dataset.TimeSeriesDataset
.. automethod:: moldynplot.Dataset.TimeSeriesDataset.downsample
.. automethod:: moldynplot.Dataset.TimeSeriesDataset.calc_pdist
.. automethod:: moldynplot.Dataset.TimeSeriesDataset.calc_error

IREDTimeSeriesDataset
_____________________
.. autoclass::  moldynplot.Dataset.IREDTimeSeriesDataset
.. automethod:: moldynplot.Dataset.IREDTimeSeriesDataset.construct_argparser

NatConTimeSeriesDataset
_______________________
.. autoclass::  moldynplot.Dataset.NatConTimeSeriesDataset

SAXSTimeSeriesDataset
_____________________
.. autoclass::  moldynplot.Dataset.SAXSTimeSeriesDataset

CorrDataset
-----------
.. autoclass::  moldynplot.Dataset.CorrDataset
.. automethod:: moldynplot.Dataset.CorrDataset.get_cache_key
.. automethod:: moldynplot.Dataset.CorrDataset.get_cache_message

SAXSDataset
-----------
.. autoclass::  moldynplot.Dataset.SAXSDataset
.. automethod:: moldynplot.Dataset.SAXSDataset.scale

SAXSExperimentDataset
_____________________
.. autoclass::  moldynplot.Dataset.SAXSExperimentDataset

SAXSTimeSeriesDataset
_____________________
.. autoclass::  moldynplot.Dataset.SAXSTimeSeriesDataset

SAXSDiffDataset
_______________
.. autoclass::  moldynplot.Dataset.SAXSDiffDataset

MDGXDataset
-----------
.. autoclass::  moldynplot.Dataset.MDGXDataset
.. automethod:: moldynplot.Dataset.MDGXDataset.get_cache_key

H5Dataset
---------
.. autoclass::  moldynplot.Dataset.H5Dataset
.. automethod:: moldynplot.Dataset.H5Dataset.load
