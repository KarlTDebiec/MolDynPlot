#!/usr/bin/env python

from distutils.core import setup

setup(
    name             = "MolDynPlot",
    version          = "0.1",
    packages         = ["moldynplot",
                        "moldynplot.dataset",
                        "moldynplot.myplotspec",
                        "moldynplot.fpblockaverager"],
    license          = "BSD 3-clause",
    long_description = open("README.rst").read()
)
