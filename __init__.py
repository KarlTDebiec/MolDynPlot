#!/usr/bin/python
# -*- coding: utf-8 -*-
#   myplotspec_sim.__init__.py
#
#   Copyright (C) 2015 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
"""
################################### MODULES ###################################
from __future__ import absolute_import,division,print_function,unicode_literals
################################## FUNCTIONS ##################################
def hsl_to_rgb(h, s, l, **kwargs):
    """
    **Arguments:**
        :*h*:   hue          (0.0 - 1.0)
        :*s*:   saturation   (0.0 - 1.0)
        :*l*:   luminescence (0.0 - 1.0)

    **Returns:**
        :*rgb*: Numpy array of equivalent red, green, blue

    .. todo::
        - Should smoothly accept separate h, s, l variables or list or
          numpy array of length 3
    """
    c   = (1. - abs(2. * l - 1.)) * s
    hp  = h * 6
    i   = c * (1. - abs(hp % 2. - 1.))
    m   = l - 0.5 * c
    if   hp >= 0. and hp <  1.:
        return np.array([c,  i,  0.]) + m
    elif hp >= 1. and hp <  2.:
        return np.array([i,  c,  0.]) + m
    elif hp >= 2. and hp <  3.:
        return np.array([0., c,   i]) + m
    elif hp >= 3. and hp <  4.:
        return np.array([0., i,   c]) + m
    elif hp >= 4. and hp <  5.:
        return np.array([i,  0.,  c]) + m
    elif hp >= 5. and hp <= 6.:
        return np.array([c,  0.,  i]) + m

def rgb_to_hsl(r, g, b):
    """
    **Arguments:**
        :*r*:   red   (0.0 - 1.0)
        :*g*:   green (0.0 - 1.0)
        :*b*:   blue  (0.0 - 1.0)

    **Returns:**
        :*hsl*: Numpy array of equivalent hue, saturation, luminescence

    .. todo::
        - Should smoothly accept separate r, g, b variables or list or
          numpy array of length 3
    """
    x   = max(r, g, b)
    n   = min(r, g, b)
    c   = x - n
    l   = 0.5 * (x + n)
    if   c == 0.:
        h = 0.
    elif x == r:
        h = (((g - b) / c) % 6.) / 6.
    elif x == g:
        h = (((b - r) / c) + 2.) / 6.
    elif x == b:
        h = (((r - g) / c) + 4.) / 6.
    if   c == 0:
        s = 0.
    else:
        s = c / (1. - np.abs(2. * l - 1.))
    return np.array([h, s, l])

def concentration(n, volume):
    """
    """
    return float(n) / 6.0221415e23 / volume
def P_bound_to_Ka(P_bound, C_mol1_total, C_mol2_total):
    """
    """
    C_complex      = P_bound * np.min(C_mol1_total, C_mol2_total)
    C_mol1_unbound = C_mol1_total - C_complex
    C_mol2_unbound = C_mol2_total - C_complex
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return C_complex / (C_mol1_unbound * C_mol2_unbound)
def Pbound_se_to_KA_se(P_bound, C_mol1_total, C_mol2_total, P_bound_se):
    """
    """
    return np.sqrt(
      (((C_mol2_total - C_mol1_total * P_bound ** 2) * P_bound_se) /
      ((P_bound - 1) ** 2 *
      (C_mol1_total * P_bound - C_mol2_total) ** 2)) ** 2)
def KA_to_Pbound(Ka, C_mol1_total, C_mol2_total):
    """
    """
    from scipy.optimize import fmin

    def model_function(P_bound, Ka, C_mol1_total, C_mol2_total):
        C_min = np.min(C_mol1_total, C_mol2_total)
        return (((P_bound * C_min) /
          ((C_mol1_total - P_bound * C_min) *
          (C_mol2_total - P_bound * C_min))) - Ka) ** 2
    P_bound   = fmin(func = model_function, x0 = 0.1,
      args = (Ka, C_mol1_total, C_mol2_total), disp = False)
    return P_bound
