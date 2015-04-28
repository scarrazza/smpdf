# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 11:32:30 2015

@author: zah
"""
import functools

import numpy as np
import matplotlib.pyplot as plt

import palettable


#TODO: use inspect to allow ax as arg?
def ax_or_gca(f):
    @functools.wraps(f)
    def _f(*args, **kwargs):
            if 'ax' not in kwargs or kwargs['ax'] is None:
                kwargs['ax'] = plt.gca()
            return f(*args, **kwargs)
    return _f

def ax_or_newfig(f):
    @functools.wraps(f)
    def _f(*args, **kwargs):
        noax = 'ax' not in kwargs or kwargs['ax'] is None
        if noax:
                plt.figure()
                kwargs['ax'] = plt.gca()
        result = f(*args, **kwargs)
        if noax:
            plt.legend(loc = 'best')
        return result
    
    return _f

def violin_stats_from_dist(coords, X):
    '''
    Like matplotlib.cbook.violin_stats, but assumes X is already a valid KDE
    at coords.

    Parameters
    ----------
    X : array-like
        Sample data that will be used to produce the gaussian kernel density
        estimates. Must have 2 or fewer dimensions.

    method : callable
        The method used to calculate the kernel density estimate for each
        column of data. When called via `method(v, coords)`, it should
        return a vector of the values of the KDE evaluated at the values
        specified in coords.

    points : scalar, default = 100
        Defines the number of points to evaluate each of the gaussian kernel
        density estimates at.

    Returns
    -------

    A list of dictionaries containing the results for each column of data.
    The dictionaries contain at least the following:

        - coords: A list of scalars containing the coordinates this particular
          kernel density estimate was evaluated at.
        - vals: A list of scalars containing the values of the kernel density
          estimate at each of the coordinates given in `coords`.
        - mean: The mean value for this column of data.
        - median: The median value for this column of data.
        - min: The minimum value for this column of data.
        - max: The maximum value for this column of data.
    '''

    # List of dictionaries describing each of the violins.
    vpstats = []

    # Want X to be a list of data sequences
    X = np.atleast_2d(X)

    for x in X:
        # Dictionary of results for this distribution
        stats = {}

        # Calculate basic stats for the distribution
        min_val = np.min(x)
        max_val = np.max(x)

        # Evaluate the kernel density estimate
        stats['vals'] = x
        stats['coords'] = coords

        # Store additional statistics for this distribution
        stats['mean'] = np.mean(x)
        stats['median'] = np.median(x)
        stats['min'] = min_val
        stats['max'] = max_val

        # Append to output
        vpstats.append(stats)

    return vpstats
    

def center_cmap_args(data,cmap=None):
    M = np.max(np.abs(data))
    return dict(c =(data),
                   cmap = cmap, vmin = -M, vmax = M)

def get_accent_colors(num_colors):
    num_colors = 3 if num_colors < 3 else 8 if num_colors > 8 else num_colors
    return palettable.colorbrewer.get_map('Accent', 'qualitative', num_colors)