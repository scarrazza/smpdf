# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 11:32:30 2015

@author: zah
"""
import functools
import itertools

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cbook import violin_stats
import matplotlib.patches
import matplotlib.mlab as mlab
import matplotlib.colors as colors



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

def extend_range(rmin, rmax, prop=0.1):
    diff = (rmax- rmin)*prop
    if diff <= 0:
        return (None, None)
    return rmin - diff, rmax + diff

def get_accent_colors(num_colors):
    import palettable
    num_colors = 3 if num_colors < 3 else 8 if num_colors > 8 else num_colors
    cmap = palettable.colorbrewer.get_map('Accent', 'qualitative', num_colors)
    return cmap

def get_spectral_colors(num_colors):
    import palettable
    num_colors = 3 if num_colors < 3 else 8 if num_colors > 11 else num_colors
    cmap = palettable.colorbrewer.get_map('Spectral', 'diverging', num_colors)
    return cmap

def get_set1_colors(num_colors):
    import palettable
    num_colors = 3 if num_colors < 3 else 9 if num_colors > 9 else num_colors
    cmap = palettable.colorbrewer.get_map('Set1', 'qualitative', num_colors)
    return cmap

def color_names_to_rgb(colorlist):
    for color in itertools.cycle(colorlist):
        yield colors.colorConverter.to_rgb(color)

@ax_or_gca
def mpl_default_color_tuple(ax=None):
    for color in color_names_to_rgb(ax._get_lines.color_cycle):
        yield color

@ax_or_gca
def violin_plot(data, normvalues=None, ax=None, bw_method=None, **kwargs):

    def _kde_method(X, coords):
            kde = mlab.GaussianKDE(X, bw_method)
            return kde.evaluate(coords)

    myargs = {}
    myargs.update(kwargs)
    if 'color' in myargs:
        color = myargs.pop('color')
    else:
        color = None
    if 'label' in myargs:
        label = myargs.pop('label')
    else:
        label = None

    if 'hatches' in myargs:
        hatches = myargs.pop('hatches')
    else:
        hatches = None

    if isinstance(data, list):
        stats = data
    else:
        stats = violin_stats(data, _kde_method)

    N = len(stats)

    if normvalues is not None:
        if np.isscalar(normvalues):
            normvalues = [normvalues] * N
        elif len(normvalues) != N:
            raise ValueError("Incorrect number of normvalues")

        widths = [normval*np.max(stat['vals']) for normval, stat
                  in zip(normvalues, stats)]
        myargs['widths'] = widths

    if 'widths' in myargs:
        widths = myargs['widths']
        if np.isscalar(widths):
            widths = [widths] * N
        elif len(widths) != N:
            raise ValueError("Incorrect number of widths")
        myargs['widths'] = widths
    else:
        myargs['widths'] = [0.5]*N



    ournorms = [w/np.max(stat['vals']) for w,stat in zip(myargs['widths'],
               stats)]

    vp = ax.violin(stats, **myargs)
    vp_edge = ax.violin(stats, **myargs)

    for pc , edge in zip(vp['bodies'], vp_edge['bodies']):
        if color:

            if len(color) == 4:
                pc.set_alpha(color[3])
                edge.set_alpha(1)
            pc.set_facecolor(color)
            pc.set_edgecolor('none')
            edge.set_edgecolor(color[:3])
            edge.set_facecolor('none')
        if hatches:
            pc.set_hatch(hatches)
    if label:
        if not color:
            color =  vp['bodies'][0].get_facecolor()[0]
        vp['bodies'][0].set_label(label)
        handle = matplotlib.patches.Patch(facecolor=color, label=label,
                                          hatch=hatches, edgecolor=color[:3])
    else:
        handle = None

    return vp, handle, ournorms
