# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 09:46:35 2015

@author: zah
"""
import functools
import os.path as osp


_initialized = set()

def _initializer(f):
    @functools.wraps(f)
    def _f(*args, **kwargs):
        if f in _initialized:
            raise RuntimeError("Initializer can only be called once: %s" %
                               f.__qualname__)
        _initialized.add(f)
        return f(*args, **kwargs)
    return _f

#@_initializer
def init_style():
    """Set style to that of our app"""
    import matplotlib.pyplot as plt
    libpath = osp.dirname(__file__)
    stylefilename = osp.join(libpath, 'small.mplstyle')
    plt.style.use(stylefilename)

@_initializer
def init_backend():
    """The Qt backend seems to cause conflicts with LHAPDF (grids appear
    corrupted for some reason). Thats why we set the 'do nothing' Agg
    backend. This is also useful in servers without X system, where it would.
    give errors otherwise"""
    import matplotlib
    matplotlib.use('Agg')

@_initializer
def init_multiprocessong():
    """Set multiprcesing method to 'spawn'. All memory is copied, and bugs
    with APPLgrid are when loading multiple grids are avoided for good.
    """
    import multiprocessing
    multiprocessing.set_start_method('spawn')

@_initializer
def init_app():
    """Global initialization. This sets the
    `spawn` method for multiprocessing and the Agg backend for matplotlib
    which are needed for APPlgrid and LHAPDF to work poperly.
    Also set the matplotlib style to
    those of SMPDF. **Must** be called **before** importing `matplotlib`.
    """
    init_multiprocessong()
    init_backend()
    init_style()


