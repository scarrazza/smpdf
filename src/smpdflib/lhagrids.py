#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" A Monte Carlo to Hessian conversion tool """

__author__ = 'Stefano Carrazza'
__license__ = 'GPL'
__version__ = '1.0.0'
__email__ = 'stefano.carrazza@mi.infn.it'

import numpy as np
import fastcache

#TODO: Drop this dependence in favour of a custop py-3 compatible wrapper
#      that provides mkpdf.
import lhapdf


def make_xgrid(xminlog=1e-5, xminlin=1e-1, xmax=1, nplog=50, nplin=50):
    """Provides the points in x to sample the PDF. `logspace` and `linspace`
    will be called with the respsctive parameters."""

    return np.append(np.logspace(np.log10(xminlog), np.log10(xminlin),
                                       num=nplog, endpoint=False),
                     np.linspace(xminlin, xmax, num=nplin, endpoint=False)
                    )

def make_flavors(nf=3):
    return np.arange(-nf,nf+1)

#This is some 5 Mb for a 1000 replica set
@fastcache.lru_cache(maxsize=128, unhashable='ignore')
def get_pdf_values(pdf, Q, xgrid=None, fl=None):
    pdf_lha = load_lhapdf(pdf)
    if xgrid is None:
        xgrid = make_xgrid()
    #Allow tuples that can be saved in cache
    elif isinstance(xgrid, tuple):
        xgrid = make_xgrid(*xgrid)

    if fl is None:
        fl = make_flavors()
    elif isinstance(fl, int):
        fl = make_flavors(fl)
    elif isinstance(fl, tuple):
        fl = make_flavors(*fl)

    #TODO: Can we implement this loop in C
    all_members = [[[r.xfxQ(f, x, Q)
                     for x in xgrid]
                     for f in fl]
                     for r in pdf_lha]
    all_members = np.array(all_members)
    mean = all_members[0]
    replicas = all_members[1:]


    return mean, replicas

#TODO: Determine how much to cache this
@fastcache.lru_cache(maxsize=4)
def load_lhapdf(pdf):
    _pdf = lhapdf.mkPDFs(str(pdf))
    return _pdf
