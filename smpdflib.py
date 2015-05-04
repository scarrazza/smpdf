# -*- coding: utf-8 -*-
from __future__ import division

"""
Created on Tue Apr 28 12:13:00 2015

@author: zah
"""
""" SMPDF """

__author__ = 'Stefano Carrazza'
__license__ = 'GPL'
__version__ = '1.0.0'
__email__ = 'stefano.carrazza@mi.infn.it'

import os
import os.path as osp
import sys
import functools
import glob
from collections import defaultdict, OrderedDict

import yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats

import lhaindex
import plotutils
try:
    sys.path.append('applwrap')
    from applwrap import initpdf, initobs, pdfreplica, convolute
except ImportError:
    os.system("make -C applwrap")
    from applwrap import initpdf, initobs, pdfreplica, convolute

ORDERS_QCD = {0: 'LO', 1: 'NLO', 2: 'NNLO'}

#for N_f = 4, LHAPDF's M_Z is actually M_{charm}
M_REF = defaultdict(lambda: 'Z', {4:'c'})

#TODO: Do we really want a subclass of dict?
class Config(dict):
    @classmethod
    def from_params(cls, **params):
        if 'pdfsets' not in params or not params['pdfsets']:
            raise ValueError("'pdfsets' not found in configuration.")

        pdfsets =  []
        for pdf in params['pdfsets']:
            pdfsets += [PDF(pdf['name'])]
        params['pdfsets'] = pdfsets


        observables = []
        for obs in params['observables']:
             names = glob.glob(obs['name'])
             if not names:
                 raise ValueError("No observables found for %s" % obs['name'])
             for name in names:
                 observables.append(Observable(name, obs['order']))
        params['observables'] = observables
        return cls(**params)

    @classmethod
    def from_yaml(cls, stream):
        return cls.from_params(**yaml.load(stream))


class TupleComp(object):
    def __hash__(self):
        return hash(self.get_key())

    def __eq__(self, other):
        return (isinstance(other, self.__class__)
                and self.get_key() == other.get_key())

class Observable(TupleComp):
    def __init__(self, name, order):
        self.filename = name
        self.order = order

    @property
    def name(self):
        return osp.splitext(osp.basename(self.filename))[0]

    def __str__(self):
        return "%s(%s)"%(self.name, ORDERS_QCD[self.order])

    def __repr__(self):
        return "<%s:%s>" % (self.__class__.__name__, self.__str__())

    def get_key(self):
        return (self.name, self.order)


class PDF(TupleComp):
    def __init__(self, name):
        self.name = name

    def get_key(self):
        return (self.name,)

    def __str__(self):
        return self.name

    def __repr__(self):
        return "<%s:%s>" % (self.__class__.__name__, self.name)

    @property
    def oqcd_str(self):
        return ORDERS_QCD[self.OrderQCD]

    @property
    def mref(self):
        return "M_%s" % M_REF[self.NumFlavors]

    @property
    def reps(self):
        return range(self.NumMembers)

    def __getattr__(self, name):
        return lhaindex.parse_info(self.name)[name]


#TODO: Decide if we really want this
def _check_central(f):
    @functools.wraps(f)
    def _f(self, *args, **kwargs):
        if not 0 in self._data:
            raise ValueError("No results for central value (Member 0) "
                             "provided")
        return f(self, *args, **kwargs)
    return _f

class Result():
    def __init__(self, obs, pdf, data):
        self.obs = obs
        self.pdf = pdf
        self._data = pd.DataFrame(data)

    @property
    @_check_central
    def central_value(self):
        return self._cv

    @property
    def _cv(self):
        return self._data[0]

    @property
    def _all_vals(self):
        return self._data.iloc[:,1:]

    @property
    def nrep(self):
        return self._all_vals.shape[1]

    def std_error(self, nsigma=1):
        raise NotImplementedError("No error computation implemented for this"
                                  "type of set")

    def sample_values(self, n):
        raise NotImplementedError("No sampling implemented for this"
                                  "type of set")

    @_check_central
    def std_interval(self, nsigma=1):
        std = self.std_error(nsigma)
        return pd.DataFrame({'min':self._cv - std,
                             'max':self._cv + std})

    def __getitem__(self, item):
        return self._data[item]

    #TODO: Should this be the default iterator?
    def iterall(self):
        return iter(self._data)

    def _violin_data(self, rel_to=None):
        absdata = pd.concat(self.sample_values(10000),axis=1)
        if rel_to is None:
            rel_to = 1
        reldata = absdata.as_matrix().T/rel_to
        return reldata

    def violin_plot(self, data=None , **kwargs):
        if data is None:
            data = self._violin_data()

        myargs = {'label': str(self.pdf)}
        myargs.update(kwargs)
        return plotutils.violin_plot(data, **myargs)



class SymHessianResult(Result):

    @_check_central
    def std_error(self, nsigma=1):
        diffsq = (self._all_vals.subtract(self._cv, axis=0))**2
        return diffsq.sum(axis=1).apply(np.sqrt)*nsigma

    def sample_values(self, n):
        diffs = self._all_vals.subtract(self._cv, axis=0)
        for _ in range(n):
            weights = np.random.normal(size=self.nrep)
            error = (diffs*weights).sum(axis=1)
            yield self._cv + error

    def _violin_data(self, rel_to = None):
        std = self.std_error()
        mean = self.central_value.as_matrix()

        if rel_to is None:
            rel_to = np.ones_like(mean)
        vpstats = []
        for m,s,r in zip(mean,std, rel_to):
            # Dictionary of results for this distribution
            stats = {}

            # Calculate basic stats for the distribution
            min_val = m - 3*s
            max_val = m + 3*s
            coords = np.linspace(m - 3*s, m+3*s, 1000)
            # Evaluate the kernel density estimate
            stats['vals'] = scipy.stats.norm(m,s).pdf(coords)*r
            stats['coords'] = coords/r

            # Store additional statistics for this distribution
            stats['mean'] = m/r
            stats['median'] = m/r
            stats['min'] = min_val/r
            stats['max'] = max_val/r

            # Append to output
            vpstats.append(stats)

        return vpstats




class MCResult(Result):
    #TODO: Is it correct to consider each bin as independant here?
    @_check_central
    def centered_interval(self, percent=68):
        n = percent*self.nrep//100
        def get_lims(row):
            s = np.argsort(np.abs(row))
            sel = row[s][:n]
            return pd.Series({'min':np.min(sel), 'max':np.max(sel)})

        diffs = self._all_vals.subtract(self._cv, axis=0)
        return diffs.apply(get_lims, axis=1).add(self._cv, axis=0)

    @property
    @_check_central
    def std_error(self, nsigma=1):
        return self._all_vals.std(axis=1)*nsigma

    def sample_values(self, n):
        for _ in range(n):
            col = np.random.choice(self._all_vals.columns)
            yield self._all_vals[col]

    def _violin_data(self, rel_to = None):
        if rel_to is None:
            rel_to = 1
        return self._all_vals.as_matrix().T/ rel_to

def aggregate_results(results):
    combined = defaultdict(lambda: {})
    for result in results:
        combined[result.obs][result.pdf] = result
    return combined

def compare_violins(results, base_pdf = None):
    if not isinstance(results, dict):
        combined = aggregate_results(results)
    else:
        combined = results
    for obs in combined:
        figure = plt.figure()
        norms = None
        handles = []
        plt.title(str(obs))
        ncolors = len(combined[obs])
        colors = iter(plotutils.get_accent_colors(ncolors).mpl_colors)
        alpha = 1
        base = combined[obs].get(base_pdf, None)
        results = sorted(combined[obs].values(), key = lambda x: x!=base)
        for result in results:

            if base is not None:
                cv = base.central_value.as_matrix()
                data = result._violin_data(rel_to=cv)
            else:
                data = data = result._violin_data()
            color = next(colors) + (alpha,)
            alpha /= 2
            plot, handle, norms = result.violin_plot(data, color=color,
                                              showextrema=False,
                                              normvalues=norms)
            handles.append(handle)
        plt.xlabel('bins')
        plt.ylabel('Rel to %s' % base_pdf)
        plt.xticks(range(1,len(result.central_value) + 1))
        plt.legend(handles=handles, loc='best')
        yield obs, figure


RESULT_TYPES = defaultdict(lambda:Result,
                           symmhessian = SymHessianResult,
                           replicas   = MCResult,
                           )

def make_result(obs, pdf, datas):
    error_type = pdf.ErrorType
    return RESULT_TYPES[error_type](obs, pdf, datas)


def make_convolution(pdf, observables):
    datas = defaultdict(lambda:OrderedDict())
    #TODO: load many replicas in C++
    #TODO: Could we loop over observables and then over memebers?
    initpdf(pdf.name)
    for obs in observables:
        initobs(obs.filename)
        for rep in pdf.reps:
            #TODO: hide this call from the api, do in convolute.
            sys.stdout.write('\r-> Computing replica %d of %s' %
                             (rep, pdf))
            sys.stdout.flush()
            pdfreplica(rep)
            res = convolute(obs.order)
            datas[obs][rep] = np.array(res)
        sys.stdout.write('\n')
    return datas

def results_from_datas(dataset):
    results = []
    for pdf in dataset:
        data = dataset[pdf]
        results += [make_result(obs, pdf, data[obs]) for obs in data]
    return results

#TODO: Refactor this after adding efficient convolution
def get_dataset(pdfsets, observables, db=None):
    def make_key(pdf, obs):
        return str((pdf.get_key(), obs.get_key()))
    dataset = OrderedDict()
    for pdf in pdfsets:
        #bool(db) == False if empty
        if db is not None:
            res = {}
            obs_to_compute = []
            for obs in observables:
                key = make_key(pdf, obs)
                if key in db:
                    res[obs] = db[key]
                else:
                    obs_to_compute.append(obs)

            computed_data = make_convolution(pdf, obs_to_compute)
            for newobs in computed_data:
                key = make_key(pdf, newobs)
                db[key] = computed_data[newobs]
            res.update(computed_data)


        else:
            res = make_convolution(pdf, observables)

        dataset[pdf] = res
    return dataset

def convolve_or_load(pdfsets, observables, db=None):
    #results = []
    results = results_from_datas(get_dataset(pdfsets, observables, db))
    return results
