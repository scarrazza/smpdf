#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" SMPDF """

__author__ = 'Stefano Carrazza'
__license__ = 'GPL'
__version__ = '1.0.0'
__email__ = 'stefano.carrazza@mi.infn.it'

import os
import sys
import yaml
import argparse
import functools
from collections import defaultdict, OrderedDict


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.cbook import violin_stats
import scipy

from lhaindex import parse_info
from plotutils import ax_or_gca, violin_stats_from_dist
try:
    sys.path.append('applwrap')
    from applwrap import loadpdf, convolute
except ImportError:
    os.system("make -C applwrap")
    from applwrap import loadpdf, convolute

#TODO: Do we really want a subclass of dict?
class Config(dict):
    @classmethod
    def from_params(cls, **params):
        if 'pdfsets' not in params or not params['pdfsets']:
            raise ValueError("'pdfsets' not found in configuration.")
        #TODO make pdf a class
        for pdf in params['pdfsets']:
            #TODO: Do we allow incomplete sets at all? Seems like asking for
            #bugs.
            if not 'reps' in pdf or pdf['reps']=='all':
                pdf['reps'] = range(parse_info(pdf['name'])['NumMembers'])
            elif isinstance(pdf['reps'], int):
                pdf['reps'] = [pdf['reps']]
            elif isinstance(pdf['reps'], dict):
                pdf['reps'] = range(pdf['reps']['min'], pdf['reps']['max'])

        return cls(**params)

    @classmethod
    def from_yaml(cls, stream):
        return cls.from_params(**yaml.load(stream))

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

    def _violin_data(self):
        absdata = pd.concat(self.sample_values(1000),axis=1)
        reldata = absdata.as_matrix().T/self._cv.as_matrix()
        return reldata
    
    @ax_or_gca
    def violin_plot(self, ax=None):
        data = self._violin_data()
        ax.violinplot(data)


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
    
    def _violin_data(self):
        absdata = scipy.stats.norm(self.central_value.as_matrix(), 
                                   self.std_error())


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

    def _violin_data(self):
        return self._all_vals.as_matrix().T/self._cv.as_matrix()


RESULT_TYPES = defaultdict(lambda:Result,
                           symmhessian = SymHessianResult,
                           replicas   = MCResult,
                           )

def make_result(obs, pdf_name, datas):
    error_type = parse_info(pdf_name)['ErrorType']
    return RESULT_TYPES[error_type](obs, pdf_name, datas)




#TODO: implement in C++?
def convolve_all(pdf, observables):
    datas = defaultdict(lambda:OrderedDict())
    #TODO: load many replicas in C++
    #TODO: Could we loop over observables and then over memebers?
    for rep in pdf['reps']:
        for obs in observables:
            #TODO: hide this call from the api, do in convolute.
            loadpdf(pdf['name'], rep)
            res = convolute(obs['name'], obs['order'])
            datas[obs['name']][rep] = np.array(res)
    results = [make_result(obs, pdf['name'], datas[obs]) for obs in datas]
    return results


def main(conf):
    # perform convolution
    results = []
    for pdf in conf['pdfsets']:
        results += convolve_all(pdf, conf['observables'])

    #TODO: load many replicas in C++
    for result in results:
        for member in result.iterall():
            print "\n-", pdf['name'], "replica", member
            print "- APPLgrid convolution results:"
            for i, val in enumerate(result[member]):
                print ("\tData bin %i: %e" % (i, val))

    print ("\n +--------+ Completed +--------+\n")
    return results


def splash():
    print " "
    print "  ███████╗███╗   ███╗██████╗ ██████╗ ███████╗ "
    print "  ██╔════╝████╗ ████║██╔══██╗██╔══██╗██╔════╝ "
    print "  ███████╗██╔████╔██║██████╔╝██║  ██║█████╗   "
    print "  ╚════██║██║╚██╔╝██║██╔═══╝ ██║  ██║██╔══╝ "
    print "  ███████║██║ ╚═╝ ██║██║     ██████╔╝██║"
    print "  ╚══════╝╚═╝     ╚═╝╚═╝     ╚═════╝ ╚═╝"
    print "  __version__:", __version__
    print "  __authors__: S. Carrazza, Z. Kassabov\n"


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('config_yml', nargs='?',
                        help = "Path to the configuration file")

    args = parser.parse_args()
    if not (args.config_yml):
        parser.error("Too few arguments: config_yml")
    splash()
    # read yml file
    with open(args.config_yml,'r') as f:
        conf = Config.from_yaml(f)
    results = main(conf)
