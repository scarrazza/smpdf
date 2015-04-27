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
from collections import defaultdict

from lhaindex import parse_info
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

class Result():
    def __init__(self, obs, pdf, data):
        self.obs = obs
        self.pdf = pdf
        self._data = data

    @property
    def central_value(self):
        try:
           return self._data[0]
        except KeyError as e: #analysis:ignore
            #In python 3 it would be raise ValueError from e
            raise ValueError("No results for central value (Member 0) "
                             "provided")

    def hessian_error(nsigma=1):
        pass

    def __getitem__(self, item):
        return self._data[item]

    #TODO: Should this be the default iterator?
    def iterall(self):
        return iter(self._data)



#TODO: implement in C++?
def convolve_all(pdf, observables):
    datas = defaultdict(lambda:{})
    #TODO: load many replicas in C++
    #TODO: Could we loop over observables and then over memebers?
    for rep in pdf['reps']:
        for obs in observables:
            #TODO: hide this call from the api, do in convolute.
            loadpdf(pdf['name'], rep)
            res = convolute(obs['name'], obs['order'])
            datas[obs['name']][rep] = res
    results = [Result(obs, pdf['name'], datas[obs]) for obs in datas]
    return results


def main(config_yml):

    # read yml file
    with open(config_yml,'r') as f:
        conf = Config.from_yaml(f)

    # perform convolution
    for pdf in conf['pdfsets']:
        results = convolve_all(pdf, conf['observables'])
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
    mainargs = vars(args)
    splash()
    results = main(**mainargs)
