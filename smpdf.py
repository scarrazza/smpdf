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

try:
    sys.path.append('applwrap')
    from applwrap import *
except ImportError:
    os.system("make -C applwrap")
    from applwrap import *

def main(config_yml):

    # read yml file
    conf = yaml.load(open(config_yml,'r'))

    # perform convolution
    for pdf in conf['pdfsets']:
        loadpdf(pdf['name'], pdf['reps'])
        for obs in conf['observables']:
            res = convolute(obs['name'], obs['order'])
            print "\n-", pdf['name'], "replica", pdf['reps']
            print "- APPLgrid convolution results:"
            for i in range(len(res)): print ("\tData bin %i: %e" % (i, res[i]))

    print "\n +--------+ Completed +--------+\n"


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


argnames = {'config_yml'}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('config_yml', nargs='?',
                        help = "Path to the configuration file")

    args = parser.parse_args()
    if not (args.config_yml):
        parser.error("Too few arguments: config_yml")
    mainargs = vars(args)
    splash()
    main(**mainargs)
