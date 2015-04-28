#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" SMPDF """
import argparse
import smpdflib as lib

import matplotlib.pyplot as plt
plt.style.use('main.mplstyle')


__author__ = 'Stefano Carrazza'
__license__ = 'GPL'
__version__ = '1.0.0'
__email__ = 'stefano.carrazza@mi.infn.it'


def main(conf):
    # perform convolution
    pdfsets, observables = conf['pdfsets'], conf['observables']
    
    results = lib.convolve_all(pdfsets, observables)
    #TODO: load many replicas in C++
    for result in results:
        for member in result.iterall():
            print "\n-", result.pdf, "replica", member
            print "- APPLgrid convolution results:"
            for i, val in enumerate(result[member]):
                print ("\tData bin %i: %e" % (i, val))

    print ("\n +--------+ Completed +--------+\n")
    #TODO: Specify base
    for obs, fig in lib.compare_violins(results, base_pdf = results[0].pdf):
        fig.savefig("%s.pdf" % obs)
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
        conf = lib.Config.from_yaml(f)
    results = main(conf)
