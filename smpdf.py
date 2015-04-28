#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" SMPDF """
import os
import os.path as osp
import argparse
import shelve

import matplotlib.pyplot as plt
plt.style.use('main.mplstyle')

import smpdflib as lib


__author__ = 'Stefano Carrazza'
__license__ = 'GPL'
__version__ = '1.0.0'
__email__ = 'stefano.carrazza@mi.infn.it'


def main(conf, db):
    # perform convolution
    pdfsets, observables = conf['pdfsets'], conf['observables']
    results = lib.convolve_or_load(pdfsets, observables, db)
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
    parser.add_argument('config_yml',
                        help = "Path to the configuration file")

    #TODO: Use db by default?
    parser.add_argument('--use-db', nargs='?', help="Use Python database"
    " file of results and do not recompute those already in there. "
    "If a file is not passed 'db/db' will be used", metavar='dbfile',
    const='db/db', default=None)

    args = parser.parse_args()

    splash()
    # read yml file
    with open(args.config_yml,'r') as f:
        conf = lib.Config.from_yaml(f)
    #TODO: handle this better
    filename = args.use_db
    if filename:

        dirname = osp.dirname(filename)
        if dirname and not osp.isdir(dirname):
            os.makedirs(dirname)
        db = shelve.open(args.use_db)
    else:
        db = None
    try:
        results = main(conf, db=db)
    finally:
        #bool(db) == False if empty
        if db is not None:
            db.close()
