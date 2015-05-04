#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" SMPDF """
from __future__ import print_function
import sys
import os
import os.path as osp
import argparse
import shelve

import matplotlib.pyplot as plt
plt.style.use('main.mplstyle')

import smpdflib as lib
import actions


__author__ = 'Stefano Carrazza'
__license__ = 'GPL'
__version__ = '1.0.0'
__email__ = 'stefano.carrazza@mi.infn.it'


def main(conf, output_dir, db):
    # perform convolution
    pdfsets, observables = conf['pdfsets'], conf['observables']
    #TODO: Use multicore
    results = lib.convolve_or_load(pdfsets, observables, db)
    #TODO: load many replicas in C++
    for result in results:
        for member in result.iterall():
            print ("\n- %s replica %d"% (result.pdf,member))
            print ("- APPLgrid convolution results:")
            for i, val in enumerate(result[member]):
                print ("\tData bin %i: %e" % (i, val))
    print ("\n +--------+ Completed +--------+\n")
    actions.save_violins(results, base_pdf=results[0], output_dir=output_dir)
    actions.save_as(results, output_dir=output_dir)

    return results

def make_output_dir(output_dir):
    if not osp.exists(args.output):
        os.makedirs(args.output)
    elif not osp.isdir(args.output):
        print("'output' is not a directory", file=sys.stderr)
        sys.exit(1)
    figpath = osp.join(output_dir,"figures")
    if not osp.isdir(figpath):
        os.mkdir(figpath)


def splash():
    s =  ("""
  ███████╗███╗   ███╗██████╗ ██████╗ ███████╗
  ██╔════╝████╗ ████║██╔══██╗██╔══██╗██╔════╝
  ███████╗██╔████╔██║██████╔╝██║  ██║█████╗
  ╚════██║██║╚██╔╝██║██╔═══╝ ██║  ██║██╔══╝
  ███████║██║ ╚═╝ ██║██║     ██████╔╝██║
  ╚══════╝╚═╝     ╚═╝╚═╝     ╚═════╝ ╚═╝
  __version__: %s"
  __authors__: S. Carrazza, Z. Kassabov
""")%__version__
    print(s)




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('config_yml',
                        help = "Path to the configuration file")

    #TODO: Use db by default?
    parser.add_argument('--use-db', nargs='?', help="Use Python database"
    " file of results and do not recompute those already in there. "
    "If a file is not passed 'db/db' will be used", metavar='dbfile',
    const='db/db', default=None)

    parser.add_argument('--output', help="Output folder where to "
                                         "store resulting plots and tables.",
                        default='output')

    args = parser.parse_args()

    splash()
    # read yml file
    with open(args.config_yml,'r') as f:
        conf = lib.Config.from_yaml(f)
    make_output_dir(args.output)
    #TODO: handle this better

    dbfolder = args.use_db
    if dbfolder:
        dirname = osp.dirname(dbfolder)
        if dirname and not osp.isdir(dirname):
            os.makedirs(dirname)
        db = shelve.open(args.use_db)
    else:
        db = None
    try:
        results = main(conf, args.output ,db=db)
    finally:
        #bool(db) == False if empty
        if db is not None:
            db.close()
