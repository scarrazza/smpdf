#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" SMPDF """
from __future__ import print_function
import sys
import os
import os.path as osp
import argparse
import shelve

import actions


__author__ = 'Stefano Carrazza'
__license__ = 'GPL'
__version__ = '1.0.0'
__email__ = 'stefano.carrazza@mi.infn.it'


def main(conf, output_dir, db, quiet=False):

    prefix_group = len(conf.actiongroups) > 1
    for group in conf.actiongroups:
        pdfsets, observables = group['pdfsets'], group['observables']
        # perform convolution
        #TODO: Use multicore
        results = lib.convolve_or_load(pdfsets, observables, db)
        data_table = lib.results_table(results)
        summed_table = lib.summed_results_table(results)
        import pandas as pd
        total = pd.concat((data_table,
                            summed_table),
                            ignore_index = True)
        if not quiet:
            print_results(results)
        kwargs = group.copy()
        kwargs.pop('actions')
        if not prefix_group:
            prefix = None
        else:
            prefix = group['prefix']
        kwargs.update({'results':results, 'output_dir':output_dir,
                       'prefix':prefix, 'data_table':data_table,
                       'summed_table':summed_table, 'total':total})
        for action, res in actions.do_actions(group['actions'], **kwargs):
            if not quiet:
                print("Finalized action '%s'." % action)

    #TODO: Think something better
    return total

def print_results(results):
    for result in results:
        for member in result.iterall():
            print ("\n- %s replica %d"% (result.pdf,member))
            print ("- APPLgrid convolution results:")
            for i, val in enumerate(result[member]):
                print ("\tData bin %i: %e" % (i, val))
    print ("\n +--------+ Completed +--------+\n")

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
    parser = argparse.ArgumentParser(
        description = "Compare phenomenology for arrays of applgrids "
         "and pdfsets. "
         "Fill yaml configuration files with the specification of the pdfsets, "
         "observables and actions (see examples for details).\n\n%s" %
         actions.gen_docs(),
       formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('config_yml',
                        help = "path to the configuration file")

    #TODO: Use db by default?
    parser.add_argument('--use-db', nargs='?', help="use Python database"
    " file of results and do not recompute those already in there"
    "If a file is not passed 'db/db' will be used", metavar='dbfile',
    const='db/db', default=None)

    parser.add_argument('-o','--output', help="output folder where to "
                                         "store resulting plots and tables",
                        default='output')

    parser.add_argument('-q','--quiet', help="Do not print the results",
                        action='store_true')

    args = parser.parse_args()

    splash()

    #Slow to import
    import matplotlib.pyplot as plt
    plt.style.use('main.mplstyle')

    import smpdflib as lib
    import config

    # read yml file
    try:
        with open(args.config_yml,'r') as f:
            conf = config.Config.from_yaml(f)
    except config.ConfigError as e:
        print("Bad configuration encountered:\n%s" % e.message,
              file=sys.stderr)
        sys.exit(1)
    except IOError as e:
        print("Cannot load configuration file:\n%s" % e.strerror,
              file=sys.stderr)
        sys.exit(1)
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
        results = main(conf, args.output ,db=db, quiet=args.quiet)
    finally:
        #bool(db) == False if empty
        if db is not None:
            db.close()
