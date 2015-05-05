# -*- coding: utf-8 -*-
"""
Created on Mon May  4 18:58:08 2015

@author: zah
"""
import os.path as  osp
import re
import inspect

import pandas as pd

import smpdflib as lib

#TODO: Specify base
def save_violins(results, output_dir, base_pdf=None):
    """Generate plots comparing the distributions obtained for the value of
    the observable using different PDF sets. If 'base_pdf' is specified, the
    values will be relative to the central value of that PDF."""
    for obs, fig in lib.compare_violins(results, base_pdf = base_pdf):
        filename = osp.join(output_dir, "figures", '%s.pdf' % obs)
        fig.savefig(filename)

def save_as(results, output_dir):
    """Generate plots showing the value of the observables as a function
    of a_s. The value is obtained by summing each bin in the applgrid."""
    if not isinstance(results, pd.DataFrame):
        results = lib.summed_results_table(results)
    for (process, nf), fig in lib.plot_alphaS(results):
        name = 'alpha_plot_%s(nf_%d).pdf'%(process,nf)
        name = re.sub(r'[^\w\.]', '_', name)
        fig.savefig(osp.join(output_dir, "figures", name))

#TODO: Implement missing methods
ACTION_DICT = {'violinplots':save_violins, 'asplots':save_as,
               'exporthtml': lambda : None,
               'exportcsv': lambda :None}


REALACTIONS = set(ACTION_DICT.keys())

METAACTION_DICT = {'all': (REALACTIONS, "Implies all other actions."),
                   'savedata': ({'exportcsv', 'exporthtml'}, "Export html and "
                                                                   "csv.")
                  }
METAACTIONS = set(METAACTION_DICT.keys())

SAVEDATA = {'exporthtml', 'exportcsv'}

ACTIONS = REALACTIONS | METAACTIONS

def parse_actions(acts):
    acts = set(acts)
    if acts - ACTIONS:
        raise ValueError("Unknown actions: %s" %  (acts - ACTIONS))
    realacts = acts & REALACTIONS
    for metaact in acts & METAACTIONS:
        realacts |= METAACTION_DICT[metaact][0]
    return realacts

def build_actions(acts, flt=None):
    if flt is None:
        return parse_actions(acts)
    return parse_actions(acts) & parse_actions(flt)

def do_actions(acts, **kwargs):
    for action in acts:
        func = ACTION_DICT[action]
        fargs = inspect.getargspec(func).args
        args = {k:v for k,v in kwargs.items() if k in fargs}
        yield action, func(**args)

def helptext(action):
    if action in ACTION_DICT:
        return ACTION_DICT[action].__doc__
    else:
        return METAACTION_DICT[action][1]

def gen_docs():
    s = "Possible actions are:\n%s" % '\n'.join('\t- %s : %s\n' %
        (action, helptext(action))
        for action in ACTION_DICT.keys())
    return s
