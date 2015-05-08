# -*- coding: utf-8 -*-
"""
Created on Mon May  4 18:58:08 2015

@author: zah
"""
import os.path as  osp
import re
from collections import OrderedDict

import pandas as pd

#http://stackoverflow.com/questions/26277757/pandas-to-html-truncates-string-contents
pd.set_option('display.max_colwidth', -1)


#TODO: Specify base
def save_violins(results, output_dir, base_pdf=None, prefix=None):
    """Generate plots comparing the distributions obtained for the value of
    the observable using different PDF sets. If 'base_pdf' is specified, the
    values will be relative to the central value of that PDF."""

    #slow to import
    import smpdflib.core as lib
    for obs, fig in lib.compare_violins(results, base_pdf = base_pdf):
        filename = "%s%s.pdf" % (prefix if prefix else '', obs)
        path = osp.join(output_dir, "figures", filename)
        fig.savefig(path)

def save_as(summed_table, output_dir, prefix = None):
    """Generate plots showing the value of the observables as a function
    of a_s. The value is obtained by summing each bin in the applgrid."""

    #slow to import
    import smpdflib.core as lib
    for (process, nf), fig in lib.plot_alphaS(summed_table):
        name = '%salpha_plot_%s(nf_%d).pdf'%(prefix if prefix else '',
                                             process,nf)
        name = re.sub(r'[^\w\.]', '_', name)
        fig.savefig(osp.join(output_dir, "figures", name))

#TODO: Refactor this so there is not so much back and forth with smpdflib
def export_html(total, output_dir, prefix = None):
    """Export results as a rich HTML table."""
    import smpdflib.core as lib
    filename = "%sresults.html" % (prefix if prefix else '')
    lib.save_html(total[lib.DISPLAY_COLUMNS], osp.join(output_dir, filename))

#TODO: Ability to import exported csv
def export_csv(total, output_dir, prefix = None):
    """Export results as a CSV so they can be processed by other tools. The
    resulting file is tab-separated."""
    filename = "%sresults.csv" % (prefix if prefix else '')
    total.to_csv(osp.join(output_dir, filename), sep='\t')

#TODO: Think how to make this better
def test_as_linearity(summed_table, diff_from_line = 0.25):
    """Test linearity of value of the value of the observable as a function of
    as. If th test fails for some point, report it in the tables and in
    subsequent plots."""
    import smpdflib.core as lib
    return lib.test_as_linearity(summed_table, diff_from_line = diff_from_line)



ACTION_DICT = OrderedDict((
               ('testas',  test_as_linearity),
               ('violinplots',save_violins),
               ('asplots',save_as),
               ('exporthtml', export_html),
               ('exportcsv', export_csv),
               ))


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
    #inspect is slow to import, and noticeable in --help.
    import inspect

    for action in ACTION_DICT.keys():
        if not action in acts:
            continue
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
        for action in list(ACTION_DICT.keys()) + list(METAACTION_DICT.keys()))
    return s
