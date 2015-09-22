# -*- coding: utf-8 -*-
"""
Created on Mon May  4 18:58:08 2015

@author: zah
"""
#TODO: Call this 'actionlib' and move actual actions to another module.
import os
import os.path as  osp
import re
import shutil
from collections import OrderedDict
import textwrap
import inspect

class ActionError(Exception):
    pass

def appends_to(resource):
    def decorator(f):
        f.appends_to = resource
        return f
    return decorator

def provides(resource):
    def decorator(f):
        f.provides = resource
        return f
    return decorator

def check(check_func):
    def decorator(f):
        if not hasattr(f, 'checks'):
            f.checks = []
        f.checks += [check_func]
        return f
    return decorator

def require_args(*args):
    def check_func(action, group, config):
        for arg in args:
            if not arg in group:
                raise ActionError("Action %s requires a parameter '%s'" %
                                  (action, arg))
    def decorator(f):
        f.required_args = args
        return check(check_func)(f)
    return decorator

def normalize_name(name):
    return re.sub(r'[\.\,]', '', str(name))


def save_figures(generator, table, output_dir, namefunc=None,
                 prefix=None, fmt ='pdf', **kwargs):

        import matplotlib.pyplot as plt
        if isinstance(fmt, str):
            fmts = [fmt]
        else:
            fmts = fmt
        prefixstr = prefix if prefix else ''
        if namefunc is None:
            namefunc = lambda *spec: ''.join(str(x) for x in spec)
        for namespec, fig in generator(table, **kwargs):
            nameresult = namefunc(*namespec)
            for fmt in fmts:
                filename = "{prefixstr}{nameresult}.{fmt}".format(**locals())
                path = osp.join(output_dir, "figures", filename)
                fig.savefig(path, bbox_inches='tight')
                plt.close(fig)

def check_know_errors(action, group, config):
    from smpdflib.core import RESULT_TYPES, PDFDoesNotExist

    bad_types = []
    for pdf in group['pdfsets']:
        try:
            error_type = pdf.ErrorType
        except PDFDoesNotExist:
            continue
        if  error_type not in RESULT_TYPES:
            bad_types.append((pdf, error_type))
    if bad_types:
        msg = ("Action '%s' cannot be executed because the following PDF sets "
              "have unknown error types:") % action
        msg += ''.join(['\n %s with ErrorType: %s' % bad for bad in bad_types])
        raise ActionError(msg)

@check(check_know_errors)
def save_violins(results, output_dir, prefix, base_pdf=None, fmt='pdf'):
    """
    Generate plots comparing the distributions obtained for the value of
    the observable using different PDF sets. If 'base_pdf' is specified, the
    values will be relative to the central value of that PDF."""
    #slow to import
    import smpdflib.plots as plots
    def namefunc(obs):
        return "violinplot_%s"%normalize_name(obs)
    return save_figures(plots.compare_violins, results, output_dir,
                        base_pdf=base_pdf,
                        prefix=prefix, fmt=fmt, namefunc=namefunc)

@check(check_know_errors)
def save_cis(results, output_dir, prefix, base_pdf=None, fmt='pdf'):
    """
    Generate plots comparing the confidence intervals for the value of
    the observable using different PDF sets. If 'base_pdf' is specified, the
    values will be relative to the central value of that PDF."""
    #slow to import
    import smpdflib.plots as plots
    def namefunc(obs):
        return "ciplot_%s"%normalize_name(obs)
    return save_figures(plots.compare_cis, results, output_dir,
                        base_pdf=base_pdf,
                        prefix=prefix, fmt=fmt, namefunc=namefunc)
@check(check_know_errors)
def save_as(summed_table, output_dir, prefix, fmt='pdf'):
    """Generate plots showing the value of the observables as a function
    of a_s. The value is obtained by summing each bin in the applgrid."""

    #slow to import
    import smpdflib.plots as plots
    def namefunc(obs, nf, bin_):
        obs = normalize_name(obs)
        return 'alphaplot_{obs}_nf_{nf}'.format(**locals())
    return save_figures(plots.plot_alphaS, summed_table, output_dir,
                        prefix=prefix, fmt=fmt, namefunc=namefunc)

def save_asq(pdfsets, output_dir, prefix, fmt='pdf'):
    import smpdflib.plots as plots
    def namefunc(nf):
        return 'alphaSQ_nf_{nf}'.format(**locals())
    return save_figures(plots.plot_asQ, pdfsets, output_dir,
                        prefix=prefix, fmt=fmt, namefunc=namefunc)

def save_nf(summed_table, output_dir, prefix, fmt='pdf'):
    """
    Generate plots showing the value of the observables as a function
    of N_f. The value is obtained by summing each bin in the applgrid."""

    #slow to import
    import smpdflib.plots as plots
    def namefunc(obs, bin_, oqcd):
        obs = normalize_name(obs)
        return 'nfplot_{obs}_{oqcd}'.format(**locals())
    return save_figures(plots.plot_nf, summed_table, output_dir,
                        prefix=prefix, fmt=fmt, namefunc=namefunc)
@check(check_know_errors)
def save_obscorrs(data_table, output_dir, prefix, base_pdf=None, fmt='pdf'):

    import smpdflib.plots as plots
    def namefunc(pdf):
        return 'obs_corrs_%s' % pdf
    return save_figures(plots.plot_observable_correlations, data_table,
                        output_dir, base_pdf = base_pdf,
                        prefix=prefix, fmt=fmt, namefunc=namefunc)


@check(check_know_errors)
def export_html(total, output_dir, prefix):
    """
    Export results as a rich HTML table."""
    import smpdflib.utils as utils
    import smpdflib.core as lib
    filename = "%sresults.html" % (prefix if prefix else '')
    utils.save_html(total[lib.DISPLAY_COLUMNS], osp.join(output_dir, filename))

@check(check_know_errors)
def export_obscorrs(data_table, output_dir, prefix, base_pdf=None):
    import pandas as pd

    from smpdflib.corrutils import observable_correlations

    for title, corrmat, labels in observable_correlations(data_table, base_pdf):
        name = prefix + "_" + normalize_name(title)
        filename = osp.join(output_dir, name + ".csv")
        pd.DataFrame(corrmat, index=labels, columns=labels).to_csv(filename, sep='\t')

#TODO: Ability to import exported csv
@check(check_know_errors)
def export_csv(total, output_dir, prefix):
    """
    Export results as a CSV so they can be processed by other tools. The
    resulting file is tab-separated."""
    filename = "%sresults.csv" % (prefix if prefix else '')
    total.to_csv(osp.join(output_dir, filename), sep='\t')

#TODO: Think how to make this better
def test_as_linearity(summed_table, diff_from_line = 0.25):
    """
    Test linearity of value of the value of the observable as a function of
    as. If th test fails for some point, report it in the tables and in
    subsequent plots."""
    import smpdflib.core as lib
    return lib.test_as_linearity(summed_table, diff_from_line = diff_from_line)

def save_correlations(results, output_dir, prefix, fmt='pdf'):
    """Plot correlations between PDFs and observables"""
    import smpdflib.plots as plots
    def namefunc(obs,pdf):
        return "smpdfplot_%s_%s"%(pdf, obs)
    return save_figures(plots.plot_correlations, results, output_dir,
                        prefix=prefix,
                        fmt=fmt, namefunc=namefunc)


def _mc2hname(prefix, pdf, group, config):
    return '_'.join((prefix, str(pdf), str(group['Neig'])))



_namemap = {'mc2hessian': _mc2hname}

def gen_gridnames(action, group, config):
    from smpdflib import lhaindex
    prefix = group['prefix']
    if 'grid_names' not in group:
        group['grid_names'] = {}
    for pdf in group['pdfsets']:
        if action in _namemap:
            grid_name = _namemap[action](prefix, pdf, group, config)
        else:
            grid_name = prefix + action + '_' + str(pdf)
        if lhaindex.isinstalled(grid_name):
            raise ActionError("The grid '%s' "
                              "that will be generated by the action '%s' "
                              "is already installed in the LHAPDF path." %
                              (grid_name, action))
        if grid_name in config._grid_names:
            raise ActionError("The grid %s would be generated multiple times."
                               % grid_name)
        config._grid_names.append(grid_name)
        group['grid_names'][(action, pdf)] = grid_name




@check(gen_gridnames)
def create_smpdf(data_table, output_dir, grid_names ,smpdf_tolerance=0.05,
                 Neig_total = 200, full_grid=False, db=None,
                 smpdf_correlation_threshold=None):
    
    from smpdflib.corrutils import DEFAULT_CORRELATION_THRESHOLD
    from smpdflib.reducedset import create_smpdf
    if smpdf_correlation_threshold is None:
        smpdf_correlation_threshold = DEFAULT_CORRELATION_THRESHOLD
    gridpaths = []
    for (pdf, pdf_table) in data_table.groupby('PDF'):
        pdf_results = pdf_table.Result.unique()
        result = create_smpdf(pdf, pdf_results, output_dir,
                                  grid_names[('smpdf', pdf)],
                                  smpdf_tolerance=smpdf_tolerance,
                                  Neig_total = Neig_total,
                                  full_grid=full_grid,
                                  correlation_threshold=smpdf_correlation_threshold,
                                  db=db)
        gridpaths.append(result)
    return gridpaths

@require_args('sample_Q', 'Neig')
@check(gen_gridnames)
def create_mc2hessian(pdfsets, Neig ,output_dir, sample_Q, grid_names,
                      db=None):
    import smpdflib.reducedset as lib
    gridpaths = []
    for pdf in pdfsets:
        result = lib.create_mc2hessian(pdf, Q=sample_Q, Neig=Neig,
                                       output_dir=output_dir,
                                       name=grid_names[('mc2hessian', pdf)],
                                       db=db)
        gridpaths.append(result)
    return gridpaths


def check_lhawrite(action, group, config):
    import smpdflib.lhaindex
    path = smpdflib.lhaindex.get_lha_datapath()
    if not os.access(path, os.W_OK):
        raise ActionError("Cannot write in LHAPDF path: %s. "
                          "Unable to complete action %s" %(path, action))

@check(check_lhawrite)
def install_grids(grid_names, output_dir):
    import smpdflib.lhaindex
    dest = smpdflib.lhaindex.get_lha_datapath()
    for name in grid_names.values():
        shutil.copytree(osp.join(output_dir, name), osp.join(dest, name))


ACTION_DICT = OrderedDict((
               ('testas',  test_as_linearity),
               ('violinplots',save_violins),
               ('ciplots',save_cis),
               ('asplots',save_as),
               ('asQplots',save_asq),
               ('nfplots',save_nf),
               ('obscorrplots', save_obscorrs),
               ('exporthtml', export_html),
               ('exportcsv', export_csv),
               ('exportobscorrs', export_obscorrs),
               ('plotcorrs', save_correlations),
               ('smpdf', create_smpdf),
               ('mc2hessian', create_mc2hessian),
               ('installgrids', install_grids),
               ))


REALACTIONS = set(ACTION_DICT.keys())

METAACTION_DICT = {'all': (REALACTIONS, "Implies all other actions."),
                   'savedata': ({'exportcsv', 'exporthtml'}, "Export html and "
                                                                   "csv.")
                  }
METAACTIONS = set(METAACTION_DICT.keys())

SAVEDATA = {'exporthtml', 'exportcsv'}

ACTIONS = REALACTIONS | METAACTIONS


#TODO: Do this better
def requires_result(action):
    args = inspect.getargspec(ACTION_DICT[action]).args
    possibles = set(['results', 'summed_table', 'data_table', 'pdfcorrlist',
                     'total'])
    return bool(possibles & set(args))

#TODO: Do this better
def requires_correlations(action):
    args = inspect.getargspec(ACTION_DICT[action]).args
    return ('pdfcorrlist' in args)

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

def do_actions(acts, resources):
    for action in ACTION_DICT.keys():
        if not action in acts:
            continue
        func = ACTION_DICT[action]
        fargs = inspect.getargspec(func).args
        args = {k:v for k,v in resources.items() if k in fargs}
        result = func(**args)
        if hasattr(func, 'provides'):
            resources[func.provides] = result
        elif hasattr(func, 'appends_to'):
            key = func.appends_to
            if key in resources:
                resources[key] += list(result)
            else:
                resources[key] = list(result)
        yield action, result

def execute_config(conf, output_dir, db):

    import pandas as pd
    import smpdflib.core as lib
    import logging

    resultset = []
    for group in conf.actiongroups:
        pdfsets, observables = group['pdfsets'], group['observables']
        resources = group.copy()
        resources.pop('actions')
        # perform convolution
        #TODO Do this better
        if any(requires_result(act) for act in group['actions']):
            results = lib.produce_results(pdfsets, observables, db)
            resultset.append(results)
            data_table = lib.results_table(results)
            summed_table = lib.summed_results_table(results)

            total = pd.concat((data_table,
                                summed_table),
                                ignore_index = True)
            if logging.getLogger().isEnabledFor(logging.DEBUG):
                ...
                #print_results(results)
            resources.update({'results':results, 'data_table':data_table,
                           'total':total, 'summed_table':summed_table})

        if any(requires_correlations(act) for act in group['actions']):
            pdfcorrlist = lib.correlations(data_table, db=db)
            resources.update({'pdfcorrlist':pdfcorrlist})


        prefix = group['prefix']
        resources.update({ 'output_dir':output_dir,
                       'prefix':prefix,
                       'pdfsets': pdfsets,
                       'db': db})
        for action, res in do_actions(group['actions'], resources):
            logging.info("Finalized action '%s'." % action)
    return resultset

#TODO: Rethink this
def print_results(results):
    for result in results:
        for member in result.iterreplicas():
            print ("\n- %s replica %d"% (result.pdf,member))
            print ("- APPLgrid convolution results:")
            for i, val in enumerate(result[member]):
                print ("\tData bin %i: %e" % (i, val))
    print ("\n +--------+ Completed +--------+\n")

def helptext(action):
    if action in METAACTION_DICT:
        return METAACTION_DICT[action][1]
    func = ACTION_DICT[action]
    doc = func.__doc__
    if doc:
        #textwrap is worse than useless...
        doc =  textwrap.fill(textwrap.dedent(func.__doc__).replace('\n', ' '),
                             subsequent_indent='\t')
    else:
        doc = ''
    spec = inspect.getargspec(func)
    if spec.defaults:
        params = spec.args[-len(spec.defaults):]
        param_doc = ' The optional parameters fot this action are:\n\t\t'
        param_doc += '\n\t\t'.join("* %s (default: %s)" % (param, default) for
                     param, default in zip(params, spec.defaults))
    else:
        param_doc = ''
    return doc + param_doc



def gen_docs():
    s = "Possible actions are:\n%s" % '\n'.join('\t- %s : %s\n' %
        (action, helptext(action))
        for action in list(ACTION_DICT.keys()) + list(METAACTION_DICT.keys()))
    return s
