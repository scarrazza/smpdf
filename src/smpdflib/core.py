# -*- coding: utf-8 -*-
from __future__ import division

"""
Created on Tue Apr 28 12:13:00 2015

@author: zah
"""
""" SMPDF """

__author__ = 'Stefano Carrazza'
__license__ = 'GPL'
__version__ = '1.0.0'
__email__ = 'stefano.carrazza@mi.infn.it'

import os
import os.path as osp
import sys
from collections import defaultdict, OrderedDict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats
from pandas.stats import ols


import smpdflib.lhaindex as lhaindex
import smpdflib.plotutils as plotutils
try:
    sys.path.append('applwrap')
    from applwrap import initpdf, initobs, pdfreplica, convolute
except ImportError:
    os.system("make -C applwrap")
    from applwrap import initpdf, initobs, pdfreplica, convolute

ORDERS_QCD = {0: 'LO', 1: 'NLO', 2: 'NNLO'}

#for N_f = 4, LHAPDF's M_Z is actually M_{charm}
M_REF = defaultdict(lambda: 'Z', {4:'c'})




class TupleComp(object):
    def __hash__(self):
        return hash(self.get_key())

    def __eq__(self, other):
        return (isinstance(other, self.__class__)
                and self.get_key() == other.get_key())

class BaseObservable(TupleComp):
    def __init__(self, name, order):
        self.name = name
        self.order = order

    def __str__(self):
        return "%s(%s)"%(self.name, ORDERS_QCD[self.order])

    def __repr__(self):
        return "<%s:%s>" % (self.__class__.__name__, self.__str__())

    def get_key(self):
        return (self.name, self.order)

class Observable(BaseObservable):
    def __init__(self, filename, order):
        self.filename = filename
        self.order = order

    @property
    def name(self):
        return osp.splitext(osp.basename(self.filename))[0]


class PDF(TupleComp):
    def __init__(self, name):
        self.name = name

    def get_key(self):
        #Convert python2 unicode to string so no u'prefix' is printed
        return (str(self.name),)

    def __str__(self):
        return self.name

    def __repr__(self):
        return "<%s:%s>" % (self.__class__.__name__, self.name)

    @property
    def oqcd_str(self):
        return ORDERS_QCD[self.OrderQCD]

    @property
    def collaboration(self):
        return lhaindex.get_collaboration(self.name)

    @property
    def as_from_name(self):
        return lhaindex.as_from_name(self.name)

    @property
    def mref(self):
        return "M_%s" % M_REF[self.NumFlavors]

    @property
    def reps(self):
        return range(self.NumMembers)

    def __getattr__(self, name):
        return lhaindex.parse_info(self.name)[name]


class Result():
    def __init__(self, obs, pdf, data):
        self.obs = obs
        self.pdf = pdf
        self._data = pd.DataFrame(data)

    @property
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

    @property
    def nbins(self):
        return self._all_vals.shape[0]

    def std_error(self, nsigma=1):
        raise NotImplementedError("No error computation implemented for this"
                                  "type of set")

    def sample_values(self, n):
        raise NotImplementedError("No sampling implemented for this"
                                  "type of set")

    def std_interval(self, nsigma=1):
        std = self.std_error(nsigma)
        return pd.DataFrame({'min':self._cv - std,
                             'max':self._cv + std})

    def rel_std_interval(self, nsigma=1):
        std = self.std_error(nsigma)
        return pd.DataFrame({'min':-std,
                             'max':std})


    def __getitem__(self, item):
        return self._data[item]

    #TODO: Should this be the default iterator?
    def iterall(self):
        return iter(self._data)

    def _violin_data(self, rel_to=None):
        absdata = pd.concat(self.sample_values(10000),axis=1)
        if rel_to is None:
            rel_to = 1
        reldata = absdata.as_matrix().T/rel_to
        return reldata

    def violin_plot(self, data=None , **kwargs):
        if data is None:
            data = self._violin_data()

        myargs = {'label': str(self.pdf)}
        myargs.update(kwargs)
        return plotutils.violin_plot(data, **myargs)

    def sumbins(self, bins = None):
        sumobs = BaseObservable(self.obs.name + '[Sum]', self.obs.order)
        data = pd.DataFrame(self._data.sum(axis=0)).T
        return self.__class__(sumobs, self.pdf, data)



class SymHessianResult(Result):

    def std_error(self, nsigma=1):
        diffsq = (self._all_vals.subtract(self._cv, axis=0))**2
        return diffsq.sum(axis=1).apply(np.sqrt)*nsigma
    @property
    def errorbar68(self):
        return self.rel_std_interval()

    def sample_values(self, n):
        diffs = self._all_vals.subtract(self._cv, axis=0)
        for _ in range(n):
            weights = np.random.normal(size=self.nrep)
            error = (diffs*weights).sum(axis=1)
            yield self._cv + error

    def _violin_data(self, rel_to = None):
        std = self.std_error()
        mean = self.central_value.as_matrix()

        if rel_to is None:
            rel_to = np.ones_like(mean)
        vpstats = []
        for m,s,r in zip(mean,std, rel_to):
            # Dictionary of results for this distribution
            stats = {}

            # Calculate basic stats for the distribution
            min_val = m - 3*s
            max_val = m + 3*s
            coords = np.linspace(m - 3*s, m+3*s, 1000)
            # Evaluate the kernel density estimate
            stats['vals'] = scipy.stats.norm(m,s).pdf(coords)*r
            stats['coords'] = coords/r

            # Store additional statistics for this distribution
            stats['mean'] = m/r
            stats['median'] = m/r
            stats['min'] = min_val/r
            stats['max'] = max_val/r

            # Append to output
            vpstats.append(stats)

        return vpstats




class MCResult(Result):
    #TODO: Is it correct to consider each bin as independant here?
    def centered_interval(self, percent=68, addcentral=True):
        n = percent*self.nrep//100
        def get_lims(row):
            row = row.as_matrix()
            s = np.argsort(np.abs(row))
            sel = row[s][:n]
            return pd.Series({'min':np.min(sel), 'max':np.max(sel)})

        diffs = self._all_vals.subtract(self._cv, axis=0)
        limits = diffs.apply(get_lims, axis=1)
        if addcentral:
            limits = limits.add(self._cv, axis=0)
        return limits

    @property
    def errorbar68(self):
        return self.centered_interval(addcentral=False)

    def std_error(self, nsigma=1):
        return self._all_vals.std(axis=1)*nsigma

    def sample_values(self, n):
        for _ in range(n):
            col = np.random.choice(self._all_vals.columns)
            yield self._all_vals[col]

    def _violin_data(self, rel_to = None):
        if rel_to is None:
            rel_to = 1
        return self._all_vals.as_matrix().T/ rel_to

def aggregate_results(results):
    combined = defaultdict(lambda: {})
    for result in results:
        combined[result.obs][result.pdf] = result
    return combined

DISPLAY_COLUMNS = ['Observable', 'PDF', 'Bin', 'CV', 'Up68', 'Down68',
                   'Remarks']

def results_table(results):
    records = pd.concat([pd.DataFrame(OrderedDict([
                ('Observable'       , result.obs),
                ('PDF'              , result.pdf),
                ('Collaboration'    , result.pdf.collaboration),
                ('as_from_name'     , result.pdf.as_from_name),
                ('alpha_sMref'      , result.pdf.AlphaS_MZ),
                ('PDF_OrderQCD'     , result.pdf.oqcd_str),
                ('NumFlavors'       , result.pdf.NumFlavors),
                ('Bin'              , np.arange(1, result.nbins + 1)),
                ('CV'               , result.central_value),
                ('Up68'             , np.abs(result.errorbar68['max'])),
                ('Down68'           , np.abs(result.errorbar68['min'])),
                ('Remarks'          , None),
                ('Result'           , result),
               ])) for result in results],
               ignore_index=True)
    #Must be an independent list for each record
    records['Remarks'] = records.apply(lambda x: [], axis=1)
    return records

def summed_results_table(results):

    if isinstance(results, pd.DataFrame):
        results = results['Result'].unique()
    table = results_table([result.sumbins() for result in results])
    table['Bin'] = 'sum'
    return table

def compare_violins(results, base_pdf = None):
    if not isinstance(results, dict):
        combined = aggregate_results(results)
    else:
        combined = results
    for obs in combined:
        figure = plt.figure()
        norms = None
        handles = []
        plt.title(str(obs))
        #ncolors = len(combined[obs])
        colorlist = ['#222222','#ff0000', '#00ff00', '#0000ff',]
        colors  = plotutils.color_names_to_rgb(colorlist)
        #colors = itertools.cycle(plotutils.get_set1_colors(
        #                         ncolors).mpl_colors)
        alpha = 1
        base = combined[obs].get(base_pdf, None)
        results = sorted(combined[obs].values(), key = lambda x: x!=base)
        for result in results:
            if base is not None:
                cv = base.central_value.as_matrix()
                data = result._violin_data(rel_to=cv)
            else:
                data = data = result._violin_data()
            color = next(colors) + (alpha,)
            alpha /= 2
            plot, handle, norms = result.violin_plot(data, color=color,
                                              showextrema=False,
                                              normvalues=norms)
            handles.append(handle)
        plt.xlabel('bins')
        if base_pdf:
            plt.ylabel('Rel to %s' % base_pdf)
        else:
            plt.ylabel("Observable value")
        plt.xticks(range(1,len(result.central_value) + 1))
        plt.legend(handles=handles)
        yield (obs,), figure

def compare_cis(results, base_pdf = None):
    if not isinstance(results, dict):
        combined = aggregate_results(results)
    else:
        combined = results
    for obs in combined:
        figure = plt.figure()
        plt.title(str(obs))
        base = combined[obs].get(base_pdf, None)
        results = sorted(combined[obs].values(), key = lambda x: x!=base)
        l = len(results)
        if l < 20:
            delta = iter(np.linspace(-0.0125*l, 0.0125*l, l))
        else:
            delta = iter(np.linspace(0.25, 0.25, l))
        #x = iter(np.arange(1, len(results) + 1))
        for result in results:
            x = np.arange(1, result.nbins+1) + next(delta)
            data_cv = result.central_value.as_matrix().copy()
            data_ci = result.errorbar68.as_matrix().copy()
            if base is not None:
                base_cv = base.central_value.as_matrix()
                data_cv /= base_cv
                data_ci /= base_cv[:,np.newaxis]
            data_ci = np.abs(data_ci)
            plt.errorbar(x, data_cv, yerr=data_ci.T,
                                     linestyle='none',
                                     label=result.pdf, elinewidth = 2,
                                     capsize=5)
        plt.xlabel('bins')
        if base_pdf:
            plt.ylabel('Rel to %s' % base_pdf)
        else:
            plt.ylabel("Observable value")
        plt.xticks(range(1,len(result.central_value) + 1))
        plt.legend()
        plt.xlim(0.5, results[0].nbins + 0.5)
        plt.grid(axis='x')
        yield (obs,), figure

@plotutils.ax_or_gca
def plot_remarks(df, ax=None):
    have_remarks = df[df['Remarks'].apply(len) > 0]

    if len(have_remarks):
            ax.plot(have_remarks['alpha_sMref'], have_remarks['CV'],
                 'ro', markersize = 20, fillstyle = 'none',
                 markeredgewidth = 5,
                 label="Problematic points")

def process_label(process, bin_):
    if bin_ == 'sum':
        return str(process)
    else:
        return "%s[bin:%d]"% (process, bin_)

#TODO: Abstract groupbyplots away? Tried, but seems too hard...
def plot_alphaS(results_table):
    df = results_table.sort('alpha_sMref')
    for (process, nf, bin_), process_df in df.groupby(['Observable',
                                                    'NumFlavors', 'Bin']):
        fig = plt.figure()


        for (oqcd,col), col_df in process_df.groupby(['PDF_OrderQCD',
                                                      'Collaboration']):
            label = "%s (%s)" % (col, oqcd)

            plt.errorbar(col_df['alpha_sMref'], col_df['CV'],
                         yerr = np.array(col_df['Down68'],
                                         col_df['Up68']),
                        label = label, linestyle='-', marker = 's')


        plot_remarks(process_df)
        plt.xlabel(r'$\alpha_S(M_%s)$' % M_REF[nf])
        plt.ylabel(r'Value of observable')
        xran = plotutils.extend_range(process_df['alpha_sMref'].min(),
                            process_df['alpha_sMref'].max())
        plt.xlim(*xran)
        plt.legend()
        plt.title("%s $N_f=$%d" % (process_label(process, bin_), nf), y = 1.08)
        plt.tight_layout()
        yield (process, nf, bin_),fig

def plot_nf(results_table):
    df = results_table.sort('NumFlavors')
    for (process, bin_, oqcd), process_df in df.groupby(['Observable',
                                                    'Bin', 'PDF_OrderQCD']):
        fig = plt.figure()


        for (col, asn), pdf_df in process_df.groupby(['Collaboration',
                                                      'as_from_name']):
            label = "%s(as: %s)"%(col, asn)

            plt.errorbar(pdf_df['NumFlavors'], pdf_df['CV'],
                         yerr = np.array(pdf_df['Down68'],
                                         pdf_df['Up68']),
                        label = label, linestyle='-', marker = 's')


        plot_remarks(process_df)
        plt.xlabel(r'$N_f$')
        plt.ylabel(r'Value of observable')
        xran = plotutils.extend_range(process_df['NumFlavors'].min(),
                            process_df['NumFlavors'].max())
        plt.xlim(*xran)
        plt.xticks(process_df['NumFlavors'].unique())
        plt.legend()
        plt.title("%s PDFs, %s" % (oqcd, process_label(process, bin_)), y=1.08)
        plt.tight_layout()
        plt.grid(axis='x')
        yield  (process, bin_, oqcd),fig

def plot_asQ(pdfsets):
    df = pd.DataFrame([{'NumFlavors':pdf.NumFlavors, 'PDF':pdf} for
                        pdf in pdfsets])
    for nf, gdf in df.groupby(['NumFlavors']):
        fig = plt.figure()
        for pdf in gdf.PDF:
            plt.plot(pdf.AlphaS_Qs, pdf.AlphaS_Vals, label=pdf.name)
            plt.ylabel(r'$\alpha_S$')
            plt.xlabel(r'Q(GeV)')
            plt.xscale('log')
            plt.title('$N_f=%d$' % nf)
            plt.legend()
        yield (nf,), fig


RESULT_TYPES = defaultdict(lambda:Result,
                           symmhessian = SymHessianResult,
                           replicas   = MCResult,
                           )

def make_result(obs, pdf, datas):
    error_type = pdf.ErrorType
    return RESULT_TYPES[error_type](obs, pdf, datas)


def make_convolution(pdf, observables):
    datas = defaultdict(lambda:OrderedDict())
    #TODO: load many replicas in C++
    #TODO: Could we loop over observables and then over memebers?
    if observables:
        initpdf(pdf.name)
    for obs in observables:
        initobs(obs.filename)
        for rep in pdf.reps:
            #TODO: hide this call from the api, do in convolute.
            sys.stdout.write('\r-> Computing replica %d of %s' %
                             (rep, pdf))
            sys.stdout.flush()
            pdfreplica(rep)
            res = convolute(obs.order)
            datas[obs][rep] = np.array(res)
        sys.stdout.write('\n')
    return datas

#TODO: Merge this with results_table
def results_from_datas(dataset):
    results = []
    for pdf in dataset:
        data = dataset[pdf]
        results += [make_result(obs, pdf, data[obs]) for obs in data]
    return results

#TODO: Refactor this after adding efficient convolution
def get_dataset(pdfsets, observables, db=None):
    def make_key(pdf, obs):
        return str((pdf.get_key(), obs.get_key()))
    dataset = OrderedDict()
    for pdf in pdfsets:
        #bool(db) == False if empty
        if db is not None:
            res = {}
            obs_to_compute = []
            for obs in observables:
                key = make_key(pdf, obs)
                if key in db:
                    res[obs] = db[key]
                else:
                    obs_to_compute.append(obs)

            computed_data = make_convolution(pdf, obs_to_compute)
            for newobs in computed_data:
                key = make_key(pdf, newobs)
                db[key] = computed_data[newobs]
            res.update(computed_data)


        else:
            res = make_convolution(pdf, observables)

        dataset[pdf] = res
    return dataset

def convolve_or_load(pdfsets, observables, db=None):
    #results = []
    results = results_from_datas(get_dataset(pdfsets, observables, db))
    return results

#TODO: Move somewhere else
def save_html(df, path):
    import jinja2
    import codecs

    env = jinja2.Environment(loader = jinja2.PackageLoader('smpdflib',
                                                           'templates'))
    template = env.get_template('results.html.template')
    def remark_formatter(remarks):
        if not remarks:
            return ''
        else:
            return '<ul>%s</ul>' % '\n'.join('<li>%s</li>' %
                   jinja2.escape(remark) for remark in remarks)

    #http://stackoverflow.com/questions/26277757/pandas-to-html-truncates-string-contents
    with pd.option_context('display.max_colwidth', -1):
        table = df.to_html(
                             formatters={'Remarks':remark_formatter},
                             escape = False)
    result = template.render(table=table)
    with codecs.open(path, 'w', 'utf-8') as f:
        f.write(result)

def test_as_linearity(summed_table, diff_from_line = 0.25):
    group_by = ('Observable','NumFlavors', 'PDF_OrderQCD', 'Collaboration')
    for (process, nf, oqcd, col), curve_df in summed_table.groupby(group_by):
        if len(curve_df) <= 2:
            continue
        fit = ols.OLS(y=curve_df['CV'], x=curve_df['alpha_sMref'],
                      weights=1/curve_df['CV']**2)

        diff = (fit.y_predict - curve_df['CV'])/curve_df['Up68']
        bad = diff > diff_from_line
        for ind in curve_df[bad].index:
                remark = (u"Point away from linear fit by %1.1fσ" %
                                diff.ix[ind])
                summed_table.loc[ind,'Remarks'].append(remark)

def corrcoeff(obs, pdf):
    return len(pdf)/(len(pdf)-1)*\
           (np.mean(obs*pdf) - np.mean(obs)*np.mean(pdf))/\
           (np.std(obs,ddof=1)*np.std(pdf,ddof=1))

def correlations(results):
    """Includes from mc2hessian library"""
    from mc2hlib.common import load_pdf, compress_X_abs
    from mc2hlib.lh import hessian_from_lincomb

    # TODO: automatization
    threshold = 0.5
    for result in results:
        #TODO: select input Q
        # allocate pdf set
        pdf, fl, xgrid = load_pdf(str(result.pdf), 100.0)
        nf = fl.n
        nx = xgrid.n

        # praparing
        figure, axarr = plt.subplots(fl.n, sharex=True, sharey=True, figsize=(8, fl.n+3))
        cc = np.zeros(shape=(fl.n, xgrid.n))
        for bin in range(len(result[:])):
            obs = np.zeros(pdf.n_rep)

            for rep in range(1,pdf.n_rep+1):
                obs[rep-1] = result[rep][bin]

            for f in range(fl.n):
                for x in range(xgrid.n):
                    cc[f,x] = corrcoeff(obs, pdf.xfxQ[:,f,x])

                axarr[f].plot(xgrid.x, cc[f])
                axarr[f].axhline(threshold, c='r', ls='--')
                axarr[f].axhline(-threshold, c='r', ls='--')
                axarr[f].set_ylim([-1,1])
                axarr[f].set_xscale('log')
                axarr[f].set_ylabel("pdg: " + str(fl.id[f]))

            axarr[0].set_title(str(result.obs) + "\n")
            plt.xlabel("x")
            figure.subplots_adjust(hspace=0)
            plt.setp([a.get_xticklabels() for a in figure.axes[:-1]], visible=False)

        yield result.obs, figure

        # Step 1: create pdf covmat
        print "\n- Building PDF covariance matrix:"
        X = (pdf.xfxQ.reshape(pdf.n_rep, nx*nf) - pdf.f0.reshape(nx*nf)).T
        print " [Done] "

        mask = (np.abs(cc.reshape(fl.n*xgrid.n)) >= threshold)
        print (" [Info] Keeping %d nf*nx with threshold = %f" % (np.count_nonzero(mask), threshold))

        X = X[mask,:]
        # Step 2: solve the system
        print "\n- Quick test:"
        for neig in range(1, pdf.n_rep):
            vec, cov = compress_X_abs(X, neig)

            stdh = iter(np.sqrt(np.diag(cov)))

            # Step 3: quick test
            rmask = mask.reshape(fl.n, xgrid.n)
            est = Norm = 0
            for f in range(fl.n):
                for x in range(xgrid.n):
                    if rmask[f,x]:
                        t0 = pdf.std[f,x]
                        t1 = next(stdh)

                        #print "1-sigma MonteCarlo (fl,x,sigma):", fl.id[f], xgrid.x[x], t0
                        #print "1-sigma Hessian    (fl,x,sigma):", fl.id[f], xgrid.x[x], t1
                        #print "Ratio:", t1/t0
                        est += abs(pdf.f0[f,x] * (1-t1/t0))
                        Norm += abs(pdf.f0[f,x])

            est /= Norm
            # TODO: is this the best estimator? and cut criteria?
            print ("Neig %3d, estimator: %e" % (neig, est))
            if est <= 1e-3: break;

        # Step 4: exporting to LHAPDF
        print "\n- Exporting new grid..."
        hessian_from_lincomb(pdf, vec, set_name="smpdf_" + str(result.pdf))
