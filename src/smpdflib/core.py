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

import os.path as osp
import sys
from collections import defaultdict, OrderedDict, namedtuple
import numbers
import multiprocessing

import numpy as np
import pandas as pd
import yaml
import scipy.stats
from pandas.stats import ols


import smpdflib.lhaindex as lhaindex
import smpdflib.plotutils as plotutils

import applwrap

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
    _meanQ = None
    _nbins = None
    def __init__(self, filename, order):
        self.filename = filename
        self.order = order

    @property
    def meanQ(self):
        return self._meanQ

    @property
    def name(self):
        return osp.splitext(osp.basename(self.filename))[0]

_selected_grid = None

class APPLGridObservable(Observable):


    @property
    def nbins(self):
        if self._nbins is not None:
            return self.nbins
        with self:
            nbins = applwrap.getnbins()
        self._nbins = nbins
        return nbins

    @property
    def meanQ(self):
        if self._meanQ is not None:
            return self._meanQ
        with self:
            meanQ = [applwrap.getobsq(self.order, i) for
                 i in range(self.nbins)]
        self._meanQ = meanQ
        return meanQ

    def __enter__(self):
        """Load observable file in memory"""
        global _selected_grid
        if _selected_grid == self.filename:
            return
        if _selected_grid is not None:
            raise RuntimeError("Contrdicting observable scope. "
                               "Was %s and trying to enter %s" %
                               (_selected_grid, self.filename))
        applwrap.initobs(self.filename)
        _selected_grid = self.filename
    #TODO: Unload Observable here
    def __exit__(self, exc_type, exc_value, traceback):
        global _selected_grid
        _selected_grid = None

class PredictionObservable(Observable):
    def __init__(self, filename):
        self.filename = filename
        with open(filename) as f:
            d = yaml.load(f)
        #TODO: All checking
        self._params = d

    def to_result(self, pdfset):
        if str(pdfset) not in self.pdf_predictions:
            raise ValueError("No predictions found for pdf %s" % pdfset)
        path = self.pdf_predictions[str(pdfset)]
        if not osp.isabs(path):
            path = osp.join(osp.dirname(self.filename), path)
        #TODO: Transpose all result dataframes
        datas = pd.DataFrame.from_csv(path, sep='\t', index_col=0).T
        return make_result(self, pdfset, datas)

    @property
    def meanQ(self):
        if isinstance(self.energy_scale, numbers.Number):
            return [self.energy_scale]*self.nbins
        else:
            return self.energy_scale

    def __getattr__(self, attr):
        return self._params[attr]



class PDF(TupleComp):
    def __init__(self, name):
        self.name = name

    def get_key(self):
        #Convert python2 unicode to string so no u'prefix' is printed
        return (str(self.name),)

    def __enter__(self):
        """Load PDF file in memory"""
        applwrap.initpdf(self.name)

    #TODO: Unload PDF here
    def __exit__(self, exc_type, exc_value, traceback):
        pass

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
        if name.startswith('__'):
            raise AttributeError()
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

    #TODO: This should really be done in Observable. Can we access nbins?
    @property
    def meanQ(self):
        return self.obs.meanQ

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

class HessianResult(SymHessianResult):

    def std_error(self, nsigma=1):
        m = self._all_vals.as_matrix()
        diffsq = (m[:, ::2] - m[:, 1::2])**2
        return np.sqrt(diffsq.sum(axis=1))/2.0*nsigma

    def sample_values(self, n):
        m = self._all_vals.as_matrix()
        plus = m[:, ::2]
        minus = m[:, 1::2]

        for _ in range(n):
            r = np.random.normal(size=len(plus))
            error = (r >=0)*r*plus - (r < 0)*r*minus
            yield self._cv + error


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
    combined = defaultdict(lambda: OrderedDict())
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

RESULT_TYPES = defaultdict(lambda:Result,
                           symmhessian = SymHessianResult,
                           hessian = HessianResult,
                           replicas   = MCResult,
                           )

def make_result(obs, pdf, datas):
    error_type = pdf.ErrorType
    return RESULT_TYPES[error_type](obs, pdf, datas)

def make_observable(name, *args, **kwargs):
    extension = osp.splitext(name)[-1]
    prediction_extensions = ('.yaml', '.yml', '.info', '.txt')
    applgrid_extensions = ('.root',)
    if extension in prediction_extensions:
        return PredictionObservable(name, *args, **kwargs)
    elif extension in applgrid_extensions:
        return APPLGridObservable(name, *args, **kwargs)
    else:
        raise ValueError("Only files with extensions: %s "
                         "are valid observables" % str(prediction_extensions
                                                   + applgrid_extensions))


def convolve_one(pdf, observable):
    import applwrap
    from smpdflib.core import PDF, APPLGridObservable
    print(pdf)
    print(observable)
    res = {}
    with pdf, observable:
        for rep in pdf.reps:
            applwrap.pdfreplica(rep)
            res[rep] = np.array(applwrap.convolute(observable.order))
    return res

def _convolve_one_args(args):
    return convolve_one(*args)


def make_convolution(pdf, observables):
    datas = defaultdict(lambda:OrderedDict())
    #TODO: load many replicas in C++
    #TODO: Could we loop over observables and then over memebers?
    if not observables:
        return {}

    with(pdf):
        for obs in observables:
            with obs:
                print(pdf)
                print(obs)
                for rep in pdf.reps:
                    sys.stdout.write('\r-> Computing replica %d of %s' %
                                     (rep, pdf))
                    sys.stdout.flush()
                    applwrap.pdfreplica(rep)
                    res = applwrap.convolute(obs.order)
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

def get_dataset_parallel(pdfsets, observables, db=None):
    def make_key(pdf, obs):
        return str((pdf.get_key(), obs.get_key()))
    n_cores = multiprocessing.cpu_count()
    dataset = OrderedDict()
    to_compute =  []
    for pdf in pdfsets:
        dataset[pdf] = OrderedDict()
        for obs in observables:
            if db is not None:
                key = make_key(pdf, obs)
                if key in db:
                    dataset[pdf][obs] = db[key]
                else:
                    to_compute.append((pdf, obs))
            else:
                to_compute.append((pdf, obs))

    def break_in(to_compute, n_cores):
       i = 0
       while i < len(to_compute):
           yield to_compute[i:i + n_cores]
           i += n_cores

    for convs in break_in(to_compute, n_cores):
        pool = multiprocessing.Pool(processes=len(convs))
        results = pool.map(_convolve_one_args, convs)
        pool.close()
        for ((pdf, obs), result) in zip(convs, results):
            dataset[pdf][obs] = result
            if db is not None:
                db[make_key(pdf, obs)] = result

    return dataset






def convolve_or_load(pdfsets, observables, db=None):
    #results = []
    results = results_from_datas(get_dataset_parallel(pdfsets, observables, db))
    return results

def produce_results(pdfsets, observables, db=None):
    predictions = [obs for obs in observables if
                   isinstance(obs, PredictionObservable)]

    applgrids = [obs for obs in observables if
                   isinstance(obs, APPLGridObservable)]


    results = (convolve_or_load(pdfsets, applgrids, db) +
               [pred.to_result(pdfset)
                for pdfset in pdfsets for pred in predictions])
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
    return (
            len(pdf)/(len(pdf)-1)*
            (np.mean(obs*pdf) - np.mean(obs)*np.mean(pdf))/
            (np.std(obs,ddof=1)*np.std(pdf,ddof=1))
            )

Corrlist = namedtuple('Corrlist', ('cc', 'threshold', 'obs', 'xgrid', 'fl'))

def compute_correlations(result, pdf,):
    """Compute correlations"""
    from mc2hlib.common import load_pdf

    lpdf, fl, xgrid = load_pdf(str(pdf), 1.0)
    cc = np.zeros(shape=(result.nbins, fl.n, xgrid.n))

    threshold = []
    for bin in range(result.nbins):

        lpdf.setQ(result.meanQ[bin])

        obs = np.array([result[rep][bin] for rep in range(1,lpdf.n_rep+1)])

        for f in range(fl.n):
            for x in range(xgrid.n):
                cc[bin,f,x] = corrcoeff(obs, lpdf.xfxQ[:,f,x])

        threshold.append( np.max(np.abs(cc[bin]))*0.5 )
    with open('log.log', 'a') as f:
        f.write("%s,%s,%s\n" % (result.pdf, result.obs, threshold))
    print("*********************")
    print(result.pdf, result.obs, threshold)
    print("*********************")



    return Corrlist(cc, threshold, result.obs, xgrid, fl)


def correlations(data_table):

    """Includes from mc2hessian library"""
    pdfcorrlist = []
    for pdf, pdf_table in data_table.groupby('PDF'):
        results = pdf_table.Result.unique()
        corrlist = []
        for result in results:
            corrlist.append(compute_correlations(result, pdf))
        pdfcorrlist += [(pdf, corrlist)]
    return  pdfcorrlist

#==============================================================================
# def compute_estimator():
#     rmask = mask.reshape(fl.n, xgrid.n)
#     est = Norm = 0
#     for f in range(fl.n):
#         for x in range(xgrid.n):
#             if rmask[f,x]:
#                 t0 = pdf.std[f,x]
#                 t1 = next(stdh)
#                 est += abs(pdf.f0[f,x] * (1-t1/t0))
#                 Norm += abs(pdf.f0[f,x])
#
#     est /= Norm
#==============================================================================
#    return est

#TODO: Integrate in mc2hessian
def optimize_hessian(X):
    from mc2hlib.common import compress_X_abs as compress_X
    std = np.std(X, axis=1)
    for neig in range(1, min(X.shape)):
        vec, cov = compress_X(X, neig)

        stdh = np.sqrt(np.diag(cov))
        # TODO: is this the best estimator? and cut criteria?
        # Step 3: quick test
        est = np.sum(np.abs(std - stdh)/std)/len(std)
        print ("Neig %3d, estimator: %e" % (neig, est))


        #est = compute_estimator()

        if est <= 1e-3:
            break
    return vec, cov

def create_smpdf(pdf, corrlist, output_dir, name,  N_eig, full_grid=False,):
    from mc2hlib.common import compress_X_abs as compress_X, load_pdf
    from mc2hlib.lh import hessian_from_lincomb
    lpdf, fl, xgrid = load_pdf(str(pdf), 1.0)

    first = corrlist[0]
    ccmax = np.zeros(shape=(first.cc.shape[1],
                            first.cc.shape[2]), dtype=bool)

    for c in corrlist:
        cc = c.cc
        threshold = c.threshold
        for bin in range(len(threshold)):
            ccmax |= (np.abs(cc[bin])>=threshold[bin])

    q2min = lpdf.pdf[0].q2Min
    lpdf.setQ(q2min)
    # Step 1: create pdf covmat
    return
    print ("\n- Building PDF covariance matrix at %f GeV:" % q2min)
    X = (lpdf.xfxQ.reshape(lpdf.n_rep, xgrid.n*fl.n) - lpdf.f0.reshape(xgrid.n*fl.n)).T
    print " [Done] "

    mask = (ccmax.reshape(fl.n*xgrid.n))
    print("[Info] Keeping %d nf*nx of %d" %
          (np.count_nonzero(mask), xgrid.n*fl.n, ))
    print("Thresholds:")
    for th in threshold:
        print ("%f"%th)

    Xm = X[mask,:]
    #TODO: Make this more efficient: compute only once the SVD.
    # Step 2: solve the system
    print "\n- Quick test:"
    vec, cov = compress_X(Xm, N_eig)

    if full_grid:
        diff = min(X.shape[0], X.shape[1]) - vec.shape[1]
        dot_vec = np.pad(vec,((0,0), (0, vec.shape[0] - vec.shape[1])),
                         mode = 'constant')
        X_comp = np.dot(X, np.eye(len(dot_vec)) - dot_vec)
        other_vec, other_cov = compress_X(X_comp, diff)
        vec = np.concatenate((vec, other_vec), axis=1)


    # Step 4: exporting to LHAPDF
    print "\n- Exporting new grid..."
    return hessian_from_lincomb(lpdf, vec, folder=output_dir,
                         set_name= name)
