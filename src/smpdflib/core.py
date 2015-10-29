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
from collections import defaultdict, OrderedDict
import contextlib
import numbers
import multiprocessing
import logging
import hashlib

import numpy as np
import pandas as pd
import yaml
import scipy.stats
import fastcache
from pandas.stats import ols

from smpdflib import lhaindex
from smpdflib import plotutils
from smpdflib.loggingutils import supress_stdout, initlogging, get_logging_queue

import applwrap


ORDERS_QCD = {0: 'LO', 1: 'NLO', 2: 'NNLO'}
NUMS_QCD = {val: key for key , val in ORDERS_QCD.items()}

#for N_f = 4, LHAPDF's M_Z is actually M_{charm}
M_REF = defaultdict(lambda: 'Z', {4:'c'})

PDG_PARTONS = {
                -3 : r"\bar{s}",
                -2 : r"\bar{u}",
                -1 : r"\bar{d}",
                 0 : r"g",
                 1 : r"d",
                 2 : r"u",
                 3 : r"s",
              }



class TupleComp(object):
    """Class whose instances compare equal if two objects have the same tuple.
    Objects also have the correct hash and can be used for dictionary keys."""
    def __hash__(self):
        return hash(self.get_key())

    def __eq__(self, other):
        return (isinstance(other, self.__class__)
                and self.get_key() == other.get_key())

#TODO: Merge with Observable
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
    """Class that represents the basic property of an Observable. Concrete
    implementations are its subslasses."""
    _meanQ = None
    _nbins = None
    def __init__(self, filename, order):
        self.filename = filename
        if not order in ORDERS_QCD:
            if order in NUMS_QCD:
                order = NUMS_QCD[order]
            else:
                raise ValueError("Invalid value for order")
        self.order = order

    @property
    def meanQ(self):
        return self._meanQ

    @property
    def name(self):
        return osp.splitext(osp.basename(self.filename))[0]

    @property
    @fastcache.lru_cache()
    def sha1hash(self):
        with open(self.filename, 'rb') as f:
            return hashlib.sha1(f.read()).digest()

_selected_grid = None

class APPLGridObservable(Observable):
    """Class that represents an APPLGrid. """

    @property
    def nbins(self):
        """Number of bins in the APPLGrid. It will be loaded
         in memory the first time this property is quiried."""
        if self._nbins is not None:
            return self._nbins
        with self:
            nbins = applwrap.getnbins()
        self._nbins = nbins
        return nbins

    @property
    def meanQ(self):
        """A list containing the mean energy of
        the nonzero weights of each bin"""
        if self._meanQ is not None:
            return self._meanQ
        with self:
            meanQ = [applwrap.getobsq(self.order, i) for
                 i in range(self.nbins)]
        self._meanQ = meanQ
        return meanQ

    def __enter__(self):
        """Load observable file in memory, using `with obs`.

        Note: Except random bugs due to APPLGrid poor implementation
        when loading several observables in the same process. In particular,pdf
        convolutions of grids made with AMCFast will not work."""
        global _selected_grid
        if _selected_grid == self.filename:
            return
        if _selected_grid is not None:
            raise RuntimeError("Contrdicting observable scope. "
                               "Was %s and trying to enter %s" %
                               (_selected_grid, self.filename))

        with contextlib.ExitStack() as stack:
            if not logging.getLogger().isEnabledFor(logging.DEBUG):
                stack.enter_context(supress_stdout())
            applwrap.initobs(self.filename)
        _selected_grid = self.filename
    def __exit__(self, exc_type, exc_value, traceback):
        global _selected_grid
        _selected_grid = None

class PredictionObservable(Observable):
    """Class representing a prediction in the custom SMPDF format."""
    def __init__(self, filename):
        self.filename = filename
        with open(filename) as f:
            d = yaml.load(f)
        #TODO: All checking
        self._params = d

    def to_result(self, pdfset):
        """Convert a prediction for the specified `pdfset`
        into a `Result` instance."""
        if str(pdfset) not in self.pdf_predictions:
            raise ValueError("No predictions found for pdf %s" % pdfset)
        path = self.pdf_predictions[str(pdfset)]
        if not osp.isabs(path):
            path = osp.join(osp.dirname(self.filename), path)
        #TODO: Transpose all result dataframes
        datas = pd.DataFrame.from_csv(path, sep='\t', index_col=0).T
        return make_result(self, pdfset, datas)

    #Raise error different from AttributeError, so it's not silented by
    # __getattr__
    @property
    def nbins(self):
        try:
            return self._params['nbins']
        except KeyError:
            raise KeyError("Description file must have nbins")


    @property
    def meanQ(self):
        if isinstance(self.energy_scale, numbers.Number):
            return [self.energy_scale]*self.nbins
        else:
            return self.energy_scale

    def __getattr__(self, attr):
        try:
            return self._params[attr]
        except KeyError:
            raise AttributeError()

class PDFDoesNotExist(AttributeError): pass

_selected_pdf = None
_context_pdf = None
class PDF(TupleComp):
    """A class representig the metadata and content of an LHAPDF grid.
    The attributes of the `.info` file can be queried directly as attribute
    of PDF objects.

    Parameters
    ----------
    name : str
           The LHAPDF name of the set. It do methoes not need to be installed
           in the
           LHAPDF path at the time the constructor is called, but the methods
           that require reading the metadata will fail.
    label : str
           A label used for plotting (instrad the gris name).
    """
    def __init__(self, name, label=None):

        self.name = name
        if label is None:
            label = name
        self.label = label


    def get_key(self):
        """Return string to indentify this object in the database"""
        #Convert python2 unicode to string so no u'prefix' is printed
        return (str(self.name),)


    def __enter__(self):
        """Load PDF in memory."""
        global _selected_pdf
        global _context_pdf
        if _selected_pdf == str(self):
            _context_pdf = str(self)
            return
        if _context_pdf is not None and _context_pdf != str(self):
            raise RuntimeError("Contrdicting PDF scope. "
                               "Was %s and trying to enter %s" %
                               (_context_pdf, self))
        _selected_pdf = str(self)
        _context_pdf = str(self)
        applwrap.initpdf(self.name)

    def __exit__(self, exc_type, exc_value, traceback):
        global _context_pdf
        _context_pdf = None

    def __str__(self):
        return self.name

    def __repr__(self):
        return "<%s:%s>" % (self.__class__.__name__, self.name)

    @property
    def oqcd_str(self):
        """String corresponging to the QCD perturbative order, such as 'LO'
        or 'NL0'."""
        return ORDERS_QCD[self.OrderQCD]

    @property
    def collaboration(self):
        """Infer the collaboration from the name of the PDF"""
        return lhaindex.get_collaboration(self.name)

    @property
    def as_from_name(self):
        r"""Infer :math:`\alpha_S(M_Z)` from the name of the PDF. This is
        useful
        when willing to group by :math:`\alpha_S` value sets of different
        number of flavours."""
        return lhaindex.as_from_name(self.name)

    @property
    def mref(self):
        """String corresponding to the reference physical quantity used to fix
        the energy reference. ($M_Z$ for $N_f \geq 5$)."""
        return "M_%s" % M_REF[self.NumFlavors]

    @property
    def reps(self):
        """Returin an iterator over the replica indexes (zero indexed, where
        0 is the mean replica)."""
        return range(self.NumMembers)

    @property
    def q2min_rep0(self):
        """Retreive the min q2 value of repica zero. NNote that this will
        load the whole grid if not already in memory."""
        with self:
            res = applwrap.q2Min()
        return res

    def make_xgrid(self, xminlog=1e-5, xminlin=1e-1, xmax=1, nplog=50, nplin=50):
        """Provides the points in x to sample the PDF. `logspace` and `linspace`
        will be called with the respsctive parameters."""

        return np.append(np.logspace(np.log10(xminlog), np.log10(xminlin),
                                           num=nplog, endpoint=False),
                         np.linspace(xminlin, xmax, num=nplin, endpoint=False)
                        )

    def make_flavors(self, nf=3):
        return np.arange(-nf,nf+1)

    @fastcache.lru_cache(maxsize=128, unhashable='ignore')
    def grid_values(self, Q, xgrid=None, fl=None):

        if Q is None:
            Q = self.q2Min
        if xgrid is None:
            xgrid = self.make_xgrid()
        #Allow tuples that can be saved in cache
        elif isinstance(xgrid, tuple):
            xgrid = self.make_xgrid(*xgrid)

        if fl is None:
            fl = self.make_flavors()
        elif isinstance(fl, int):
            fl = self.make_flavors(fl)
        elif isinstance(fl, tuple):
            fl = self.make_flavors(*fl)

        with self:
            all_members = [[[applwrap.xfxQ(r, f, x, Q)
                             for x in xgrid]
                             for f in fl]
                             for r in range(len(self))]

            all_members = np.array(all_members)
            mean = all_members[0]
            replicas = all_members[1:]

        return mean, replicas

    def xfxQ(self, rep, fl, x, Q):
        with self:
            res = applwrap.xfxQ(rep, fl, x, Q)
        return res

    @property
    def infopath(self):
        return lhaindex.infofilename(self.name)

    @property
    @fastcache.lru_cache()
    def sha1hash(self):
        """Currently hashs the info file. This could not be sufficient for
        as proof since the grid inside could be modified."""
        with open(self.infopath, 'rb') as f:
            return hashlib.sha1(f.read()).digest()


    def __getattr__(self, name):
        #next is for pandas not to get confused
        if name.startswith('__') or name == 'next':
            raise AttributeError()
        try:
            return lhaindex.parse_info(self.name)[name]
        except KeyError:
            raise AttributeError()
        except IOError:
            raise PDFDoesNotExist(self.name)


    def __len__(self):
        return self.NumMembers


class Result():
    """A class representing a result of the computation of an observable for
    each member of a PDF set. `pd.DataFrame` will be called on `data`. The
    result must have the bins as the columns and each replica as the rows.
    Subclasses of `Result` provide specialized methods to compute uncertainty.

    #TODO: This docstring is **WRONG**, must traspose evey signle DataFrame!

    Parameters
    ----------
    obs :
        `Observable`


    pdf :
        `PDF`


    data :
        `DataFrame`-like

    """
    def __init__(self, obs, pdf, data):
        self.obs = obs
        self.pdf = pdf
        self._data = pd.DataFrame(data)

    @property
    def central_value(self):
        """Return a `Series` containing the central value for each bin."""
        return self._cv

    @property
    def _cv(self):
        return self._data[0]

    @property
    def _all_vals(self):
        return self._data.iloc[:,1:]

    @property
    def nrep(self):
        """Number of PDF members"""
        return self._all_vals.shape[1]

    @property
    def nbins(self):
        """Number of bins in the preduction."""
        return self._all_vals.shape[0]

    @property
    def meanQ(self):
        """Mean energy of the observable, for each bin"""
        return self.obs.meanQ

    def std_error(self, nsigma=1):
        raise NotImplementedError("No error computation implemented for this"
                                  "type of set")

    def sample_values(self, n):
        raise NotImplementedError("No sampling implemented for this"
                                  "type of set")

    def errorbar68(self):
        raise NotImplementedError("Error computation for PDF set"
                                   " not implmented")

    def std_interval(self, nsigma=1):
        std = self.std_error(nsigma)
        return pd.DataFrame({'min':self._cv - std,
                             'max':self._cv + std})

    def rel_std_interval(self, nsigma=1):
        std = self.std_error(nsigma)
        return pd.DataFrame({'min':-std,
                             'max':std})

    #TODO: Start from 1 by default?
    @property
    def binlabels(self):
        "Return the labels of the bins for this observable"
        return self._data.index


    def __getitem__(self, item):
        return self._data[item]

    def iterreplicas(self):
        """Iterate over all data, first being the central prediction"""
        return iter(self._data)

#==============================================================================
#     def __iter__(self):
#         """Give the prdictions for each bin"""
#         return self._all_vals.iterrows()
#
#     def __len__(self):
#         return len(self._all_vals)
#==============================================================================


    def _violin_data(self, rel_to=None):
        absdata = pd.concat(self.sample_values(10000),axis=1)
        if rel_to is None:
            rel_to = 1
        reldata = absdata.as_matrix().T/rel_to
        return reldata

    def violin_plot(self, data=None , **kwargs):
        if data is None:
            data = self._violin_data()

        myargs = {'label': str(self.pdf.label)}
        myargs.update(kwargs)
        return plotutils.violin_plot(data, **myargs)

    def sumbins(self, bins = None):
        sumobs = BaseObservable(self.obs.name + '[Sum]', self.obs.order)
        data = pd.DataFrame(self._data.sum(axis=0)).T
        return self.__class__(sumobs, self.pdf, data)

    @property
    def sha1hash(self):
        return hashlib.sha1(self.pdf.sha1hash + self.obs.sha1hash).digest()



class SymHessianResult(Result):
    """Result obtained from a symmetric Hessain PDF set"""
    
    def rescale_ci(self):
        if hasattr(self.pdf, "ErrorConfLevel"):
            return scipy.stats.norm.isf((1 - self.pdf.ErrorConfLevel)/2)
        else:
            return 1

    def std_error(self, nsigma=1):
        diffsq = (self._all_vals.subtract(self._cv, axis=0))**2
        return diffsq.sum(axis=1).apply(np.sqrt)*nsigma*self.rescale_ci()

    @property
    def errorbar68(self):
        """Compute the errorbars from the one sigma error"""
        return self.rel_std_interval()

    def sample_values(self, n):
        """Sample n random values from th resulting Gaussian distribution"""
        diffs = self._all_vals.subtract(self._cv, axis=0)
        for _ in range(n):
            weights = np.random.normal(scale=self.rescale_ci(), size=self.nrep)
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
    """Result obtained from an asymmetric Hessian PDF set"""

    def std_error(self, nsigma=1):
        m = self._all_vals.as_matrix()
        diffsq = (m[:, ::2] - m[:, 1::2])**2
        return np.sqrt(diffsq.sum(axis=1))/2.0*nsigma*self.rescale_ci()

    def sample_values(self, n):
        """Sample n random values from the resulting asymmetric
        distribution"""
        m = self._all_vals.as_matrix()
        plus = m[:, ::2]
        minus = m[:, 1::2]

        for _ in range(n):
            r = np.random.normal(scale=self.rescale_ci(), size=len(plus))
            error = (r >=0)*r*plus - (r < 0)*r*minus
            yield self._cv + error


class MCResult(Result):
    """Result obtained from a Monte Carlo PDF set"""
    def centered_interval(self, percent=68, addcentral=True):
        """Compute the 69% prediction gor each bin in the following way:
        Sort all results by the absolute value of the distance from the mean,
        and select """
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
        """Sample n random values from the results for the replicas"""
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



def convolve_one(pdf, observable, logger=None):
    import applwrap
    from smpdflib.core import PDF, APPLGridObservable #analysis:ignore
    res = {}
    import os
    logging.debug("Convolving in PID: %d" % os.getpid())
    logging.info("Convolving %s with %s" % (observable, pdf))
    with pdf, observable:
        for rep in pdf.reps:
            applwrap.pdfreplica(rep)
            res[rep] = np.array(applwrap.convolute(observable.order))
    return res


def make_convolution(pdf, observables):
    datas = defaultdict(lambda:OrderedDict())
    if not observables:
        return {}

    with(pdf):
        for obs in observables:
            with obs:
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
    """Convolve a set of pdf with a set of observables. Note that to get rid of
    issues arising from applgrid poor design, the multiprocessing start method
    must be 'spawn', ie:

    .. code:: python

        multiprocessing.set_start_method('spawn')

    Only once at the beginning of the program. This only works in Python 3.4+.
    """
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


    nprocesses = min((n_cores, len(to_compute)))
    if nprocesses:
        q = get_logging_queue()
        loglevel = logging.getLogger().level
        #http://stackoverflow.com/questions/30943161/multiprocessing-pool-with-maxtasksperchild-produces-equal-pids#30943161
        pool = multiprocessing.Pool(processes=nprocesses,
                                    maxtasksperchild=1,
                                    initializer=initlogging,
                                    initargs=(q, loglevel))
        results = pool.starmap(convolve_one, to_compute, chunksize=1)
        pool.close()
        for ((pdf, obs), result) in zip(to_compute, results):
            dataset[pdf][obs] = result
            if db is not None:
                key = make_key(pdf, obs)
                logging.debug("Appending result for %s to db" % key)
                db[key] = result

    return dataset


def convolve_or_load(pdfsets, observables, db=None):
    #results = []
    results = results_from_datas(get_dataset_parallel(pdfsets, observables, db))
    return results

def produce_results(pdfsets, observables, db=None):
    if isinstance(pdfsets, PDF) or isinstance(pdfsets, str):
        pdfsets = [pdfsets]
    pdfsets = [PDF(pdf) if isinstance(pdf,str) else pdf for pdf in pdfsets]

    if isinstance(observables, Observable):
        observables = [observables]
    predictions = [obs for obs in observables if
                   isinstance(obs, PredictionObservable)]

    applgrids = [obs for obs in observables if
                   isinstance(obs, APPLGridObservable)]


    results = (convolve_or_load(pdfsets, applgrids, db) +
               [pred.to_result(pdfset)
                for pdfset in pdfsets for pred in predictions])
    return results

def test_as_linearity(summed_table, diff_from_line = 0.25):
    group_by = ('Observable','NumFlavors', 'PDF_OrderQCD', 'Collaboration')
    for (process, nf, oqcd, col), curve_df in summed_table.groupby(group_by,
                                                                   sort=False):
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




def match_spec(corrlist, smpdf_spec):
    corr_obs = {item.obs:item for item in corrlist}
    result = {}
    for label in smpdf_spec:
        result[label] = [corr_obs[obs] for obs in smpdf_spec[label]]

    return result


def get_X(pdf, Q=None,  reshape=False, xgrid=None, fl=None):
    # Step 1: create pdf covmat
    if Q is None:
        Q = pdf.q2min_rep0
    logging.debug("Building PDF matrix at %f GeV:" % Q)
    mean, replicas = pdf.grid_values(Q, xgrid, fl)
    Xt = (replicas - mean)
    if reshape:
        Xt = Xt.reshape(Xt.shape[0], Xt.shape[1]*Xt.shape[2])
    return Xt.T

