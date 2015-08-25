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
import numpy.linalg as la
import pandas as pd
import yaml
import scipy.stats
import fastcache
from pandas.stats import ols

from smpdflib import lhaindex
from smpdflib import plotutils
from smpdflib.loggingutils import supress_stdout, initlogging, get_logging_queue
from smpdflib.lhio import hessian_from_lincomb

import applwrap

_DEFAULT_CORRELATION_THRESHOLD = 0.9

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
        if getattr(self, 'CUSTOM_XFXQ', False):
            applwrap.set_custom_xfxQ(self.xfxQ)
        else:
            applwrap.set_custom_xfxQ(None)
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

    @fastcache.lru_cache(maxsize=2**16, unhashable='ignore')
    def grid_values(self, Q, xgrid=None, fl=None):

        if Q is None:
            Q = self.q2Min
        if xgrid is None:
            xgrid = self.make_xgrid()
        #Allow tuples that can be saved in cache
        elif isinstance(xgrid, tuple):
            xgrid = self.make_xgrid(*xgrid)
        elif isinstance(xgrid, numbers.Number):
            xgrid = [xgrid]

        if fl is None:
            fl = self.make_flavors()
        elif isinstance(fl, int):
            fl = [fl]
        elif isinstance(fl, tuple):
            fl = self.make_flavors(*fl)

        with self:
            all_members = [[[self.xfxQ(r, f, x, Q)
                             for x in xgrid]
                             for f in fl]
                             for r in range(len(self))]

            all_members = np.array(all_members)
            mean = all_members[0]
            replicas = all_members[1:]

        return mean, replicas

    @fastcache.lru_cache()
    def grid_values_from_tuple(self, fxQtuple):
        print("Len of tuple is")
        print(len(fxQtuple))
        replist = []
        for r in range(len(self)):
            l = []
            for val in fxQtuple:
                l.append(self.xfxQ(r, *val))
            replist.append(l)
        return pd.DataFrame(replist, columns = fxQtuple)

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

class LincombPDF(PDF):
    CUSTOM_XFXQ = True
    def __init__(self, base, lincomb):
        self.base  = base
        self.lincomb = lincomb
        self._grid_points = []
        super().__init__(self.base.name)


    def sha1hash(self):
        raise RuntimeError("Lincomb PDFs should not be hashed")

    def __getattr__(self, name):
        #This is needed for pickle to work
        if name == 'base':
            raise AttributeError(name)
        return getattr(self.base, name)

    @fastcache.lru_cache()
    def evaluate_points(self):
        vals = self.base.grid_values_from_tuple(tuple(self._grid_points))
        X = vals.ix[1:] - vals.ix[0]
        return vals.ix[0] + pd.DataFrame(np.dot(self.lincomb.T,X),
                                         columns=X.columns)

    def xfxQ(self, rep, fl, x, Q):
        if rep == 0:
            self._grid_points.append((fl, x, Q))
            return self.base.xfxQ(rep, fl, x, Q)

        return self.evaluate_points()[(fl, x, Q)][rep - 1]

    @property
    def NumMembers(self):
        return self.lincomb.shape[1] + 1

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

    def std_error(self, nsigma=1):
        diffsq = (self._all_vals.subtract(self._cv, axis=0))**2
        return diffsq.sum(axis=1).apply(np.sqrt)*nsigma

    @property
    def errorbar68(self):
        """Compute the errorbars from the one sigma error"""
        return self.rel_std_interval()

    def sample_values(self, n):
        """Sample n random values from th resulting Gaussian distribution"""
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
    """Result obtained from an asymmetric Hessian PDF set"""

    def std_error(self, nsigma=1):
        m = self._all_vals.as_matrix()
        diffsq = (m[:, ::2] - m[:, 1::2])**2
        return np.sqrt(diffsq.sum(axis=1))/2.0*nsigma

    def sample_values(self, n):
        """Sample n random values from the resulting asymmetric
        distribution"""
        m = self._all_vals.as_matrix()
        plus = m[:, ::2]
        minus = m[:, 1::2]

        for _ in range(n):
            r = np.random.normal(size=len(plus))
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
    logging.info("Convloving %s with %s" % (observable, pdf))
    with pdf, observable:
        for rep in pdf.reps:
            applwrap.pdfreplica(rep)
            res[rep] = np.array(applwrap.convolute(observable.order))
            logging.debug("Did bin in PID: %d" % os.getpid())
    logging.debug("Finished with PID: %d" % os.getpid())

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
    if isinstance(pdfsets, PDF):
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

def corrcoeff(prediction, pdf_val):
    return (
            len(pdf_val)/(len(pdf_val)-1)*
            (np.mean(prediction*pdf_val) - np.mean(prediction)*np.mean(pdf_val))/
            (np.std(prediction,ddof=1)*np.std(pdf_val,ddof=1))
            )


def bin_corrs_from_X(bin_val, X,
                     correlation_threshold=_DEFAULT_CORRELATION_THRESHOLD):
    nxf, nrep = X.shape
    cc = np.zeros(shape=(nxf))
    #TODO: Optimize this
    for i in range(nxf):
        cc[i] = corrcoeff(bin_val, X[i,:])
    #Replace nan with zero when std is zero, ie when pdf value at x===0
    cc[np.isnan(cc)] = 0
    threshold = np.max(np.abs(cc))*correlation_threshold
    return cc, threshold


def match_spec(corrlist, smpdf_spec):
    corr_obs = {item.obs:item for item in corrlist}
    result = {}
    for label in smpdf_spec:
        result[label] = [corr_obs[obs] for obs in smpdf_spec[label]]

    return result

@fastcache.lru_cache(2**16)
def get_X(pdf, Q=None, xgrid=None, fl=None, *, reshape=False):
    # Step 1: create pdf covmat
    if Q is None:
        Q = pdf.q2min_rep0
    logging.debug("Building PDF matrix at %f GeV:" % Q)
    mean, replicas = pdf.grid_values(Q, xgrid, fl)
    Xt = (replicas - mean)
    if reshape:
        Xt = Xt.reshape(Xt.shape[0], Xt.shape[1]*Xt.shape[2])
    return Xt.T

def decompose_eigenvectors(X, predictions, target_estimator):
    target_value = target_estimator

    U,s,Vt = la.svd(X)

    newrot = np.dot(Vt, predictions)
    total = np.dot(predictions, predictions)
    s = 0
    logging.debug("Target value: %.4f" % target_value)
    for i in range(len(newrot)):
        s += newrot[i]**2
        value = s/total
        logging.debug("Added new eigenvector. Value: %.4f" % value)

        if value >= target_value:
            neig = i + 1
            break
    else: #for .. else is no break
        neig = len(newrot)

    Pt = Vt[:neig,:]
    #Wt = np.zeros_like(Vt)
    Rt = Vt[neig:,:]
    return Pt.T, Rt.T

def compress_X(X, neig):
    U, s, V = np.linalg.svd(X, full_matrices=False)
    norm = np.sqrt(X.shape[1] - 1)
    sn = s/norm
    u = U[:,:neig]
    vec = V[:neig,:].T/norm
    cov = np.dot(u, np.dot(np.diag(sn[:neig]**2), u.T))

    return vec, cov


def _get_error(rotated_diffs, original_diffs):
    """Rotated diffs are the residual differences, that are not reproduced in
    the hessian. They are to be minimized. It holds that
    `norm(rotated_diffs)**2+norm(<reproduced>)**2==norm(original_diffs)**2`
    """

    rotsqnorm = np.dot(rotated_diffs, rotated_diffs)
    origsqnorm = np.dot(original_diffs, original_diffs)
    estimator = rotsqnorm/origsqnorm
    error = 1 - np.sqrt(1 - estimator)

    return error

def _mask_X(X, diffs, correlation_threshold=_DEFAULT_CORRELATION_THRESHOLD):
     cc, threshold = bin_corrs_from_X(diffs, X, correlation_threshold=
                                                correlation_threshold)
     mask = np.abs(cc) > threshold
     Xm = X[mask]
     logging.debug("Masked shape is %s" % (Xm.shape,))
     return Xm

def _pop_eigenvector(X):
    U,s,Vt = la.svd(X)
    Pt = Vt[:1,:]
    Rt = Vt[1:,:]
    return Pt.T, Rt.T

def _pdf_normalization(pdf):
    nrep = len(pdf) - 1
    if pdf.ErrorType == 'replicas':
        norm = np.sqrt(nrep - 1)
    elif pdf.ErrorType in ('hessian', 'symmhessian'):
        norm = 1
    else:
        raise NotImplementedError("SMPDF is not implemented for this type of "
                                  "PDF error: %s" % pdf.ErrorType)
    return norm

def get_smpdf_lincomb(pdf, pdf_results, full_grid = False,
                      target_error = 0.1,
                      correlation_threshold=_DEFAULT_CORRELATION_THRESHOLD,
                      nonlinear_estimate = True
                      ):
    #Estimator= norm**2(rotated)/norm**2(total) which is additive when adding
    #eigenvecotors
    #Error = (1 - sqrt(1-estimator))
    #TODO: Optimize by calculating estimator instead of error?
    #target_estimator = 1 - (1-target_error)**2
    Rold = None


    nxf = len(pdf.make_xgrid())*len(pdf.make_flavors())
    nrep = len(pdf) - 1
    max_neig = np.min([nxf, nrep])
    #We must divide by norm since we are reproducing the covmat and not XX.T
    norm = _pdf_normalization(pdf)

    lincomb = np.zeros(shape=(nrep,max_neig))
    all_observables = [result.obs for result in pdf_results]
    desc = []

    index = 0

    for result in pdf_results:
        obs_desc = {}
        desc.append({str(result.obs) : obs_desc})
        if result.pdf != pdf:
            raise ValueError("PDF results must be for %s" % pdf)
        for b in range(result.nbins):
            Xreal = get_X(pdf, Q=result.meanQ[b], reshape=True)
            prediction = result._all_vals.iloc[b,:]
            original_diffs = prediction - np.mean(prediction)
            if Rold is not None:
                X = np.dot(Xreal,Rold)
                rotated_diffs = np.dot(original_diffs, Rold)
            else:
                rotated_diffs = original_diffs
                X = Xreal

            eigs_for_bin = 0
            while _get_error(rotated_diffs, original_diffs) > target_error:
                X = _mask_X(X, rotated_diffs, correlation_threshold=
                                              correlation_threshold)
                P, R = _pop_eigenvector(X)
                if Rold is not None:
                    P = np.dot(Rold, P)
                    R = np.dot(Rold, R)
                Rold = R

                rotated_diffs = np.dot(original_diffs, Rold)
                X = np.dot(Xreal,Rold)
                lincomb[:,index:index+1] = P
                if nonlinear_estimate:
                    logging.info("Start wasting time")
                    lpdf = LincombPDF(pdf, P/norm)
                    produce_results(lpdf, all_observables)
                    logging.info("Stop wasting time")
                index += 1
                eigs_for_bin += 1
            if eigs_for_bin:
                logging.info("Obtained %d eigenvector%s for observable %s, "
                             "bin %d" % (eigs_for_bin, 's'*(eigs_for_bin>1),
                                         result.obs, b+1))
            else:
                logging.debug("Observable %s, "
                             "bin %d is already well reproduced."
                             % (result.obs, b+1))
            obs_desc[b+1] = index
    lincomb = lincomb[:,:index]
    logging.info("Final linear combination has %d eigenvectors" %
                 lincomb.shape[1])


    return lincomb, norm, desc

def complete_smpdf_description(desc, pdf ,pdf_results, full_grid,
                      target_error ):

    input_hash = smpdf_input_hash(pdf, pdf_results, full_grid,
                                  target_error)

    desc = {'neig': desc,
           'input_hash' : input_hash,
           'target_tolarance' : target_error,
           'full_grid' : full_grid,
           'input_hash' : input_hash,
           'input_set' : str(pdf),
           }

    return desc


#TODO: Add smpdf version info here
def smpdf_input_hash(pdf, pdf_results, full_grid,
                      target_error):

    hashstr = b''.join(r.obs.sha1hash for r in pdf_results)
    hashstr += target_error.hex().encode()
    hashstr += bytes(full_grid)
    input_hash = hashlib.sha1(hashstr).hexdigest()
    return input_hash


def create_mc2hessian(pdf, Q, Neig, output_dir, name=None, db=None):
    X = get_X(pdf, Q, reshape=True)
    vec, _ = compress_X(X, Neig)
    return hessian_from_lincomb(pdf, vec, folder=output_dir,
                         set_name= name, db=db)


def save_lincomb(lincomb, norm, description, output_dir, name):

    nrep, neig = lincomb.shape
    columns = [description['input_hash'][:8]+'_%d'%i for i in range(1,neig+1)]
    rows = range(1, nrep+1)

    frame = pd.DataFrame(lincomb, columns=columns, index=rows)
    direct = frame/norm
    dirname = name + "_lincomb.csv"
    direct.to_csv(osp.join(output_dir, dirname), sep='\t', float_format='%e')

    invname= name + "_lincomb_inverse.csv"
    inverse = frame.T*norm

    inverse.to_csv(osp.join(output_dir, invname), sep='\t', float_format='%e')


def create_smpdf(pdf, pdf_results, output_dir, name,  smpdf_tolerance=0.05,
                 Neig_total = 200,
                 full_grid=False, db = None,
                 correlation_threshold= _DEFAULT_CORRELATION_THRESHOLD):

    lincomb, norm ,description = get_smpdf_lincomb(pdf, pdf_results,
                                               full_grid=full_grid,
                                               target_error=smpdf_tolerance,
                                               correlation_threshold=correlation_threshold)

    vec = lincomb/norm

    description = complete_smpdf_description(description, pdf, pdf_results,
                                             full_grid=full_grid,
                                             target_error=smpdf_tolerance)
    #We have do do this because LHAPDF seems to not parse complex structures
    parsed_desc = {'smpdf_description':yaml.dump(description,
                                                 default_flow_style=False)}

    save_lincomb(lincomb, norm, description, output_dir, name)

    with open(osp.join(output_dir, name + '_description.yaml'), 'w') as f:
        yaml.dump(description, f, default_flow_style=False)


    return hessian_from_lincomb(pdf, vec, folder=output_dir,
                         set_name= name, db=db, extra_fields=parsed_desc)

def observable_correlations(results_table, base_pdf=None):

    if base_pdf is not None:
        base_results = results_table[results_table.PDF == base_pdf].Result.unique()
        M = np.concatenate([result._all_vals.as_matrix() for
                            result in base_results])
        base_corr = np.corrcoef(M)


    for pdf, pdf_table in results_table.groupby('PDF', sort=False):
        if base_pdf == pdf:
            continue
        results = pdf_table.Result.unique()
        M = np.concatenate([result._all_vals.as_matrix() for result in results])
        corrmat = np.corrcoef(M)
        title = "Observables correlation\n%s" % pdf
        if base_pdf:
            corrmat -= base_corr
            title += "-%s" % base_pdf


        def make_labels(results):
            for result in results:
                obs = result.obs
                obslabel = str(obs)
                if len(obslabel) > 10:
                    obslabel = obslabel[:10] + '...'
                if result.nbins == 1:
                    yield obslabel
                else:
                    for b in range(result.nbins):
                        yield "%s (Bin %d)" % (obslabel, b+1)

        labels = list(make_labels(results))
        yield title, corrmat, labels
