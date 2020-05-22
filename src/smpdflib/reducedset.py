# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 10:25:33 2015

@author: zah


This module containts the functionality to compute reduced set using the
`mc2hessian` (Appendix of `1505.06736 <http://arxiv.org/abs/1505.06736>`_)
and `SMPDF` <paper> algorithms.
"""
import os.path as osp
import logging
import hashlib
import tempfile
import uuid
import collections
from collections import OrderedDict, namedtuple
import numbers
import itertools

import numpy as np
import numpy.linalg as la
import pandas as pd
import yaml

from smpdflib.core import get_X, produce_results, PDF
from smpdflib.lhio import hessian_from_lincomb
from smpdflib.corrutils import DEFAULT_CORRELATION_THRESHOLD, bin_corrs_from_X
import applwrap


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
    vec = V[:neig,:].T

    return vec

def merge_lincombs(lincomb1, lincomb2, desc1, desc2):
    """Merge `lincomb2` into `lincomb1` according to the specifications given
    by `desc1` and `desc2`.

    Parameters
    ----------
    lincomb1 : array nrep x neig1
           The original linear combination.
    lincomb1 : array nrep x neig2
           The refined linear combination we are merging. It contains
           eigenvectors to be inserted into the appropiate places.
    desc1: OrderedDict
           Description in the format returned by `get_smpdf_lincomb`: An
           ordered dictionary containing the number of eigenvectors needed to
           reproduce each bin of each observable (another ordered dict,
           {bin:neig}).
    desc2: OrderedDict
           The corresponding spec for lincomb2.
    """

    nrep = lincomb1.shape[0]
    if nrep != lincomb2.shape[0]:
        raise ValueError("Incompatible lincombs. "
                         "Need to have the same number replicas.")
    neig_tot = lincomb1.shape[1] + lincomb2.shape[1]
    final_lincomb = np.empty((nrep, neig_tot))
    last_neig2 = 0
    last_neig1 = 0
    total_index = 0
    desc = OrderedDict()
    for obs, obs_desc1 in desc1.items():
        obs_desc = OrderedDict()
        desc[obs] = obs_desc
        for b, neig1 in obs_desc1.items():
            #Catch both inner and outer KeyError
            try:
                neig2 = desc2[obs][b]
            except KeyError:
                neig2 = last_neig2


            nold =  neig1 - last_neig1
            final_lincomb[:, total_index:total_index+nold] = lincomb1[:,
                             last_neig1:last_neig1+nold]
            total_index += nold

            to_add = neig2 - last_neig2
            final_lincomb[:, total_index:total_index+to_add] = lincomb2[:,
                             last_neig2:last_neig2+to_add]
            total_index += to_add


            last_neig2 = neig2
            last_neig1 = neig1

            obs_desc[b] = total_index

    return final_lincomb, desc



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

def _mask_X(X, diffs, correlation_threshold=DEFAULT_CORRELATION_THRESHOLD):
     cc, threshold = bin_corrs_from_X(diffs, X, correlation_threshold=
                                                correlation_threshold)
     mask = np.abs(cc) > threshold
     Xm = X[mask]
     logging.debug("Masked shape is %s" % (Xm.shape,))
     return Xm

def _pop_eigenvector(X):
    """Extract the biggest eigenvector from X"""
    U,s,Vt = la.svd(X)
    Pt = Vt[:1,:]
    Rt = Vt[1:,:]
    return Pt.T, Rt.T

def _pdf_normalization(pdf):
    """Extract the quantity by which we have to divide the eigenvectors to
    get the correct errors, depending on the `ErrorType` of `pdf`."""
    nrep = len(pdf) - 1
    if pdf.ErrorType == 'replicas':
        norm = np.sqrt(nrep - 1)
    elif pdf.ErrorType in ('hessian', 'symmhessian'):
        norm = 1
    else:
        raise NotImplementedError("SMPDF is not implemented for this type of "
                                  "PDF error: %s" % pdf.ErrorType)
    return norm

SMPDFLincombResult = namedtuple('SMPDFLincombResult',
                                ('lincomb', 'norm', 'desc',
                                'errors', 'Rold'))


class TooMuchPrecision(Exception):
    def __init__(self, obs, b):
        super().__init__(("A SMPDF cannot be calculated with the requested "
        "precision for observable %s, bin %s. You need to increase the "
        "tolerance")%(obs,b))

def get_smpdf_lincomb(pdf, pdf_results,
                      target_error, full_grid = False,
                      correlation_threshold=DEFAULT_CORRELATION_THRESHOLD,
                      Rold = None):
    """Extract the linear combination that describes the linar part of
    the error of the given results with at least `target_error` precision`.
    See <paper> for details.
    `Rold` is returned so computation can be resumed iteratively (and then
    merged) with for example `merge_lincombs`."""
    #Estimator= norm**2(rotated)/norm**2(total) which is additive when adding
    #eigenvecotors
    #Error = (1 - sqrt(1-estimator))
    #TODO: Optimize by calculating estimator instead of error?
    #target_estimator = 1 - (1-target_error)**2

    nxf = len(pdf.make_xgrid())*len(pdf.make_flavors())
    nrep = len(pdf) - 1
    max_neig = np.min([nxf, nrep])
    #We must divide by norm since we are reproducing the covmat and not XX.T
    norm = _pdf_normalization(pdf)

    if isinstance(target_error, collections.Container):
        total_bins = sum(r.nbins for r in pdf_results)
        if len(target_error) != total_bins:
            raise ValueError("Incorrect target error specification")
        target_error = iter(target_error)
    elif isinstance(target_error, numbers.Real):
        target_error = itertools.repeat(target_error)
    elif not isinstance(target_error, collections.Iterator):
        raise ValueError("Target error not understood")


    lincomb = np.zeros(shape=(nrep,max_neig))

    desc = OrderedDict()
    errors = OrderedDict()

    index = 0

    for result in pdf_results:
        obs_desc = OrderedDict()
        obs_errors = OrderedDict()
        desc[str(result.obs)] = obs_desc
        errors[str(result.obs)] = obs_errors
        if result.pdf != pdf:
            raise ValueError("PDF results must be for %s" % pdf)
        for b in result.binlabels:
            Xreal = get_X(pdf, Q=result.meanQ[b], reshape=True)
            prediction = result._all_vals.loc[b]
            original_diffs = prediction - np.mean(prediction)
            if Rold is not None:
                X = np.dot(Xreal,Rold)
                rotated_diffs = np.dot(original_diffs, Rold)
            else:
                rotated_diffs = original_diffs
                X = Xreal

            eigs_for_bin = 0
            error_val = next(target_error)

            #Would be
            #while _get_error(rotated_diffs, original_diffs) > error_val
            #except that we want to capture current_error
            while True:
                current_error = _get_error(rotated_diffs, original_diffs)
                if current_error < error_val:
                    break
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
                index += 1
                if index == max_neig:
                    raise TooMuchPrecision(result.obs, b+1)
                eigs_for_bin += 1
            if eigs_for_bin:
                logging.info("Obtained %d eigenvector%s for observable %s, "
                             "bin %d" % (eigs_for_bin, 's'*(eigs_for_bin>1),
                                         result.obs, b+1))
            else:
                logging.debug("Observable %s, "
                             "bin %d is already well reproduced."
                             % (result.obs, b+1))
            obs_desc[int(b+1)] = index
            obs_errors[int(b+1)] = current_error

    #Prune extra zeros
    lincomb = lincomb[:,:index]
    logging.debug("Linear combination has %d eigenvectors" %
                 lincomb.shape[1])


    return SMPDFLincombResult(lincomb=lincomb, norm=norm, desc=desc,
                              errors=errors, Rold=Rold,
                             )

def complete_smpdf_description(desc, pdf ,pdf_results, full_grid,
                      target_error, total_neig ):

    input_hash = smpdf_input_hash(pdf, pdf_results, full_grid,
                                  target_error)

    desc = {'neig': desc,
           'input_hash' : input_hash,
           'target_tolarance' : target_error,
           'full_grid' : full_grid,
           'input_hash' : input_hash,
           'input_set' : str(pdf),
           'total_neig': total_neig,
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

def mc2h_input_hash(pdf, Q, Neig):
    hashstr = b''.join([pdf.sha1hash, float(Q).hex().encode(),
                        hex(Neig).encode()])
    return hashlib.sha1(hashstr).hexdigest()

def create_mc2hessian(pdf, Q, Neig, output_dir, name=None, db=None,
        photon=False):
    X = get_X(pdf, Q, reshape=True, photon=photon)
    vec = compress_X(X, Neig)
    norm = _pdf_normalization(pdf)
    description = {'input_hash': mc2h_input_hash(pdf,Q,Neig)}
    save_lincomb(vec, norm, description=description,
                 output_dir=output_dir, name=name)

    return hessian_from_lincomb(pdf, vec/norm, folder=output_dir,
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

def get_smpdf_params(pdf, pdf_results, smpdf_tolerance, full_grid=False,
                    db=None,
                    correlation_threshold=DEFAULT_CORRELATION_THRESHOLD,
                    nonlinear_correction=True):


    first_res = get_smpdf_lincomb(pdf, pdf_results,
                                  full_grid=full_grid,
                                  target_error=smpdf_tolerance,
                                  correlation_threshold=correlation_threshold)
    norm = first_res.norm
    lincomb = first_res.lincomb
    description = first_res.desc
    vec = first_res.lincomb/norm

    if nonlinear_correction:
        logging.info("Estimating nonlinear correction")
        with tempfile.TemporaryDirectory() as td:
            applwrap.setlhapdfpath(td)
            tempname = str(uuid.uuid1())
            logging.info("Creating temporary PDF %s" % tempname)
            hessian_from_lincomb(pdf, vec, folder=td,
                         set_name= tempname, db=db)
            observables = [r.obs for r in pdf_results]
            temppdf = PDF(tempname)
            temppdf.infopath
            real_results    = produce_results(temppdf, observables)
            logging.info("Real results obtained")

        results_to_refine = []
        newtols = []
        for smpdf_res, prior_res in zip(real_results, pdf_results):
            real_error = 1 - smpdf_res.std_error()/prior_res.std_error()
            #pandas indexing is broken, so have to call as_matrix....
            bad_bins = np.array(real_error > smpdf_tolerance, copy=False)
            if bad_bins.any():
                lin_errors = list(first_res.errors[str(smpdf_res.obs)].values())

                newtol = smpdf_tolerance - (real_error[bad_bins] -
                                             np.array(lin_errors)[bad_bins])

                impossible = np.argwhere(newtol < 0)
                if len(impossible):
                    raise TooMuchPrecision(smpdf_res.obs,
                                           impossible[0] + 1)
                newtols +=  list(newtol)
                logging.debug("New tolerances for observable %s: %s" %
                              (prior_res.obs, newtol))
                #Create result with the same type as prior, and only
                #bad_bins.
                newres = type(prior_res)(prior_res.obs, prior_res.pdf,
                                         prior_res._data.ix[bad_bins])

                results_to_refine.append(newres)
        if results_to_refine:
            logging.info("Calculating eigenvectors to refine")
            ref_res = get_smpdf_lincomb(pdf, results_to_refine,
                              full_grid=full_grid,
                              target_error=newtols,
                              correlation_threshold=correlation_threshold,
                              Rold=first_res.Rold)

            lincomb, description = merge_lincombs(first_res.lincomb,
                                                  ref_res.lincomb,
                                                  description,
                                                  ref_res.desc)

        else:
            logging.info("All results are within tolerance")

    return lincomb, norm, description



def create_smpdf(pdf, pdf_results, output_dir, name,
                 smpdf_tolerance,
                 full_grid=False, db=None,
                 correlation_threshold=DEFAULT_CORRELATION_THRESHOLD,
                 nonlinear_correction=True):

    lincomb, norm, description = get_smpdf_params(pdf, pdf_results,
                                     smpdf_tolerance,
                                     full_grid=full_grid,
                                     db=db,
                                     correlation_threshold=correlation_threshold,
                                     nonlinear_correction=nonlinear_correction)


    vec = lincomb/norm

    description = complete_smpdf_description(description, pdf, pdf_results,
                                             full_grid=full_grid,
                                             target_error=smpdf_tolerance,
                                             total_neig=lincomb.shape[1])
    #We have do do this because LHAPDF seems to not parse complex structures
    parsed_desc = {'smpdf_description':yaml.dump(description,
                                                 default_flow_style=False)}

    save_lincomb(lincomb, norm, description, output_dir, name)

    with open(osp.join(output_dir, name + '_description.yaml'), 'w') as f:
        yaml.dump(description, f, default_flow_style=False)

    logging.info("Final linear combination has %d eigenvectors" %
                 lincomb.shape[1])



    return hessian_from_lincomb(pdf, vec, folder=output_dir,
                         set_name= name, db=db, extra_fields=parsed_desc)
