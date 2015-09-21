# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 10:25:33 2015

@author: zah
"""
import os.path as osp
import logging
import hashlib
import tempfile
import uuid
import collections
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

def _mask_X(X, diffs, correlation_threshold=DEFAULT_CORRELATION_THRESHOLD):
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
                      correlation_threshold=DEFAULT_CORRELATION_THRESHOLD,
                      Rold = None):
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

    desc = []

    index = 0

    for result in pdf_results:
        obs_desc = {}
        desc.append({str(result.obs) : obs_desc})
        if result.pdf != pdf:
            raise ValueError("PDF results must be for %s" % pdf)
        for b in result.binlabels:
            Xreal = get_X(pdf, Q=result.meanQ[b], reshape=True)
            prediction = result._all_vals.ix[b]
            original_diffs = prediction - np.mean(prediction)
            if Rold is not None:
                X = np.dot(Xreal,Rold)
                rotated_diffs = np.dot(original_diffs, Rold)
            else:
                rotated_diffs = original_diffs
                X = Xreal

            eigs_for_bin = 0
            error_val = next(target_error)
            while _get_error(rotated_diffs, original_diffs) > error_val:
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
    lincomb = lincomb[:,:index]
    logging.info("Final linear combination has %d eigenvectors" %
                 lincomb.shape[1])


    return lincomb, norm, desc, Rold

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

def filter_results(real_results, prior_results, target_error):
    ...


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
                 Neig_total=200, full_grid=False, db =None,
                 correlation_threshold=DEFAULT_CORRELATION_THRESHOLD,
                 nonlinear_correction=True):

    lincomb, norm ,description, Rold = get_smpdf_lincomb(pdf, pdf_results,
                                               full_grid=full_grid,
                                               target_error=smpdf_tolerance,
                                               correlation_threshold=correlation_threshold)

    vec = lincomb/norm

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
            real_results = produce_results(temppdf, observables)
        logging.info("Real results obtained")





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


