# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 10:27:38 2015

@author: zah
"""
import numpy as np

DEFAULT_CORRELATION_THRESHOLD = 0.9


def corrcoeff(prediction, pdf_val):
    return (
            len(pdf_val)/(len(pdf_val)-1)*
            (np.mean(prediction*pdf_val) - np.mean(prediction)*np.mean(pdf_val))/
            (np.std(prediction,ddof=1)*np.std(pdf_val,ddof=1))
            )

def corrcoeffhessian(pdf_val):
    size = pdf_val.shape[0]
    cc = np.zeros(shape=(size,size))
    for i in range(size):
        for j in range(i, size):
            v1 = pdf_val[i,::2] - pdf_val[i,1::2]
            v2 = pdf_val[j,::2] - pdf_val[j,1::2]
            cc[i,j] = cc[j,i] = np.sum(v1*v2)/np.sum(v1**2)**0.5/np.sum(v2**2)**0.5
    return cc

def corrcoeffsymmhessian(central, pdf_val):
    size = pdf_val.shape[0]
    cc = np.zeros(shape=(size,size))
    for i in range(size):
        for j in range(i, size):
            v1 = pdf_val[i,:] - central[i]
            v2 = pdf_val[j,:] - central[j]
            cc[i,j] = cc[j,i] = np.sum(v1*v2)/np.sum(v1**2)**0.5/np.sum(v2**2)**0.5
    return cc

def bin_corrs_from_X(bin_val, X,
                     correlation_threshold=DEFAULT_CORRELATION_THRESHOLD):
    nxf, nrep = X.shape
    cc = np.zeros(shape=(nxf))
    #TODO: Optimize this
    for i in range(nxf):
        cc[i] = corrcoeff(bin_val, X[i,:])
    #Replace nan with zero when std is zero, ie when pdf value at x===0
    cc[np.isnan(cc)] = 0
    threshold = np.max(np.abs(cc))*correlation_threshold
    return cc, threshold

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
        #Without atleast_2d would return a scalar if M has only one element.
        corrmat = np.atleast_2d(np.corrcoef(M))
        title = "Observables correlation\n%s" % pdf.label
        filename = "%s" % pdf
        if base_pdf:
            corrmat -= base_corr
            title += "-%s" % base_pdf.label
            filename += "-%s" % base_pdf


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
        yield title, corrmat, labels, filename
