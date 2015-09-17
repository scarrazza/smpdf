# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 18:04:42 2015

@author: zah
"""
import shelve
from smpdflib.initialization import init_app
if __name__ == '__main__':
    try:
        init_app()
    except RuntimeError:
        pass

import numpy as np
import numpy.linalg as la
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib.ticker import MaxNLocator


from smpdflib.core import PDF, make_observable, produce_results

thresholds = [0, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99]

def to_percent(y, position):
    # Ignore the passed in position. This has the effect of scaling the default
    # tick locations.
    s = str(100 * y)

    # The percent symbol needs escaping in latex
    return '$' + s + r'\%$'

formatter = FuncFormatter(to_percent)

obs_name = 'data/z/z_13tev.root'

prior_name = "MC900_nnlo"

prefixes = ['z2tr'  + x for x in '0 10 25 50 75 90 99'.split()]

smpdf_names = [prefix + 'smpdf_' + prior_name for prefix in prefixes]

lincoef_paths = ['output/%s_lincomb.csv' % name for name in smpdf_names]

if __name__ == '__main__':
    prior = PDF(prior_name)
    smpdfs = [PDF(name) for name in smpdf_names]
    pdfs = [prior] + smpdfs
    obs = make_observable(obs_name, order='NLO')
    with shelve.open('db/db') as db:
        res_prior, *res_smpdfs = produce_results(pdfs, [obs], db=db)
    coefs = [pd.DataFrame.from_csv(path, sep='\t') for path in lincoef_paths]

    prior_std = res_prior.std_error()

    real_tols = [1 - res.std_error()/prior_std for res in res_smpdfs]

    prior_diffs = (res_prior._all_vals.T - res_prior.central_value).as_matrix().ravel()

    neig = [len(pdf) - 1 for pdf in smpdfs]

    rotated_tols = []
    for path in lincoef_paths:
        coefs = pd.DataFrame.from_csv(path, sep='\t')
        rotated_diffs = np.dot(prior_diffs, coefs)
        rotated_std = la.norm(rotated_diffs)
        tol = 1 - rotated_std/prior_std
        rotated_tols.append(tol)


    def _align():
        yield 'left'
        while True:
            yield 'right'
    align = _align()
    plt.plot(thresholds, real_tols, marker='s', label="Real T")
    plt.plot(thresholds, rotated_tols, marker='o', label="Linear T")
    plt.xticks(thresholds)
    for text, x, y in zip(neig, thresholds, real_tols):
        plt.annotate(text, (x, y), horizontalalignment=next(align), verticalalignment='bottom',)
    plt.gca().yaxis.set_major_formatter(formatter)

    plt.xlabel('t')
    plt.ylabel('T')

    plt.title(str(obs))

    color = next(plt.gca()._get_lines.color_cycle)
    plt.axhline(0.01, linestyle='--', label="Requested T", color=color)

    plt.gca().yaxis.set_major_locator(MaxNLocator(prune='lower'))

    plt.grid(axis='x')


    plt.legend()

    plt.savefig('tolerancesz2.pdf')


