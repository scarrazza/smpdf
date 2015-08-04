# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 16:29:18 2015

@author: zah
"""
import shelve
import logging
import glob
import itertools
import contextlib

from smpdflib.initialization import init_app
if __name__ == '__main__':
    with contextlib.suppress(RuntimeError):
        init_app()

import matplotlib.pyplot as plt

from smpdflib.core import (PDF, make_observable, produce_results,
                           get_smpdf_lincomb)


pdf = PDF("MC900_nnlo")
obs = [make_observable(path, order='NLO') for path in itertools.chain(
       glob.glob("data/z/*.root"),
       glob.glob("data/w/*.root"),
       glob.glob("data/higgs/*.root"),
       glob.glob("data/ttbar/*.root"),
       )]



tolerance = 0.05

thresholds = [0.0, 0.25, 0.5, 0.75, 0.9, 0.99]

if __name__ == '__main__':
    logging.basicConfig(level=logging.WARNING)
    db = shelve.open('db/db')

    try:
        results = produce_results([pdf], obs,  db=db)
    finally:
        db.close()
    for t in thresholds:
        V, desc = get_smpdf_lincomb(pdf, results, target_error=tolerance,
                          correlation_threshold=t)
        #Sprint(V.shape)
        print("For thresholf %.2f we get:" %t)
        print(desc)
