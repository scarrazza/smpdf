# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 16:29:18 2015

@author: zah
"""
import json
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
                           )

from smpdflib.reducedset import get_smpdf_params

pdf = PDF("MC900_nlo")
obs = [make_observable(path, order='NLO') for path in itertools.chain(
       sorted(glob.glob("data/z/*.root")),
       sorted(glob.glob("data/w/*.root")),
       sorted(glob.glob("data/higgs/*.root")),
       sorted(glob.glob("data/ttbar/*.root")),
       )]



tolerance = 0.05

thresholds = [0.0, 0.25, 0.5, 0.75, 0.9, 0.99]

if __name__ == '__main__':
    neig = []
    logging.basicConfig(level=logging.INFO)
    db = shelve.open('db/db')

    results = produce_results([pdf], obs,  db=db)
    for t in thresholds:
        V, _ ,desc = get_smpdf_params(pdf, results, smpdf_tolerance=tolerance,
                          correlation_threshold=t, db=db)
        neig.append(V.shape[1])
        print("For thresholf %.2f we get:" %t)
        print(desc)
    with open("thresholdsladder.json", 'w') as f:
        json.dump([thresholds, neig], f)

    db.close()
