# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 16:29:18 2015

@author: zah
"""
import shelve
import logging

from smpdflib.initialization import init_app
if __name__ == '__main__':
    init_app()

import matplotlib.pyplot as plt

from smpdflib.core import (PDF, make_observable, produce_results,
                           get_smpdf_lincomb)


pdf = PDF("1000rep")
obs = make_observable("data/higgs/ggh_13tev.root", order='NLO')

tolerance = 0.05

thresholds = [0.0, 0.25, 0.5, 0.75, 0.9]

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    db = shelve.open('db/db')
    results = produce_results([pdf], [obs],  db=db)
    for t in thresholds:
        V, desc = get_smpdf_lincomb(pdf, results, target_error=tolerance,
                          correlation_threshold=t)
        #Sprint(V.shape)
        print(desc)
    db.close()