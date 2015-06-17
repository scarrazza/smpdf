# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 15:11:04 2015

@author: zah
"""
import glob
import shelve

from smpdflib import core

pdf = core.PDF('1000rep')

observables = [core.make_observable(path, order = 1)
               for path in glob.glob("data/higgs/ggH_pt_13tev.root")]

db = shelve.open("db/db")

results = core.produce_results([pdf], observables, db = db)
results_table = core.results_table(results)

core.create_smpdf(pdf, results_table, '.', 'xpt', 1231, None, full_grid=False,
                  )

db.close()