# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 15:11:04 2015

@author: zah
"""
#TODO: Convert this into proper unit test
import glob
import shelve

from smpdflib import core

pdf = core.PDF('1000rep')

observables = [core.make_observable(path, order = 1)
               for path in glob.glob("data/higgs/*.root")]

db = shelve.open("db/db")

results = core.produce_results([pdf], observables, db = db)
results_table = core.results_table(results)

core.create_smpdf(pdf, results_table, '.', 'hall2', 1231, None, full_grid=False,
                  db = db)

db.close()