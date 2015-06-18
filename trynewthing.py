# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 15:11:04 2015

@author: zah
"""
#TODO: Convert this into proper unit test
import glob
import shelve
import logging

from smpdflib import core

logging.basicConfig(format='%(levelname)s: %(message)s')
logging.getLogger().setLevel(logging.INFO)

pdf = core.PDF('1000rep')

observables = [core.make_observable(path, order = 1)
               for path in glob.glob("data/higgs/*.root")]

db = shelve.open("db/db")

results = core.produce_results([pdf], observables, db = db)

core.create_smpdf(pdf, results, '.', 'hallN', full_grid=False,
                  db = db,  smpdf_tolerance=0.05)

db.close()