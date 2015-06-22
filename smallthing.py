import shelve

import numpy as np

from smpdflib.core import make_observable, PDF, produce_results

db = shelve.open("db/db")
obs = make_observable("data/higgs/ggH_pt_13tev.root", order=1)
pdf = PDF("1000rep")
result = produce_results([pdf], [obs], db=db)[0]
ratios = result._all_vals.mean(axis=1)/result[0]
print(ratios)

