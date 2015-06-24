import os.path as osp
import shelve

from smpdflib.core import make_observable, PDF, produce_results, results_table

from smpdflib.plots import plot_bindist

import matplotlib.pyplot as plt

import smpdflib.core as lib


libpath = osp.dirname(lib.__file__)
stylefilename = osp.join(libpath, 'main.mplstyle')
plt.style.use(stylefilename)

db = shelve.open("db/db")
obs = make_observable("data/higgs/ggH_pt_13tev.root", order=1)
pdf = PDF("1000rep")
pdf2 = PDF('new-120-abs')
results = produce_results([pdf, pdf2], [obs], db=db)

obs_table = results_table(results)

for (obs,b), fig in plot_bindist(obs_table, 6, base_pdf = pdf):
    path = "%s_bin_%s.pdf" % (obs, b+1)
    fig.savefig(path)
    plt.close(fig)