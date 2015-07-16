import os.path as osp
import shelve

import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot as plt

from smpdflib.core import make_observable, PDF, produce_results, results_table

from smpdflib.plots import plot_bindist



import smpdflib.core as lib
libpath = osp.dirname(lib.__file__)
stylefilename = osp.join(libpath, 'main.mplstyle')
plt.style.use(stylefilename)

db = shelve.open('db/db')


to_plot = {make_observable('data/applgrid/CMSWCHARM-WpCb-eta4.root', order=1)                      : 4,
           make_observable('data/applgrid/APPLgrid-LHCb-Z0-ee_arXiv1212.4260-eta34.root', order=1) : 8,
          }

pdfs = [PDF('MC900_nnlo', label="legendtext"),
        PDF('CMC100_nnlo', label="other"),
        PDF('MCH_nnlo_100', label="another")
       ]


for obs, bin in to_plot.items():
    results = produce_results(pdfs, [obs], db=db)

    obs_table = results_table(results)
    for (obs,b), fig in plot_bindist(obs_table, bin, base_pdf = pdfs[0]):
        ax = fig.axes[0]
        ax.set_xlabel("My X label")
        ax.set_ylabel("MY Y label")
        ax.set_title("My title")
        path = "%s_bin_%s.pdf" % (obs, b+1)

        fig.savefig(path)
        plt.close(fig)
