# -*- coding: utf-8 -*-
"""
Created on Mon May  4 18:58:08 2015

@author: zah
"""
import os.path as  osp
import re

import pandas as pd

import smpdflib as lib


#TODO: Specify base
def save_violins(results, base_pdf, output_dir):
    for obs, fig in lib.compare_violins(results, base_pdf = results[0].pdf):
        filename = osp.join(output_dir, "figures", '%s.pdf' % obs)
        fig.savefig(filename)

def save_as(results, output_dir):
    if not isinstance(results, pd.DataFrame):
        results = lib.summed_results_table(results)
    for (process, nf), fig in lib.plot_alphaS(results):
        name = 'alpha_plot_%s(nf_%d).pdf'%(process,nf)
        name = re.sub(r'[^\w\.]', '_', name)
        fig.savefig(osp.join(output_dir, "figures", name))