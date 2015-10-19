# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 15:41:26 2015

@author: zah
"""
import pandas as pd
import numpy as np

def save_html(df, path):
    import jinja2

    env = jinja2.Environment(loader = jinja2.PackageLoader('smpdflib',
                                                           'templates'))
    template = env.get_template('results.html.template')
    def remark_formatter(remarks):
        if not remarks:
            return ''
        else:
            return '<ul>%s</ul>' % '\n'.join('<li>%s</li>' %
                   jinja2.escape(remark) for remark in remarks)

    #http://stackoverflow.com/questions/26277757/pandas-to-html-truncates-string-contents
    with pd.option_context('display.max_colwidth', -1):
        table = df.to_html(
                             formatters={'Remarks':remark_formatter},
                             escape = False)
    result = template.render(table=table)
    with open(path, 'w') as f:
        f.write(result)

def break_bins(results):
    for result in results:
        for bin in range(result.nbins):
            name = str(result.obs)
            if result.nbins > 1:
                name += " Bin (%d)" % bin+1
            yield name, result._all_vals.iloc[bin,:]

def split_ranges(a,cond=None,*, filter_falses=False):
    """Split ``a`` so that each range has the same
    value for ``cond`` . If ``filter_falses`` is true, only the ranges
    for which the
    condition is true will be returned."""
    if cond is None:
        cond = a
    cond = cond.astype(bool)
    d = np.r_[False, np.diff(cond)]
    split_at = np.argwhere(d)
    splits = np.split(a, split_at)
    if filter_falses:
        #Evaluate condition at split points
        it = iter(cond[np.r_[0, np.ravel(split_at)]])
        return [s for s in splits if next(it)]
    else:
        return splits

