# -*- coding: utf-8 -*-
"""
Created on Fri May 15 15:35:45 2015

@author: zah
"""
from __future__ import division
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib.cbook import violin_stats
import matplotlib.patches
import matplotlib.mlab as mlab

from smpdflib import plotutils
from smpdflib.core import (aggregate_results, M_REF,
                           MCResult, get_X, bin_corrs_from_X,
                           observable_correlations, PDG_PARTONS)


colorlist = ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854',
                     '#ffd92f']

def compare_violins(results, base_pdf = None):
    if not isinstance(results, dict):
        combined = aggregate_results(results)
    else:
        combined = results
    for obs in combined:
        figure = plt.figure()
        norms = None
        handles = []
        plt.title(str(obs))
        colors  = plotutils.color_names_to_rgb(colorlist)
        alpha = 1
        base = combined[obs].get(base_pdf, None)
        results = sorted(combined[obs].values(), key = lambda x: x!=base)
        for result in results:
            if base is not None:
                cv = base.central_value.as_matrix()
                data = result._violin_data(rel_to=cv)
            else:
                data = data = result._violin_data()
            color = next(colors) + (alpha,)
            alpha /= 2
            plot, handle, norms = result.violin_plot(data, color=color,
                                              showextrema=False,
                                              normvalues=norms)
            handles.append(handle)
        plt.xlabel('bins')
        if base_pdf:
            plt.ylabel('Ratio to %s' % base_pdf.label)
        else:
            plt.ylabel("Observable value")
        plt.xticks(range(1,len(result.central_value) + 1))
        plt.legend(handles=handles)
        yield (obs,), figure

def compare_cis(results, base_pdf = None):
    if not isinstance(results, dict):
        combined = aggregate_results(results)
    else:
        combined = results
    for obs in combined:
        figure = plt.figure()
        colors  = plotutils.color_names_to_rgb(colorlist)
        plt.title(str(obs))
        base = combined[obs].get(base_pdf, None)
        results = sorted(combined[obs].values(), key = lambda x: x!=base)
        l = len(results)
        if l < 5:
            delta = iter(np.linspace(-0.05*l, 0.05*l, l))
        else:
            delta = iter(np.linspace(-0.25, 0.25, l))
        #x = iter(np.arange(1, len(results) + 1))
        for result in results:
            x = np.arange(1, result.nbins+1) + next(delta)
            data_cv = result.central_value.as_matrix().copy()
            data_ci = result.errorbar68.as_matrix().copy()
            if base is not None:
                base_cv = base.central_value.as_matrix()
                data_cv /= base_cv
                data_ci /= base_cv[:,np.newaxis]
            data_ci = np.abs(data_ci)
            color = next(colors)
            plt.errorbar(x, data_cv, yerr=data_ci.T,
                                     linestyle='none',
                                     color=color,
                                     label=result.pdf.label, elinewidth = 2,
                                     capsize=10)
            if isinstance(result, MCResult):
                data_std = result.rel_std_interval().as_matrix().copy()
                if base is not None:
                    data_std /= base_cv[:,np.newaxis]
                data_std = np.abs(data_std)
                plt.errorbar(x, data_cv, yerr=data_std.T,
                                         linestyle='none',
                                         color=color,
                                         elinewidth = 2,
                                         capsize=15)

        plt.xlabel('bins')
        if base_pdf:
            plt.ylabel('Ratio to %s' % base_pdf.label)
        else:
            plt.ylabel("Observable value")
        plt.xticks(range(1,len(result.central_value) + 1))
        plt.legend()
        plt.xlim(0.5, results[0].nbins + 0.5)
        plt.grid(axis='x')
        yield (obs,), figure

@plotutils.ax_or_gca
def plot_remarks(df, x, ax=None):
    have_remarks = df[df['Remarks'].apply(len) > 0]

    if len(have_remarks):
            ax.plot(have_remarks[x], have_remarks['CV'],
                 'ro', markersize = 20, fillstyle = 'none',
                 markeredgewidth = 5,
                 label="Problematic points")

def process_label(process, bin_):
    if bin_ == 'sum':
        return str(process)
    else:
        return "%s[bin:%d]"% (process, bin_)

#TODO: Abstract groupbyplots away? Tried, but seems too hard...
def plot_alphaS(results_table):
    x = 'alpha_sMref'
    df = results_table.sort(x)
    for (process, nf, bin_), process_df in df.groupby(['Observable',
                                                    'NumFlavors', 'Bin'],
                                                    sort = False):
        fig = plt.figure()


        for (oqcd,col), col_df in process_df.groupby(['PDF_OrderQCD',
                                                      'Collaboration']):
            label = "%s (%s)" % (col, oqcd)

            plt.errorbar(col_df[x], col_df['CV'],
                         yerr = np.array(col_df['Down68'],
                                         col_df['Up68']),
                        label = label, linestyle='-', marker = 's')


        plot_remarks(process_df, x)
        plt.xlabel(r'$\alpha_S(M_%s)$' % M_REF[nf])
        plt.ylabel(r'Value of observable')
        xran = plotutils.extend_range(process_df[x].min(),
                            process_df[x].max())
        plt.xlim(*xran)
        plt.legend()
        plt.title("%s $N_f=$%d" % (process_label(process, bin_), nf), y = 1.08)
        plt.tight_layout()
        yield (process, nf, bin_),fig

def plot_nf(results_table):
    x = 'NumFlavors'
    df = results_table.sort(x)
    for (process, bin_, oqcd), process_df in df.groupby(['Observable',
                                                    'Bin', 'PDF_OrderQCD'],
                                                    sort=False):
        fig = plt.figure()


        for (col, asn), pdf_df in process_df.groupby(['Collaboration',
                                                      'as_from_name'],
                                                      sort=False):
            label = "%s(as: %s)"%(col, asn)

            plt.errorbar(pdf_df[x], pdf_df['CV'],
                         yerr = np.array(pdf_df['Down68'],
                                         pdf_df['Up68']),
                        label = label, linestyle='-', marker = 's')


        plot_remarks(process_df, x)
        plt.xlabel(r'$N_f$')
        plt.ylabel(r'Value of observable')
        xran = plotutils.extend_range(process_df[x].min(),
                            process_df[x].max())
        plt.xlim(*xran)
        plt.xticks(process_df[x].unique())
        plt.legend()
        plt.title("%s PDFs, %s" % (oqcd, process_label(process, bin_)), y=1.08)
        plt.tight_layout()
        plt.grid(axis='x')
        yield  (process, bin_, oqcd),fig

def plot_asQ(pdfsets):
    df = pd.DataFrame([{'NumFlavors':pdf.NumFlavors, 'PDF':pdf} for
                        pdf in pdfsets])
    for nf, gdf in df.groupby(['NumFlavors'], sort=False):
        fig = plt.figure()
        for pdf in gdf.PDF:
            plt.plot(pdf.AlphaS_Qs, pdf.AlphaS_Vals, label=pdf.label)
            plt.ylabel(r'$\alpha_S$')
            plt.xlabel(r'Q(GeV)')
            plt.xscale('log')
            plt.title('$N_f=%d$' % nf)
            plt.legend()
        yield (nf,), fig

def plot_bindist(obs_table, b, base_pdf=None):

    def _kde_method(X, coords):
            kde = mlab.GaussianKDE(X, None)
            return kde.evaluate(coords)

    obs = obs_table.Observable.unique()
    if len(obs) != 1:
        raise ValueError("Must be only one observable")
    obs = obs[0]
    figure = plt.figure()
    plt.title("%s [Bin %d]" % (obs, b+1))
    colors  = plotutils.color_names_to_rgb(colorlist)
    alpha = 1
    if base_pdf:
        base = obs_table[obs_table.PDF.get_values()==base_pdf].Result[0]
    else:
        base = None
    results = obs_table.Result.unique(
    )
    for result in results:
        if base is not None:
            cv = base.central_value.as_matrix()
            data = result._violin_data(rel_to=cv)
        else:
            data = data = result._violin_data()

        if isinstance(data, list):
            stats = data[b]
        else:
            stats = violin_stats(data, _kde_method)[b]

        color = next(colors)
        alphacolor = color + (alpha,)
        plt.plot(stats['coords'], stats['vals'], color=color,label=result.pdf.label)
        plt.fill_between(stats['coords'], 0, stats['vals'], color=alphacolor,
                 )

        alpha /= 2


    plt.ylabel("Distribution")
    if base_pdf:
        plt.xlabel('Ratio to %s' % base_pdf.label)
    else:
        plt.xlabel("Observable value")
    plt.legend()
    yield (obs, b), figure


def plot_correlations(results):

    for result in results:
        pdf = result.pdf
        obs = result.obs

        Qs = obs.meanQ
        xgrid = pdf.make_xgrid()

        fl = pdf.make_flavors()

        figure, axarr = plt.subplots(len(fl), sharex=True,
                                     sharey=True,
                                     figsize=(8, len(fl)+3))

        for b in range(obs.nbins):
            Q = Qs[b]
            X = get_X(pdf, Q=Q, xgrid=xgrid, fl=fl, reshape=True)
            values, threshold = bin_corrs_from_X(result._all_vals.ix[b], X)
            ind = 0
            for f in range(len(fl)):
                step = len(xgrid)
                current_vals = values[ind:ind+step]
                ind+=step
                axarr[f].plot(xgrid, current_vals)
                axarr[f].set_ylim([-1,1])
                axarr[f].set_xscale('log')
                axarr[f].set_ylabel("$%s$"%PDG_PARTONS[fl[f]])

                axarr[f].axhline(threshold, c='r', ls='--')
                axarr[f].axhline(-threshold, c='r', ls='--')

        axarr[0].set_title(str(obs) + "\n")
        plt.xlabel("x")
        figure.subplots_adjust(hspace=0)
        plt.setp([a.get_xticklabels() for a in figure.axes[:-1]], visible=False)

        yield (obs,pdf), figure


def plot_observable_correlations(results_table, base_pdf=None):

    for title, corrmat, labels in observable_correlations(results_table, base_pdf):
        if base_pdf is not None:
            rlim = np.max(np.abs(corrmat))
        else:
            rlim = 1
        ranges = (-rlim, rlim)



        fig = plt.figure()
        plt.imshow(corrmat, cmap=plotutils.spectral_cm, vmin=ranges[0],
                   vmax=ranges[1], interpolation='none')
        plt.grid(False)
        if len(corrmat) < 20:
            ticks = np.arange(len(corrmat)), labels
            plt.xticks(*ticks, rotation=90)
            plt.yticks(*ticks)
        else:
            ...
            #empty = matplotlib.patches.Rectangle((0,0), 1, 1, fill=False,
            #                                     edgecolor='none',
            #                     visible=False)
            #plt.legend([empty]*len(corrmat), labels)

        plt.title(title)
        plt.tight_layout()
        plt.colorbar()

        yield (title.split()[-1],) , fig



