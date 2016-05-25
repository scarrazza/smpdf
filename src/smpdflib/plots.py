# -*- coding: utf-8 -*-
"""
Created on Fri May 15 15:35:45 2015

@author: zah
"""
#TODO: Enable this ASAP
#from __future__ import generator_stop
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib.cbook import violin_stats
import matplotlib.mlab as mlab
from matplotlib.ticker import MaxNLocator, LinearLocator, FixedFormatter


from smpdflib import plotutils
from smpdflib.core import (aggregate_results, M_REF, MCResult,
                           get_X, PDG_PARTONS, make_pdf_results, make_pdfcorr_results, PDF)

from smpdflib.corrutils import (bin_corrs_from_X, observable_correlations,
                                corrcoeff,)

from smpdflib.utils import split_ranges


colorlist = ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854',
                     '#ffd92f']

def plot_pdfs(pdfsets, Q, base_pdf = None, flavors=None, photon=False):

    if flavors is None:
        flavors = PDF.make_flavors(photon=photon)

    xgrid = PDF.make_xgrid()

    #Squeeze=False is for the case where we have one flavour
    w,h = plt.rcParams["figure.figsize"]

    resultlists = []
    for pdf in pdfsets:
        r = make_pdf_results(pdf, Q, flavors=flavors, xgrid=xgrid)
        if pdf == base_pdf:
            base_results = iter(r)
        resultlists.append(r)



    for fl, flresults in zip(flavors, zip(*resultlists)):
        fig, ax = plt.subplots()
        hatchiter = plotutils.hatch_iter()
        if base_pdf:
            base_result = next(base_results)

        view_min = float("inf")
        view_max = -float("inf")

        labels = []
        handles = []

        for result in flresults:
            #Without copy the /= operator affects the underlying data.
            central = result.central_value.as_matrix().copy()
            if isinstance(result, MCResult):
                stderror = result.rel_std_interval().as_matrix().copy()
            error = result.errorbar68.as_matrix().copy()
            if base_pdf:
                base_cv = base_result.central_value.as_matrix()
                central /= base_cv
                error /=  base_cv[:,np.newaxis]
                if isinstance(result, MCResult):
                    stderror /=  base_cv[:,np.newaxis]

            line, = ax.plot(xgrid, central)

            color = line.get_color()
            alpha = 0.2
            hatch = next(hatchiter)
            band_up = central + error[:,1]
            band_down = central+ error[:,0]
            ax.fill_between(xgrid, band_up,
                            band_down, color=color,
                            alpha=alpha, hatch=hatch, zorder=-1)
            if isinstance(result, MCResult):
                outer = True
                label = "%s ($68\%%$ c.l.+$1\sigma$)" % result.pdf.label
                ax.plot(xgrid, central + stderror[:,1], color=color,
                    linestyle="--")
                ax.plot(xgrid, central + stderror[:,0], color=color,
                    linestyle="--")
            else:
                outer = False
                label = "%s ($68\%%$ c.l.)" % result.pdf.label

            handle = plotutils.HandlerSpec(color=color, alpha=alpha,
                                           hatch=hatch, outer=outer)
            handles.append(handle)
            labels.append(label)

            #Adjust view to middle region
            lx = np.log10(xgrid)
            extrema = lx[[0,-1]]
            middle = np.mean(extrema)
            diff = np.diff(extrema)
            mask =  np.abs(lx - middle) < diff/4
            ext_max = np.percentile((central + 2*error[:,0])[mask], 90)

            ext_min = np.percentile((central + 2*error[:,1])[mask], 10)
            if ext_max > view_max:
                view_max = ext_max
            if ext_min < view_min:
                view_min = ext_min

        amin, amax = ax.get_ylim()
        ax.set_ylim(max(view_min, amin), min(view_max, amax))

        ax.set_xscale('log')

        ax.set_title("$%s$ at $Q=%.2f$ GeV" % (PDG_PARTONS[fl], Q))

        ax.legend(handles, labels, handler_map={plotutils.HandlerSpec:
                                                plotutils.ComposedHandler()})

        if base_pdf:
            ax.set_ylabel('Ratio to %s' % base_pdf.label)
        else:
            ax.set_ylabel("$x%s(x)$" % PDG_PARTONS[fl])
        ax.set_xlabel("$x$")

        fig.tight_layout()

        yield (fl,),fig

def plot_pdfcorr(pdfsets, Q, base_pdf=None, flavors=None, photon=False):

    if flavors is None:
        flavors = PDF.make_flavors(photon=photon)

    xgrid = PDF.make_xgrid()

    resultlists = []
    for pdf in pdfsets:
        r = make_pdfcorr_results(pdf, Q, flavors=flavors, xgrid=xgrid)
        if pdf == base_pdf:
            base_results = r
        resultlists.append(r)

    for res in resultlists:
        fig, frame = plt.subplots()
        corrmat = res.correlations()
        if base_pdf: corrmat -= base_results.correlations()

        plt.imshow(corrmat, cmap=plotutils.spectral_cm, vmin=-1, vmax=1)
        plt.grid(False)
        if base_pdf:
            plt.title("Correlations at $Q=%.2f$ GeV for\n%s-%s" %
                      (Q, res.pdf.label, base_pdf.label), fontsize=15)
        else:
            plt.title("Correlations at $Q=%.2f$ GeV for\n%s" % (Q, res.pdf.label),
                      fontsize=15)

        ylabels = [ '\n\n' + r'$%s$' % PDG_PARTONS[fl] for fl in flavors]
        xlabels = [ '\t' + r'$%s$' % PDG_PARTONS[fl] for fl in flavors]
        frame = plt.gca()
        frame.axes.get_yaxis().set_major_locator(LinearLocator(len(flavors)+1))
        frame.get_yaxis().set_major_formatter(FixedFormatter(ylabels))
        frame.axes.get_xaxis().set_major_locator(LinearLocator(len(flavors)+1))
        frame.get_xaxis().set_major_formatter(FixedFormatter(xlabels))

        plt.tight_layout()
        plt.colorbar()

        yield (res.pdf.name,),fig

def compare_violins(results, base_pdf = None):
    if not isinstance(results, dict):
        combined = aggregate_results(results)
    else:
        combined = results
    for obs in combined:
        figure = plt.figure()
        norms = None
        handles = []
        plt.title(str(obs),  y=1.05)
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

        #Small style
        dfs = plt.yticks()[0] - 1
        l = len(dfs) // 2  + 1 - ((len(dfs) // 2) % 2)
        mdiff = np.max(dfs)

        plt.yticks(np.linspace(-mdiff, mdiff, l) + 1)
        #

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
        plt.title(str(obs), y=1.05)
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
            #Without copy the /= operator affects the underlying data.
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
                #Without copy the /= operator affects the underlying data.
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
        #Small style
        dfs = plt.yticks()[0] - 1
        l = len(dfs) // 2  + 1 - ((len(dfs) // 2) % 2)
        mdiff = np.max(dfs)

        plt.yticks(np.linspace(-mdiff, mdiff, l) + 1)

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
    figure, ax = plt.subplots()

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

    ax.yaxis.set_major_locator(MaxNLocator(nbins=10, prune="both"))
    ax.xaxis.set_major_locator(MaxNLocator(nbins=8, prune="both"))

    plt.legend()
    yield (obs, b), figure


def plot_correlations(results):

    for result in results:
        pdf = result.pdf
        obs = result.obs

        Qs = iter(obs.meanQ)
        xgrid = pdf.make_xgrid()

        fl = pdf.make_flavors()

        figure, axarr = plt.subplots(len(fl), sharex=True,
                                     sharey=True,
                                     figsize=(8, len(fl)+3))

        for b in result.binlabels:
            Q = next(Qs)
            X = get_X(pdf, Q=Q, xgrid=xgrid, fl=fl, reshape=True)
            values, threshold = bin_corrs_from_X(result._all_vals.ix[b], X)
            ind = 0
            for f, axis in zip(fl, axarr):
                step = len(xgrid)
                current_vals = values[ind:ind+step]
                ind+=step
                line, = axis.plot(xgrid, current_vals)
                stacked = np.array([xgrid, current_vals]).T
                sel_ranges = split_ranges(stacked, abs(current_vals)>threshold,
                             filter_falses=True)
                for arr in sel_ranges:
                    x,y = arr.T
                    axis.plot(x,y, linewidth=3, color=line.get_color())
                    axis.axvspan(np.min(x), np.max(x), color="#eeeeff")
                axis.set_ylim([-1,1])
                axis.set_xscale('log')
                axis.set_ylabel("$%s$"%PDG_PARTONS[f])

                axis.axhline(threshold, c='r', ls='--')
                axis.axhline(-threshold, c='r', ls='--')

        axarr[0].set_title(str(obs) + "\n")
        plt.xlabel("$x$")
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
